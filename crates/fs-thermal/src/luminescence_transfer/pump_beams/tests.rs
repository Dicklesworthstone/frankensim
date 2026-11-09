use super::*;
use super::super::pump_column::{SaturablePumpCell, SaturablePumpColumn};

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-290),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn model() -> SingletTripletDecay {
    SingletTripletDecay::try_new([4e-19, 3e-19], [0.0; 2], [0.0; 2], [0.0; 2]).unwrap()
}

fn path(indices: &[usize], depth: f64, energy: f64) -> SaturablePumpPath<1> {
    SaturablePumpPath {
        photon_energies_j: [energy],
        segments: indices.iter().map(|&cell_index| SaturablePumpSegment {
            cell_index,
            optical_depth_per_ground: [[depth, 0.0]],
            passive_optical_depth: [0.0],
        }).collect(),
    }
}

fn limits() -> PumpBeamsLimits {
    PumpBeamsLimits {
        column: PumpColumnLimits {
            relative_tolerance: 2e-5,
            maximum_cell_evaluations: 500_000,
            ..PumpColumnLimits::default()
        },
        ..PumpBeamsLimits::default()
    }
}

#[test]
fn g1_identical_path_split_preserves_one_shared_history_and_original_column() {
    let count = 1e12;
    let authored = path(&[0, 1], 1.0 / count, 5e-19);
    let single = SaturablePumpBeams::try_new(
        &[model(); 2], &[authored.clone()], 2, 2, 4, || false,
    ).unwrap();
    let split = SaturablePumpBeams::try_new(
        &[model(); 2], &[authored.clone(), authored], 2, 2, 4, || false,
    ).unwrap();
    assert!(single.topological.is_some());
    assert!(split.topological.is_some());
    let initial = [[count, 0.0, 0.0]; 2];
    let one = single.advance(&initial, &[[5.0 * count]], 1.0, limits(), || false).unwrap();
    let two = split.advance(
        &initial, &[[1.25 * count], [3.75 * count]], 1.0, limits(), || false,
    ).unwrap();
    let column = SaturablePumpColumn::try_new(
        &[SaturablePumpCell {
            model: model(),
            optical_depth_per_ground: [[1.0 / count, 0.0]],
            passive_optical_depth: [0.0],
        }; 2],
        [5e-19], 2, || false,
    ).unwrap().advance(&initial, [5.0 * count], 1.0, limits().column, || false).unwrap();
    for cell in 0..2 {
        for state in 0..3 {
            close(two.cells[cell].populations[state], one.cells[cell].populations[state], 2e-10);
            close(one.cells[cell].populations[state], column.cells[cell].kinetics.populations[state], 4e-5);
        }
        close(two.cells[cell].pump_events[0], one.cells[cell].pump_events[0], 2e-10);
        close(two.cells[cell].pump_relaxation_heat_j[0], one.cells[cell].pump_relaxation_heat_j[0], 2e-10);
    }
    for (beam, fraction) in two.beams.iter().zip([0.25, 0.75]) {
        close(beam.transmitted_photons[0], fraction * one.beams[0].transmitted_photons[0], 2e-10);
        close(beam.absorbed_photons[0], fraction * one.beams[0].absorbed_photons[0], 2e-10);
    }
    assert_eq!(two.coupling_iterations, 0);
    assert_eq!(initial, [[count, 0.0, 0.0]; 2]);
}

#[test]
fn g1_overlapping_pumps_off_decay_each_original_cell_once_including_unvisited_cells() {
    let model = SingletTripletDecay::try_new(
        [4e-19, 3e-19], [2.0, 0.1], [0.3, 0.2], [0.4, 0.02],
    ).unwrap();
    let network = SaturablePumpBeams::try_new(
        &[model; 3],
        &[path(&[0, 1], 1e-12, 5e-19), path(&[1, 0], 1e-12, 5e-19)],
        3, 2, 4, || false,
    ).unwrap();
    assert!(network.topological.is_none());
    let initial = [
        [1e12, 2e11, 3e11],
        [7e11, 4e11, 1e11],
        [3e11, 1e11, 5e11],
    ];
    let duration = 4.0;
    let step = network.advance(&initial, &[[0.0]; 2], duration, limits(), || false).unwrap();
    let exact = SaturableSingletTripletPlan::prepare(
        model, [0.0; 2], [0.0; 2], duration, limits().column.kinetic, || false,
    ).unwrap();
    for (index, before) in initial.into_iter().enumerate() {
        let reference = exact.advance(before).unwrap();
        for state in 0..3 {
            close(step.cells[index].populations[state], reference.populations[state], 3e-12);
        }
        for target in 0..2 {
            close(step.cells[index].radiative_events[target], reference.radiative_events[target], 3e-12);
            close(step.cells[index].nonradiative_events[target], reference.nonradiative_events[target], 3e-12);
        }
    }
    assert_eq!(step.substeps, 1);
    assert_eq!(step.cell_evaluations, 3);
    assert_eq!(step.coupling_iterations, 0);
    for beam in &step.beams {
        assert_eq!(beam.incident_photons, [0.0]);
        assert_eq!(beam.transmitted_photons, [0.0]);
        assert_eq!(beam.absorbed_photons, [0.0]);
    }
}

#[test]
fn g1_passive_opposed_spectra_heat_shared_cells_once_and_empty_paths_transmit() {
    let paths: [SaturablePumpPath<2>; 2] = std::array::from_fn(|beam| {
        let indices = if beam == 0 { [0, 1] } else { [1, 0] };
        SaturablePumpPath {
            photon_energies_j: if beam == 0 { [5e-19, 4e-19] } else { [6e-19, 3e-19] },
            segments: indices.map(|cell_index| SaturablePumpSegment {
                cell_index,
                optical_depth_per_ground: [[0.0; 2]; 2],
                passive_optical_depth: if cell_index == 0 { [0.2, 0.3] } else { [0.7, 0.1] },
            }).into_iter().collect(),
        }
    });
    let incident = [[1e12, 2e12], [3e12, 4e12]];
    let duration = 0.3;
    let network = SaturablePumpBeams::try_new(&[model(); 2], &paths, 2, 2, 4, || false).unwrap();
    let step = network.advance(&[[1e12, 0.0, 0.0]; 2], &incident, duration, limits(), || false).unwrap();
    let mut expected_heat = 0.0;
    for beam in 0..2 {
        for lane in 0..2 {
            let input = incident[beam][lane] * duration;
            let depth: f64 = if lane == 0 { 0.9 } else { 0.4 };
            let output = input * (-depth).exp();
            close(step.beams[beam].transmitted_photons[lane], output, 2e-12);
            close(step.beams[beam].passive_absorbed_photons[lane], input - output, 2e-12);
            assert_eq!(step.beams[beam].absorbed_photons[lane], 0.0);
            expected_heat += (input - output) * paths[beam].photon_energies_j[lane];
        }
    }
    close(step.passive_heat_j.iter().sum(), expected_heat, 2e-12);
    assert_eq!(step.substeps, 1);
    assert_eq!(step.cell_evaluations, 2);
    assert_eq!(step.coupling_iterations, 0);
    let empty = SaturablePumpBeams::<1>::try_new(
        &[], &[path(&[], 0.0, 5e-19)], 1, 1, 1, || false,
    ).unwrap().advance(&[], &[[2e12]], 0.25, limits(), || false).unwrap();
    assert_eq!(empty.beams[0].incident_photons, [5e11]);
    assert_eq!(empty.beams[0].transmitted_photons, [5e11]);
    assert!(empty.cells.is_empty());
    assert_eq!(empty.cell_evaluations, 0);
    assert_eq!(empty.residual_j, 0.0);
    let no_paths = SaturablePumpBeams::<1>::try_new(
        &[], &[], 1, 1, 1, || false,
    ).unwrap().advance(&[], &[], 1.0, limits(), || false).unwrap();
    assert!(no_paths.beams.is_empty());
}

#[test]
fn g3_weak_energetic_beam_retains_heat_below_the_mean_photon_ulp() {
    let network = SaturablePumpBeams::try_new(
        &[model()],
        &[path(&[0], 1e-12, 4e-19), path(&[0], 1e-12, 8e-19)],
        1, 2, 2, || false,
    ).unwrap();
    let step = network.advance(
        &[[1e12, 0.0, 0.0]], &[[1e12], [1e-6]], 0.1, limits(), || false,
    ).unwrap();
    let weak = &step.beams[1];
    assert!(weak.absorbed_photons[0] > 0.0);
    assert!(step.cells[0].pump_relaxation_heat_j[0] > 0.0);
    close(
        step.cells[0].pump_relaxation_heat_j[0],
        weak.absorbed_photons[0] * 4e-19,
        2e-12,
    );
    close(
        weak.absorbed_photons[0],
        step.beams[0].absorbed_photons[0] * 1e-18,
        2e-12,
    );
    assert!(weak.photon_residual[0].abs() <= 2e-12 * weak.incident_photons[0]);
}

#[test]
fn g0_refuses_invalid_cells_bad_quanta_and_bounded_work_without_mutation() {
    let models = [model(); 2];
    assert!(SaturablePumpBeams::try_new(
        &models, &[path(&[0, 2], 1e-12, 5e-19)], 2, 1, 2, || false,
    ).is_err());
    // Returning to a cell is valid, but its second segment still consumes
    // the original aggregate geometry allowance.
    assert!(SaturablePumpBeams::try_new(
        &models, &[path(&[0, 1, 0], 1e-12, 5e-19)], 2, 1, 2, || false,
    ).is_err());
    assert!(SaturablePumpBeams::try_new(
        &models, &[path(&[0], 1e-12, 3e-19)], 2, 1, 2, || false,
    ).is_err());
    assert!(SaturablePumpBeams::try_new(
        &models, &[path(&[0, 1], 1e-12, 5e-19)], 2, 1, 1, || false,
    ).is_err());
    let network = SaturablePumpBeams::try_new(
        &models, &[path(&[0, 1], 1e-12, 5e-19), path(&[1, 0], 1e-12, 5e-19)],
        2, 2, 4, || false,
    ).unwrap();
    let initial = [[1e12, 0.0, 0.0]; 2];
    assert!(network.advance(&initial, &[[1e12]; 2], 1.0, limits(), || true).is_err());
    let mut bounded = limits();
    bounded.column.maximum_cell_evaluations = 1;
    assert!(network.advance(&initial, &[[1e12]; 2], 1.0, bounded, || false).is_err());
    bounded = limits();
    bounded.maximum_coupling_iterations = 1;
    bounded.column.maximum_substeps = 2;
    assert!(network.advance(&initial, &[[1e12]; 2], 1.0, bounded, || false).is_err());
    assert_eq!(initial, [[1e12, 0.0, 0.0]; 2]);
}
