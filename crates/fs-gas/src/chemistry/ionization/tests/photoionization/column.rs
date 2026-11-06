//! Analytic spatial shielding and temporal composition of photon fluence.
use super::*;
use crate::chemistry::photoionization::column::{
    PhotoionizationCell, PhotoionizationColumnControls, PhotoionizationColumnError,
    advance_monochromatic_column,
};

fn column_controls(maximum_cells: usize) -> PhotoionizationColumnControls {
    PhotoionizationColumnControls {
        cell: photo_controls(),
        maximum_cells,
        relative_photon_tolerance: 1e-12,
    }
}

#[test]
fn stratified_cells_match_independent_prefix_solution_and_close_photons_and_heat() {
    let (upstream, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let (downstream, _) = initial(0.5, 30000.0, 6000.0, 0.0, 0.0);
    let source = photon_source(atomic);
    let beam = exposure(1.0, 2.0, 2.5);
    let result = advance_monochromatic_column(
        &[
            PhotoionizationCell {
                initial: &upstream,
                source: &source,
            },
            PhotoionizationCell {
                initial: &downstream,
                source: &source,
            },
        ],
        beam,
        column_controls(2),
        || false,
    )
    .unwrap();
    // With one cross section, every prefix obeys the scalar conservation law
    // at its *total* optical depth, regardless of how its gas is partitioned.
    let first_depth = (2.0_f64.exp_m1() * (-2.5_f64).exp()).ln_1p();
    let total_depth = (3.0_f64.exp_m1() * (-2.5_f64).exp()).ln_1p();
    let second_depth = total_depth - first_depth;
    relative(result.cells()[0].optical_depths()[1], first_depth, 5e-14);
    relative(result.cells()[1].optical_depths()[1], second_depth, 5e-14);
    assert!(first_depth / 2.0 < second_depth);
    let expected_absorbed = (3.0 - total_depth) * beam.illuminated_area.value() / SIGMA;
    let expected_transmitted = (2.5 - 3.0 + total_depth) * beam.illuminated_area.value() / SIGMA;
    let b = result.balance();
    relative(b.absorbed_photons, expected_absorbed, 5e-14);
    relative(b.transmitted_photons, expected_transmitted, 5e-14);
    relative(
        b.absorbed_photons + b.transmitted_photons,
        b.incident_photons,
        5e-14,
    );
    relative(
        result.cells()[1].balance().incident_photons,
        result.cells()[0].balance().transmitted_photons,
        5e-14,
    );
    relative(
        b.chemical_storage_change_j,
        expected_absorbed * (CHI / NA),
        5e-14,
    );
    relative(
        b.electron_sensible_energy_change_j,
        expected_absorbed * (16.0 * EV - CHI / NA),
        5e-14,
    );
    let material_energy: f64 = result
        .cells()
        .iter()
        .map(|cell| {
            relative(
                cell.state().mass().value(),
                cell.initial_state().mass().value(),
                5e-14,
            );
            assert_eq!(cell.state().amounts()[1], cell.state().amounts()[2]);
            assert_eq!(
                cell.state().heavy_temperature(),
                cell.initial_state().heavy_temperature()
            );
            cell.state().internal_energy().value() - cell.initial_state().internal_energy().value()
        })
        .sum();
    assert!((material_energy - b.absorbed_energy_j).abs() < 1e-7);
    assert!(b.energy_residual_j.abs() < 1e-7);
    assert!(result.photon_balance_residual().abs() < 1e-12 * b.incident_photons);
    assert!(
        result.cells()[1].state().electron_temperature().value()
            < downstream.electron_temperature().value()
    );
}

#[test]
fn different_cell_cross_sections_match_high_precision_fluence_oracle() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = photon_source(atomic.clone());
    let stronger = photon_source_with_cross_section(atomic, 4.0 * SIGMA);
    let result = advance_monochromatic_column(
        &[
            PhotoionizationCell {
                initial: &initial,
                source: &source,
            },
            PhotoionizationCell {
                initial: &initial,
                source: &stronger,
            },
        ],
        exposure(1.0, 1.0, 3.0),
        column_controls(2),
        || false,
    )
    .unwrap();
    // Independent 80-digit Decimal exp/log evaluation at tau0=[1,4], dose=3.
    // Second dose = 4*(3-1+ln(1+(e-1)*exp(-3))).
    relative(
        result.cells()[0].state().amounts()[0].value(),
        0.08208512654808950,
        5e-14,
    );
    relative(
        result.cells()[1].state().amounts()[0].value(),
        0.0032161844389138644,
        5e-14,
    );
    relative(
        result.balance().transmitted_photons / NA,
        4.3412052439480135 / 4.0,
        5e-14,
    );
    assert_eq!(result.cells()[0].source(), &source);
    assert_eq!(result.cells()[1].source(), &stronger);
}

#[test]
fn layered_exposure_composes_in_time_despite_changing_upstream_opacity() {
    let (upstream, atomic) = initial(0.1, 30000.0, 4000.0, 0.0, 0.0);
    let (downstream, _) = initial(0.0, 12000.0, 6000.0, 0.0, 0.0);
    let source = photon_source(atomic.clone());
    let other = photon_source_with_cross_section(atomic, 4.0 * SIGMA);
    let cells = [
        PhotoionizationCell {
            initial: &upstream,
            source: &source,
        },
        PhotoionizationCell {
            initial: &downstream,
            source: &other,
        },
    ];
    let full_exposure = exposure(0.9, 2.0, 4.0);
    let full =
        advance_monochromatic_column(&cells, full_exposure, column_controls(2), || false).unwrap();
    let mut partial = full_exposure;
    partial.duration = Time::new(0.25);
    let first =
        advance_monochromatic_column(&cells, partial, column_controls(2), || false).unwrap();
    let next_cells: Vec<_> = first
        .cells()
        .iter()
        .map(|cell| PhotoionizationCell {
            initial: cell.state(),
            source: cell.source(),
        })
        .collect();
    partial.duration = Time::new(0.75);
    let second =
        advance_monochromatic_column(&next_cells, partial, column_controls(2), || false).unwrap();
    for (whole, split) in full.cells().iter().zip(second.cells()) {
        for (a, b) in whole.state().amounts().iter().zip(split.state().amounts()) {
            relative(a.value(), b.value(), 2e-13);
        }
        relative(
            whole.state().electron_temperature().value(),
            split.state().electron_temperature().value(),
            2e-13,
        );
    }
    relative(
        full.balance().absorbed_photons,
        first.balance().absorbed_photons + second.balance().absorbed_photons,
        2e-13,
    );
    relative(
        full.balance().transmitted_photons,
        first.balance().transmitted_photons + second.balance().transmitted_photons,
        2e-13,
    );
    relative(
        full.balance().absorbed_energy_j,
        first.balance().absorbed_energy_j + second.balance().absorbed_energy_j,
        2e-13,
    );
    // Ionization opens the first cell substantially: the later downstream mean
    // flux must increase, rather than using the same initial Beer-Lambert factor.
    assert!(
        second.cells()[1].exposure().beam.number_flux_per_m2_s
            > first.cells()[1].exposure().beam.number_flux_per_m2_s
    );
}

#[test]
fn small_transmitted_population_still_ionizes_and_heats_a_shielded_cell() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = photon_source(atomic.clone());
    let weaker = photon_source_with_cross_section(atomic, SIGMA / 40.0);
    let beam = exposure(1.0, 40.0, 1.0);
    let result = advance_monochromatic_column(
        &[
            PhotoionizationCell {
                initial: &initial,
                source: &source,
            },
            PhotoionizationCell {
                initial: &initial,
                source: &weaker,
            },
        ],
        beam,
        column_controls(2),
        || false,
    )
    .unwrap();
    // Upstream transmission is below one ulp of its incident inventory, yet
    // about 110,000 photons reach the second cell in this macroscopic column.
    let entering = (NA / 40.0) * (-40.0_f64).exp() * 1.0_f64.exp_m1();
    let absorbed = entering * -(-1.0_f64).exp_m1();
    let second = &result.cells()[1];
    relative(second.balance().incident_photons, entering, 1e-13);
    relative(second.balance().absorbed_photons, absorbed, 1e-13);
    relative(
        second.balance().transmitted_photons,
        entering * (-1.0_f64).exp(),
        1e-13,
    );
    relative(second.reaction_extent().value(), absorbed / NA, 1e-13);
    assert!(second.state().amounts()[2].value() > 0.0);
    relative(
        second.state().electron_temperature().value(),
        (16.0 * EV * NA - CHI) / (1.5 * R),
        5e-14,
    );
    relative(
        second.balance().absorbed_energy_j,
        absorbed * 16.0 * EV,
        1e-13,
    );
}

#[test]
fn empty_dark_budgeted_and_late_refusals_preserve_caller_states() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let before = initial.clone();
    let source = photon_source(atomic);
    let beam = exposure(1.0, 1.0, 1.0);
    let cell = PhotoionizationCell {
        initial: &initial,
        source: &source,
    };
    let clear = advance_monochromatic_column(&[], beam, column_controls(0), || false).unwrap();
    assert!(clear.cells().is_empty());
    assert_eq!(
        clear.balance().incident_photons,
        clear.balance().transmitted_photons
    );
    assert_eq!(clear.balance().absorbed_energy_j, 0.0);
    assert_eq!(
        advance_monochromatic_column(&[cell; 2], beam, column_controls(1), || false),
        Err(PhotoionizationColumnError::CellLimit {
            requested: 2,
            maximum: 1
        })
    );
    let mut dark = beam;
    dark.beam.number_flux_per_m2_s = 0.0;
    let dark_result =
        advance_monochromatic_column(&[cell; 2], dark, column_controls(2), || false).unwrap();
    for result in dark_result.cells() {
        assert_eq!(result.state(), &initial);
    }
    assert_eq!(dark_result.balance().incident_photons, 0.0);
    let (_, mismatched) = initial_with_shift();
    let mismatched = photon_source(mismatched);
    // Even an exactly dark downstream cell must retain the original source
    // binding. A valid prefix is never published when a later cell refuses.
    assert!(matches!(
        advance_monochromatic_column(
            &[
                cell,
                PhotoionizationCell {
                    initial: &initial,
                    source: &mismatched
                },
            ],
            exposure(1.0, 1e6, 1.0),
            column_controls(2),
            || false
        ),
        Err(PhotoionizationColumnError::Cell {
            index: 1,
            error: PhotoionizationError::Invalid(_)
        })
    ));
    let mut cell_polls = 0;
    advance_monochromatic_photoionization(&initial, &source, beam, photo_controls(), || {
        cell_polls += 1;
        false
    })
    .unwrap();
    let mut polls = 0;
    assert_eq!(
        advance_monochromatic_column(&[cell; 2], beam, column_controls(2), || {
            polls += 1;
            polls > cell_polls + 1
        }),
        Err(PhotoionizationColumnError::Cancelled)
    );
    assert_eq!(initial, before);
}
