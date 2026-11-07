//! Spatial shielding, spectral hardening and conservative gas heating.
use super::*;
use crate::chemistry::photoionization::column::{
    PhotoionizationCell, PhotoionizationColumnControls, advance_monochromatic_column,
};
use crate::chemistry::photoionization::spectrum::column::{
    SpectralPhotoionizationColumnError as ColumnError, SpectralPhotoionizationColumnResult,
    SpectralPhotoionizationColumnSolve, advance_spectral_column,
};

fn column_controls() -> SpectralPhotoionizationColumnSolve {
    SpectralPhotoionizationColumnSolve {
        cell: solve_controls(),
        maximum_cells: 16,
        maximum_cell_groups: 64,
    }
}
fn column_run(
    cells: &[PhotoionizationCell<'_>],
    beams: &[PhotonNumberFluxV1],
    controls: SpectralPhotoionizationColumnSolve,
) -> SpectralPhotoionizationColumnResult {
    advance_spectral_column(cells, beams, area(), Time::new(1.0), controls, || false).unwrap()
}

// Independent RK4 of the simultaneous Beer-Lambert column ODE, not split
// monochromatic updates. State is [n0,n1,abs00,abs01,abs10,abs11] in mol.
fn column_reference() -> [f64; 6] {
    let rate = |state: [f64; 6]| {
        let mut result = [0.0; 6];
        for (g, (incoming, opacity)) in [(0.5, 3.0), (4.0 / 3.0, 0.6)].into_iter().enumerate() {
            let mut photons = incoming;
            for c in 0..2 {
                let tau = opacity * state[c];
                let absorbed = photons * -(-tau).exp_m1();
                result[c] -= absorbed;
                result[2 + 2 * c + g] = absorbed;
                photons *= (-tau).exp();
            }
        }
        result
    };
    let add = |a: [f64; 6], b: [f64; 6], h: f64| std::array::from_fn(|i| a[i] + h * b[i]);
    let mut state = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let dt = 1.0 / 8192.0;
    for _ in 0..8192 {
        let k1 = rate(state);
        let k2 = rate(add(state, k1, dt * 0.5));
        let k3 = rate(add(state, k2, dt * 0.5));
        let k4 = rate(add(state, k3, dt));
        state = std::array::from_fn(|i| {
            state[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    state
}

#[test]
fn one_band_matches_exact_heterogeneous_column_in_thin_thick_and_depleted_limits() {
    let (front, atomic) = initial(0.1, 18000.0, 4000.0, 0.0, 0.0);
    let (back, shifted) = initial(0.0, 12000.0, 6000.0, 1e12, 1e8);
    let source = photon_source(atomic);
    let other = photon_source_with_cross_section(shifted, 0.3 * SIGMA);
    let cells = [
        PhotoionizationCell {
            initial: &front,
            source: &source,
        },
        PhotoionizationCell {
            initial: &back,
            source: &other,
        },
    ];
    for (tau, dose) in [(1.0, 2.0), (1e-12, 50.0), (1e6, 100.0)] {
        let beam = exposure(0.9, tau, dose);
        let exact = advance_monochromatic_column(
            &cells,
            beam,
            PhotoionizationColumnControls {
                cell: photo_controls(),
                maximum_cells: 2,
                relative_photon_tolerance: 1e-12,
            },
            || false,
        )
        .unwrap();
        let mut controls = column_controls();
        controls.cell.initial_step_s = 1.0;
        let result = advance_spectral_column(
            &cells,
            &[beam.beam],
            beam.illuminated_area,
            beam.duration,
            controls,
            || false,
        )
        .unwrap();
        assert_eq!(result.accepted_steps(), 1);
        assert_eq!(result.rejected_steps(), 0);
        for (actual, expected) in result.cells().iter().zip(exact.cells()) {
            for (a, b) in actual
                .state()
                .amounts()
                .iter()
                .zip(expected.state().amounts())
            {
                relative(a.value(), b.value(), 2e-11);
            }
            relative(
                actual.state().electron_temperature().value(),
                expected.state().electron_temperature().value(),
                2e-13,
            );
            relative(
                actual.balance().absorbed_photons,
                expected.balance().absorbed_photons,
                2e-13,
            );
            relative(
                actual.balance().transmitted_photons,
                expected.balance().transmitted_photons,
                2e-11,
            );
            assert_eq!(actual.source(), expected.source());
            assert_eq!(
                actual.state().thermodynamics(),
                actual.initial_state().thermodynamics()
            );
        }
        relative(
            result.balance().absorbed_energy_j,
            exact.balance().absorbed_energy_j,
            2e-13,
        );
    }
}

#[test]
fn one_cell_matches_the_original_spectral_integrator() {
    let (initial, atomic) = initial(0.1, 18000.0, 4000.0, 0.0, 0.0);
    let source = spectral_source(atomic);
    let mut controls = column_controls();
    controls.cell.relative_tolerance = 0.1;
    controls.cell.initial_step_s = 0.125;
    controls.cell.maximum_step_s = 0.125;
    let expected = spectral_run(&initial, &source, &groups(), controls.cell);
    let result = column_run(
        &[PhotoionizationCell {
            initial: &initial,
            source: &source,
        }],
        &groups(),
        controls,
    );
    let actual = &result.cells()[0];
    for (a, b) in actual
        .state()
        .amounts()
        .iter()
        .zip(expected.state().amounts())
    {
        relative(a.value(), b.value(), 5e-14);
    }
    relative(
        actual.state().electron_temperature().value(),
        expected.state().electron_temperature().value(),
        5e-14,
    );
    for (a, b) in actual
        .group_balances()
        .iter()
        .zip(expected.group_balances())
    {
        relative(a.absorbed_photons, b.absorbed_photons, 5e-14);
        relative(a.transmitted_photons, b.transmitted_photons, 5e-14);
    }
    assert_eq!(result.accepted_steps(), expected.accepted_steps());
}

#[test]
fn evolving_shielding_and_downstream_hardening_match_independent_coupled_ode() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = spectral_source(atomic);
    let cell = PhotoionizationCell {
        initial: &initial,
        source: &source,
    };
    let beams = groups();
    let result = column_run(&[cell; 2], &beams, column_controls());
    let reference = column_reference();
    let mut material_energy = 0.0;
    for (c, cell) in result.cells().iter().enumerate() {
        relative(cell.state().amounts()[0].value(), reference[c], 2e-6);
        let absorbed: f64 = reference[2 + 2 * c..4 + 2 * c].iter().sum();
        let deposited: f64 = (0..2)
            .map(|g| reference[2 + 2 * c + g] * (beams[g].photon_energy.value() * NA - CHI))
            .sum();
        relative(cell.reaction_extent().value(), absorbed, 2e-6);
        relative(
            cell.state().electron_temperature().value(),
            deposited / (1.5 * R * absorbed),
            2e-6,
        );
        for (g, b) in cell.group_balances().iter().enumerate() {
            relative(b.absorbed_photons / NA, reference[2 + 2 * c + g], 2e-6);
            relative(
                b.absorbed_photons + b.transmitted_photons,
                b.incident_photons,
                2e-14,
            );
            relative(
                b.absorbed_energy_j,
                b.absorbed_photons * beams[g].photon_energy.value(),
                2e-14,
            );
            if c > 0 {
                assert_eq!(
                    b.incident_photons,
                    result.cells()[c - 1].group_balances()[g].transmitted_photons
                );
            }
        }
        assert_eq!(cell.state().amounts()[1], cell.state().amounts()[2]);
        relative(cell.state().mass().value(), initial.mass().value(), 2e-13);
        assert_eq!(
            cell.state().heavy_temperature(),
            initial.heavy_temperature()
        );
        assert_eq!(cell.state().thermodynamics(), initial.thermodynamics());
        material_energy +=
            cell.state().internal_energy().value() - initial.internal_energy().value();
    }
    for (g, b) in result.group_balances().iter().enumerate() {
        relative(
            b.absorbed_photons / NA,
            reference[2 + g] + reference[4 + g],
            2e-6,
        );
        relative(
            b.absorbed_photons + b.transmitted_photons,
            b.incident_photons,
            2e-14,
        );
    }
    assert!((material_energy - result.balance().absorbed_energy_j).abs() < 1e-7);
    assert!(result.balance().energy_residual_j.abs() < 1e-7);
    assert!(result.maximum_normalized_local_error() <= 1.0);
    // The front becomes transparent during the exposure. Its original opacity
    // would deliver less than a third of the actual downstream soft photons.
    let frozen_soft = beams[0].number_flux_per_m2_s * area().value() * (-3.0_f64).exp();
    let back = &result.cells()[1];
    assert!(back.group_balances()[0].incident_photons > 3.0 * frozen_soft);
    let hard_ratio =
        back.group_balances()[1].incident_photons / back.group_balances()[0].incident_photons;
    assert!(hard_ratio > beams[1].number_flux_per_m2_s / beams[0].number_flux_per_m2_s);
    assert!(
        back.state().electron_temperature().value()
            > result.cells()[0].state().electron_temperature().value() + 5000.0
    );
}

#[test]
fn spatial_band_competition_converges_at_second_order_and_loses_order_bias() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = spectral_source(atomic);
    let cells = [PhotoionizationCell {
        initial: &initial,
        source: &source,
    }; 2];
    let forward = groups();
    let reverse = [forward[1], forward[0]];
    let reference = column_reference();
    let mut errors = Vec::new();
    let mut biases = Vec::new();
    for step in [0.25, 0.125, 0.0625] {
        let mut controls = column_controls();
        controls.cell.relative_tolerance = 0.1;
        controls.cell.initial_step_s = step;
        controls.cell.maximum_step_s = step;
        let a = column_run(&cells, &forward, controls);
        let b = column_run(&cells, &reverse, controls);
        assert_eq!(a.rejected_steps(), 0);
        assert_eq!(b.rejected_steps(), 0);
        errors.push(
            (0..2)
                .map(|c| (a.cells()[c].state().amounts()[0].value() - reference[c]).abs())
                .sum::<f64>(),
        );
        biases.push(
            (0..2)
                .map(|c| {
                    (a.cells()[c].state().amounts()[0].value()
                        - b.cells()[c].state().amounts()[0].value())
                    .abs()
                })
                .sum::<f64>(),
        );
    }
    for values in [&errors, &biases] {
        for pair in values.windows(2) {
            assert!(
                pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2,
                "{values:?}"
            );
        }
    }
}

#[test]
fn each_band_retains_tiny_transmission_that_still_creates_downstream_electrons() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = photon_source(atomic.clone());
    let weak = photon_source_with_cross_section(atomic, SIGMA / 40.0);
    let cells = [
        PhotoionizationCell {
            initial: &initial,
            source: &source,
        },
        PhotoionizationCell {
            initial: &initial,
            source: &weak,
        },
    ];
    let beams = [16.0, 20.0].map(|ev| PhotonNumberFluxV1 {
        photon_energy: Energy::new(ev * EV),
        number_flux_per_m2_s: 1.0 / SIGMA,
    });
    let mut controls = column_controls();
    controls.cell.relative_tolerance = 1e-9;
    let result = advance_spectral_column(
        &cells,
        &beams,
        Area::new(SIGMA * NA / 40.0),
        Time::new(1.0),
        controls,
        || false,
    )
    .unwrap();
    // Equal cross sections and fluxes absorb equal band populations in the
    // simultaneous ODE. Each escaped count is below one ulp of its input.
    let entering = 0.5 * (NA / 40.0) * (-40.0_f64).exp() * 2.0_f64.exp_m1();
    let absorbed = entering * -(-1.0_f64).exp_m1();
    let back = &result.cells()[1];
    for g in 0..2 {
        let b = back.group_balances()[g];
        relative(b.incident_photons, entering, 2e-6);
        relative(b.absorbed_photons, absorbed, 2e-6);
        relative(b.transmitted_photons, entering * (-1.0_f64).exp(), 2e-6);
        assert_eq!(
            b.incident_photons,
            result.cells()[0].group_balances()[g].transmitted_photons
        );
        assert!(b.incident_photons < f64::EPSILON * result.group_balances()[g].incident_photons);
    }
    relative(back.reaction_extent().value(), 2.0 * absorbed / NA, 2e-6);
    assert!(back.state().amounts()[2].value() > 0.0);
    relative(
        back.state().electron_temperature().value(),
        (18.0 * EV * NA - CHI) / (1.5 * R),
        2e-6,
    );
}

#[test]
fn combined_birth_temperatures_and_original_references_survive_spatial_transport() {
    let beams = [14.0, 22.0].map(|ev| PhotonNumberFluxV1 {
        photon_energy: Energy::new(ev * EV),
        number_flux_per_m2_s: 1.0 / SIGMA,
    });
    let mut results = Vec::new();
    for (elemental, charge) in [(0.0, 0.0), (1e12, 1e8)] {
        let (initial, atomic) = initial(0.0, 12000.0, 4000.0, elemental, charge);
        let source = photon_source(atomic);
        let mut controls = column_controls();
        controls.cell.minimum_electron_temperature_k = 5000.0;
        controls.cell.maximum_electron_temperature_k = 45000.0;
        controls.cell.relative_tolerance = 0.1;
        controls.cell.initial_step_s = 0.125;
        controls.cell.maximum_step_s = 0.125;
        // Individual subflows have birth temperatures outside this interval;
        // the physical combined pool is inside it in both cells throughout.
        let result = column_run(
            &[PhotoionizationCell {
                initial: &initial,
                source: &source,
            }; 2],
            &beams,
            controls,
        );
        for cell in result.cells() {
            assert!((5000.0..=45000.0).contains(&cell.state().electron_temperature().value()));
            assert_eq!(cell.state().thermodynamics(), initial.thermodynamics());
        }
        assert!(result.balance().energy_residual_j.abs() < controls.cell.energy_tolerance_j);
        results.push(result);
    }
    for (a, b) in results[0].cells().iter().zip(results[1].cells()) {
        assert_eq!(a.state().amounts(), b.state().amounts());
        relative(
            a.state().electron_temperature().value(),
            b.state().electron_temperature().value(),
            1e-9,
        );
    }
}

#[test]
fn empty_dark_subthreshold_and_zero_duration_columns_preserve_states_and_counts() {
    let (initial, atomic) = initial(0.1, 12000.0, 4000.0, 0.0, 0.0);
    let source = spectral_source(atomic);
    let cells = [PhotoionizationCell {
        initial: &initial,
        source: &source,
    }; 2];
    let dark = PhotonNumberFluxV1 {
        photon_energy: Energy::new(16.0 * EV),
        number_flux_per_m2_s: 0.0,
    };
    let below = PhotonNumberFluxV1 {
        photon_energy: Energy::new(10.0 * EV),
        number_flux_per_m2_s: 1e20,
    };
    for beams in [vec![], vec![dark], vec![dark, below]] {
        let result = column_run(&cells, &beams, column_controls());
        assert_eq!(result.accepted_steps(), 0);
        for cell in result.cells() {
            assert_eq!(cell.state(), &initial);
            assert_eq!(cell.reaction_extent().value(), 0.0);
            for b in cell.group_balances() {
                assert_eq!(b.absorbed_photons, 0.0);
                assert_eq!(b.incident_photons, b.transmitted_photons);
            }
        }
    }
    let mut controls = column_controls();
    controls.maximum_cells = 0;
    controls.maximum_cell_groups = 0;
    let clear = column_run(&[], &groups(), controls);
    assert!(clear.cells().is_empty());
    for b in clear.group_balances() {
        assert_eq!(b.incident_photons, b.transmitted_photons);
        assert_eq!(b.absorbed_energy_j, 0.0);
    }
    let zero = advance_spectral_column(
        &cells,
        &groups(),
        area(),
        Time::new(0.0),
        column_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(zero.accepted_steps(), 0);
    assert_eq!(zero.balance().incident_photons, 0.0);
    for cell in zero.cells() {
        assert_eq!(cell.state(), &initial);
    }
}

#[test]
fn work_limits_source_refusals_and_late_cancellation_publish_no_partial_column() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let before = initial.clone();
    let source = spectral_source(atomic.clone());
    let cell = PhotoionizationCell {
        initial: &initial,
        source: &source,
    };
    let mut controls = column_controls();
    controls.maximum_cells = 1;
    let mut polls = 0;
    assert!(matches!(
        advance_spectral_column(
            &[cell; 2],
            &groups(),
            area(),
            Time::new(1.0),
            controls,
            || {
                polls += 1;
                false
            }
        ),
        Err(ColumnError::CellLimit { .. })
    ));
    assert_eq!(polls, 1);
    controls = column_controls();
    controls.maximum_cell_groups = 3;
    polls = 0;
    assert_eq!(
        advance_spectral_column(
            &[cell; 2],
            &groups(),
            area(),
            Time::new(1.0),
            controls,
            || {
                polls += 1;
                false
            }
        ),
        Err(ColumnError::CellGroupLimit {
            cells: 2,
            groups: 2,
            maximum: 3
        })
    );
    assert_eq!(polls, 1);
    controls = column_controls();
    controls.cell.maximum_steps = 1;
    controls.cell.initial_step_s = 1e-4;
    assert_eq!(
        advance_spectral_column(
            &[cell; 2],
            &groups(),
            area(),
            Time::new(1.0),
            controls,
            || false
        ),
        Err(ColumnError::StepBudget)
    );
    let broad = photon_source(atomic);
    let invalid_dark = [PhotonNumberFluxV1 {
        photon_energy: Energy::new(21.0 * EV),
        number_flux_per_m2_s: 0.0,
    }];
    assert!(matches!(
        advance_spectral_column(
            &[
                PhotoionizationCell {
                    initial: &initial,
                    source: &broad
                },
                cell
            ],
            &invalid_dark,
            area(),
            Time::new(1.0),
            column_controls(),
            || false
        ),
        Err(ColumnError::Cell {
            index: 1,
            error: SpectralPhotoionizationError::Group { index: 0, .. }
        })
    ));
    let mut completed_polls = 0;
    advance_spectral_column(
        &[cell; 2],
        &groups(),
        area(),
        Time::new(1.0),
        column_controls(),
        || {
            completed_polls += 1;
            false
        },
    )
    .unwrap();
    polls = 0;
    assert_eq!(
        advance_spectral_column(
            &[cell; 2],
            &groups(),
            area(),
            Time::new(1.0),
            column_controls(),
            || {
                polls += 1;
                polls >= completed_polls - 1
            }
        ),
        Err(ColumnError::Cancelled)
    );
    assert_eq!(initial, before);
}
