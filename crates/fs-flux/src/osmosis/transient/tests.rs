use super::*;
use crate::osmosis::OsmoticPressureVolume;
use fs_material::osmosis::IdealDiluteSolution;

#[path = "interface_tests.rs"]
mod interfaces;
#[path = "kinetic_tests.rs"]
mod kinetic;
#[path = "sorption_tests.rs"]
mod sorption;

fn bulk(concentrations: &[f64], volume: f64, pressure: f64) -> OsmoticCompartment {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(
            298.15,
            18e-6,
            volume,
            concentrations.iter().map(|c| c * volume).collect(),
        )
        .unwrap(),
        OsmoticPressureVolume::new(volume, pressure, 0.0).unwrap(),
        volume * 1e-6,
    )
    .unwrap()
}
fn membrane(lp: f64, sigma: f64, permeability: f64, count: usize) -> MembraneTransport {
    MembraneTransport::new(
        HydraulicPermeance::new(lp).unwrap(),
        vec![
            (
                ReflectionCoefficient::new(sigma).unwrap(),
                SolutePermeance::new(permeability).unwrap()
            );
            count
        ],
    )
    .unwrap()
}
fn film(cells: usize, concentration: f64) -> TransientFilm {
    TransientFilm::new(1e-4, cells, vec![1e-9], vec![concentration; cells]).unwrap()
}
fn retaining_pair(cells: usize, volume: f64, right_pressure: f64) -> TransientOsmoticPair {
    TransientOsmoticPair::new(
        [
            bulk(&[50.0], volume, 101_325.0),
            bulk(&[100.0], volume, right_pressure),
        ],
        1e-4,
        membrane(2e-11, 1.0, 0.0, 1),
        Some(film(cells, 50.0)),
        Some(film(cells, 100.0)),
        1 << 24,
    )
    .unwrap()
}
fn close(left: f64, right: f64, relative: f64, absolute: f64) {
    assert!(
        (left - right).abs() <= absolute + relative * left.abs().max(right.abs()),
        "{left:.16e} != {right:.16e}"
    );
}
fn snapshot(pair: &TransientOsmoticPair) -> Vec<u64> {
    let mut values = vec![pair.time_s().to_bits()];
    for b in pair.bulks() {
        values.push(b.solution().solvent_volume_m3().to_bits());
        values.extend(b.solution().solute_moles().iter().map(|n| n.to_bits()));
    }
    for side in 0..2 {
        for cell in 0..pair.film_cell_count(side).unwrap() {
            values.extend(
                pair.film_solute_moles(side, cell)
                    .unwrap()
                    .iter()
                    .map(|n| n.to_bits()),
            );
        }
    }
    values
}

fn diffusing_reservoirs() -> TransientOsmoticPair {
    TransientOsmoticPair::new(
        [bulk(&[100.0], 1e-6, 0.0), bulk(&[0.0], 1e-6, 0.0)],
        1e-4,
        membrane(0.0, 0.0, 1e-3, 1),
        None,
        None,
        1 << 20,
    )
    .unwrap()
}

fn layered_startup() -> TransientOsmoticPair {
    TransientOsmoticPair::new(
        [
            bulk(&[50.0], 1e-6, 101_325.0),
            bulk(&[100.0], 1e-6, 101_325.0),
        ],
        1e-4,
        membrane(2e-11, 1.0, 0.0, 1),
        Some(
            TransientFilm::from_cells(
                vec![80e-6, 15e-6, 5e-6],
                1,
                vec![1e-9, 3e-10, 5e-10],
                vec![50.0; 3],
            )
            .unwrap(),
        ),
        Some(
            TransientFilm::from_cells(
                vec![5e-6, 15e-6, 80e-6],
                1,
                vec![5e-10, 3e-10, 1e-9],
                vec![100.0; 3],
            )
            .unwrap(),
        ),
        1 << 20,
    )
    .unwrap()
}

#[test]
fn adaptive_reservoir_diffusion_converges_to_the_analytical_exponential() {
    // Equal, fixed-volume reservoirs: c_L=50+50*exp(-2*A*P*t/V).
    let exact = 50.0 + 50.0 * (-2.0_f64).exp();
    let mut errors = Vec::new();
    let mut accepted = Vec::new();
    for relative_tolerance in [1e-2, 1e-5] {
        let mut pair = diffusing_reservoirs();
        let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
        let report = pair
            .advance(
                10.0,
                TransientAdvanceControls {
                    maximum_step_s: 10.0,
                    relative_tolerance,
                    ..Default::default()
                },
                &mut workspace,
                || false,
            )
            .unwrap();
        let observed = pair.bulks()[0].solution().solute_moles()[0] / 1e-6;
        errors.push((observed - exact).abs());
        accepted.push(report.accepted_steps);
        assert!(report.rejected_steps > 0);
        assert!(report.maximum_scaled_error <= 1.0);
        assert_eq!(
            report.step_evaluations,
            3 * (report.accepted_steps + report.rejected_steps)
        );
        assert_eq!(report.root_evaluations, report.step_evaluations);
        assert_eq!(pair.time_s().to_bits(), 10.0_f64.to_bits());
        close(pair.total_amount(0), 1e-4, 1e-12, 1e-18);
        assert_eq!(report.solvent_transfer_m3, 0.0);
        assert_eq!(report.water_balance_residual_m3, 0.0);
    }
    assert!(accepted[1] > accepted[0]);
    assert!(
        errors[1] < 0.15 * errors[0],
        "refinement errors: {errors:?}"
    );
    assert!(errors[1] < 0.06, "tight error: {}", errors[1]);
}

#[test]
fn adaptive_membrane_faces_resolve_layered_film_startup_when_amount_scales_hide_it() {
    let mut reference = layered_startup();
    let mut reference_endpoint = None;
    for _ in 0..2048 {
        reference_endpoint = Some(
            reference
                .step(1.0 / 2048.0, TransientStepControls::default(), || false)
                .unwrap(),
        );
    }
    let reference_endpoint = reference_endpoint.unwrap();
    let mut errors = Vec::new();
    let mut counts = Vec::new();
    for face_tolerance in [1e6, 1e-5] {
        let mut pair = layered_startup();
        let before: Vec<f64> = (0..pair.rows).map(|row| pair.amount(row, 0)).collect();
        let initial_water = pair.total_solvent_volume_m3();
        let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
        let report = pair
            .advance(
                1.0,
                TransientAdvanceControls {
                    maximum_step_s: 1.0,
                    // These large absolute scales deliberately disable amount/volume
                    // control: surface concentrations must independently see startup.
                    absolute_amount_tolerance_mol: 1e3,
                    absolute_water_tolerance_m3: 1.0,
                    absolute_face_concentration_tolerance_mol_m3: face_tolerance,
                    ..Default::default()
                },
                &mut workspace,
                || false,
            )
            .unwrap();
        let error = (report.endpoint.membrane_left_concentration_mol_m3[0]
            - reference_endpoint.membrane_left_concentration_mol_m3[0])
            .abs()
            .max(
                (report.endpoint.membrane_right_concentration_mol_m3[0]
                    - reference_endpoint.membrane_right_concentration_mol_m3[0])
                    .abs(),
            );
        errors.push(error);
        counts.push(report.accepted_steps);
        let direct_residual: f64 = before
            .iter()
            .enumerate()
            .map(|(row, amount)| pair.amount(row, 0) - amount)
            .sum();
        close(
            report.solute_balance_residual_mol[0],
            direct_residual,
            0.0,
            1e-25,
        );
        close(pair.total_solvent_volume_m3(), initial_water, 1e-13, 1e-20);
        assert!(report.solute_balance_residual_mol[0].abs() < 1e-16);
        assert!(report.solvent_transfer_m3 > 0.0);
        assert!(report.maximum_scaled_error <= 1.0);
        for row in 0..pair.rows {
            assert!(pair.amount(row, 0) >= 0.0);
        }
    }
    assert_eq!(counts[0], 1);
    assert!(counts[1] > 1);
    assert!(
        errors[1] < errors[0] * 0.2,
        "face refinement errors: {errors:?}"
    );
}

#[test]
fn adaptive_late_cancellation_restores_the_entire_window_and_reuses_storage() {
    let controls = TransientAdvanceControls {
        maximum_step_s: 0.03,
        ..Default::default()
    };
    let mut expected = layered_startup();
    let mut workspace = expected.adaptive_workspace(1 << 20).unwrap();
    let mut polls = 0;
    let success = expected
        .advance(0.1, controls, &mut workspace, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(success.accepted_steps > 1);
    let expected_state = snapshot(&expected);
    let mut actual = layered_startup();
    let original = snapshot(&actual);
    let mut calls = 0;
    assert!(matches!(
        actual.advance(0.1, controls, &mut workspace, || {
            calls += 1;
            calls == polls
        }),
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::Cancelled
        ))
    ));
    assert_eq!(calls, polls);
    assert_eq!(snapshot(&actual), original);
    let retried = actual
        .advance(0.1, controls, &mut workspace, || false)
        .unwrap();
    assert_eq!(snapshot(&actual), expected_state);
    assert_eq!(retried.root_evaluations, success.root_evaluations);
    // The next exact target is based on the already accepted floating-point time.
    let target = actual.time_s() + 0.2;
    let next = actual
        .advance(0.2, controls, &mut workspace, || false)
        .unwrap();
    assert_eq!(actual.time_s().to_bits(), target.to_bits());
    assert_eq!(next.interval_s[1].to_bits(), target.to_bits());
    assert_eq!(next.endpoint.interval_s[1].to_bits(), target.to_bits());
}

#[test]
fn adaptive_equilibrium_folds_a_one_ulp_final_remainder_before_solving() {
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[50.0], 1e-6, 0.0), bulk(&[50.0], 1e-6, 0.0)],
        1e-4,
        membrane(0.0, 0.0, 1e-3, 1),
        None,
        None,
        1 << 20,
    )
    .unwrap();
    let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
    let controls = TransientAdvanceControls {
        maximum_step_s: 0.01,
        ..Default::default()
    };
    let report = pair
        .advance(0.1, controls, &mut workspace, || false)
        .unwrap();
    assert_eq!(pair.time_s().to_bits(), 0.1_f64.to_bits());
    assert_eq!(report.endpoint.interval_s[1].to_bits(), 0.1_f64.to_bits());
    assert_eq!(report.accepted_steps, 10);
    assert_eq!(report.rejected_steps, 0);
    let endpoint_ulp = 0.1_f64.next_up() - 0.1;
    assert!(report.maximum_accepted_step_s <= controls.maximum_step_s + endpoint_ulp);
    assert_eq!(report.solvent_transfer_m3, 0.0);
    close(pair.total_amount(0), 1e-4, 1e-13, 0.0);
}

#[test]
fn adaptive_error_accounts_for_unequal_representable_half_steps() {
    let mut pair = diffusing_reservoirs();
    // The autonomous transport law is unchanged by a clock offset. At this
    // epoch one ULP is 0.25 s, so a 0.75 s interval splits into 0.5 + 0.25 s.
    pair.time_s = 2.0_f64.powi(50);
    let start = pair.time_s();
    let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
    let controls = TransientAdvanceControls {
        maximum_step_s: 0.75,
        relative_tolerance: 0.1,
        ..Default::default()
    };
    let report = pair
        .advance(0.75, controls, &mut workspace, || false)
        .unwrap();
    assert_eq!(report.accepted_steps, 1);
    assert_eq!(report.endpoint.interval_s, [start + 0.5, start + 0.75]);
    // Exact backward-Euler solutions of the original two-reservoir equations.
    // The initially empty right reservoir sets the strictest scaled error.
    let coarse = 50.0 * (1.0 - 1.0 / (1.0 + 0.2 * 0.75));
    let fine = 50.0 * (1.0 - 1.0 / ((1.0 + 0.2 * 0.5) * (1.0 + 0.2 * 0.25)));
    let expected = 1.25 * (fine - coarse) * 1e-6
        / (controls.absolute_amount_tolerance_mol + controls.relative_tolerance * fine * 1e-6);
    close(report.maximum_scaled_error, expected, 1e-12, 1e-14);
    close(
        pair.bulks()[1].solution().solute_moles()[0] / 1e-6,
        fine,
        1e-13,
        1e-13,
    );
}

#[test]
fn adaptive_budget_and_temporal_refusals_are_atomic() {
    let mut pair = diffusing_reservoirs();
    let initial = snapshot(&pair);
    let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
    assert!(matches!(
        pair.adaptive_workspace(workspace.workspace_bytes() - 1),
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::Workspace { .. }
        ))
    ));
    assert!(matches!(
        pair.advance(
            1.0,
            TransientAdvanceControls {
                maximum_step_s: 0.001,
                maximum_steps: 2,
                ..Default::default()
            },
            &mut workspace,
            || false
        ),
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::WorkBudget
        ))
    ));
    assert_eq!(snapshot(&pair), initial);
    assert!(matches!(
        pair.advance(
            10.0,
            TransientAdvanceControls {
                maximum_step_s: 10.0,
                minimum_step_s: 10.0,
                ..Default::default()
            },
            &mut workspace,
            || false
        ),
        Err(TransientOsmosisError::TemporalAccuracy { .. })
    ));
    assert_eq!(snapshot(&pair), initial);
    assert!(matches!(
        pair.advance(
            10.0,
            TransientAdvanceControls {
                maximum_step_s: 10.0,
                maximum_rejections: 0,
                ..Default::default()
            },
            &mut workspace,
            || false
        ),
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::WorkBudget
        ))
    ));
    assert_eq!(snapshot(&pair), initial);
    let mut different_shape = layered_startup();
    let different_initial = snapshot(&different_shape);
    assert!(
        different_shape
            .advance(
                1.0,
                TransientAdvanceControls::default(),
                &mut workspace,
                || false
            )
            .is_err()
    );
    assert_eq!(snapshot(&different_shape), different_initial);
}

#[test]
fn explicit_geometry_layout_and_workspace_admission() {
    assert!(TransientFilm::new(0.0, 2, vec![1e-9], vec![1.0; 2]).is_err());
    assert!(TransientFilm::new(1e-4, 0, vec![1e-9], vec![]).is_err());
    assert!(TransientFilm::new(1e-4, 2, vec![0.0], vec![1.0; 2]).is_err());
    assert!(TransientFilm::new(1e-4, 2, vec![1e-9], vec![1.0]).is_err());
    let pair = retaining_pair(4, 1e-6, 101_325.0);
    assert_eq!(pair.film_cell_count(0), Some(4));
    assert_eq!(pair.film_cell_count(2), None);
    assert!(pair.film_solute_moles(1, 4).is_none());
    close(pair.total_solvent_volume_m3(), 2.02e-6, 1e-15, 0.0);
    let mut total = [0.0];
    pair.total_solute_moles_into(&mut total).unwrap();
    close(total[0], 151.5e-6, 1e-15, 0.0);
    assert!(pair.total_solute_moles_into(&mut []).is_err());
    let refused = TransientOsmoticPair::new(
        [bulk(&[50.0], 1e-6, 0.0), bulk(&[100.0], 1e-6, 0.0)],
        1e-4,
        membrane(2e-11, 1.0, 0.0, 1),
        Some(film(4, 50.0)),
        Some(film(4, 100.0)),
        pair.workspace_bytes() - 1,
    );
    assert!(matches!(
        refused,
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::Workspace { .. }
        ))
    ));
}

#[test]
fn layered_cells_retain_local_storage_and_conserve_positive_inventories() {
    for (widths, diffusivities, profile) in [
        (vec![0.0], vec![1e-9], vec![1.0]),
        (vec![f64::NAN], vec![1e-9], vec![1.0]),
        (vec![1e-4], vec![0.0], vec![1.0]),
        (vec![1e-4, 2e-4], vec![1e-9], vec![1.0, 2.0]),
    ] {
        assert!(TransientFilm::from_cells(widths, 1, diffusivities, profile).is_err());
    }
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[0.0], 1e-6, 0.0), bulk(&[50.0], 1e-6, 0.0)],
        1e-4,
        membrane(0.0, 0.0, 2e-6, 1),
        Some(
            TransientFilm::from_cells(
                vec![1e-5, 3e-5, 6e-5],
                1,
                vec![1e-9, 1e-13, 2e-9],
                vec![100.0, 20.0, 0.0],
            )
            .unwrap(),
        ),
        Some(
            TransientFilm::from_cells(vec![2e-5, 8e-5], 1, vec![4e-10, 1e-9], vec![0.0, 50.0])
                .unwrap(),
        ),
        1 << 20,
    )
    .unwrap();
    close(pair.film_cell_volume_m3(0, 1).unwrap(), 3e-9, 1e-15, 0.0);
    close(pair.film_cell_volume_m3(1, 1).unwrap(), 8e-9, 1e-15, 0.0);
    assert_eq!(pair.film_cell_volume_m3(0, 3), None);
    assert_eq!(pair.film_cell_volume_m3(2, 0), None);
    let mut amount = [0.0];
    pair.total_solute_moles_into(&mut amount).unwrap();
    let original = 50e-6 + 1e-4 * (1e-5 * 100.0 + 3e-5 * 20.0 + 8e-5 * 50.0);
    close(amount[0], original, 1e-15, 0.0);
    close(pair.total_solvent_volume_m3(), 2.02e-6, 1e-15, 0.0);
    let before = snapshot(&pair);
    let mut polls = 0;
    assert!(matches!(
        pair.step(100.0, TransientStepControls::default(), || {
            polls += 1;
            polls == 10
        }),
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::Cancelled
        ))
    ));
    assert_eq!(snapshot(&pair), before);
    pair.step(100.0, TransientStepControls::default(), || false)
        .unwrap();
    pair.total_solute_moles_into(&mut amount).unwrap();
    close(amount[0], original, 1e-12, 1e-19);
    for side in 0..2 {
        for cell in 0..pair.film_cell_count(side).unwrap() {
            let c = pair.film_solute_moles(side, cell).unwrap()[0]
                / pair.film_cell_volume_m3(side, cell).unwrap();
            assert!((0.0..=100.0).contains(&c));
        }
    }
}

#[test]
fn heterogeneous_steady_flux_matches_series_resistance_at_both_membrane_faces() {
    let widths = [vec![0.5e-4, 1.5e-4], vec![0.2e-4, 1.8e-4]];
    let diffusivities = [vec![1e-9, 1e-11], vec![2e-10, 4e-9]];
    let resistance = [
        0.5e-4 / 1e-9 + 1.5e-4 / 1e-11,
        0.2e-4 / 2e-10 + 1.8e-4 / 4e-9,
    ];
    let permeability = 2e-6;
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[100.0], 1.0, 0.0), bulk(&[0.0], 1.0, 0.0)],
        1e-4,
        membrane(0.0, 0.0, permeability, 1),
        Some(
            TransientFilm::from_cells(widths[0].clone(), 1, diffusivities[0].clone(), vec![0.0; 2])
                .unwrap(),
        ),
        Some(
            TransientFilm::from_cells(widths[1].clone(), 1, diffusivities[1].clone(), vec![0.0; 2])
                .unwrap(),
        ),
        1 << 20,
    )
    .unwrap();
    let mut report = pair
        .step(1e4, TransientStepControls::default(), || false)
        .unwrap();
    for _ in 0..19 {
        report = pair
            .step(1e4, TransientStepControls::default(), || false)
            .unwrap();
    }
    let cl = pair.bulks()[0].solution().solute_moles()[0];
    let cr = pair.bulks()[1].solution().solute_moles()[0];
    let flux = (cl - cr) / (resistance[0] + 1.0 / permeability + resistance[1]);
    close(
        report.membrane_left_concentration_mol_m3[0],
        cl - flux * resistance[0],
        2e-7,
        1e-8,
    );
    close(
        report.membrane_right_concentration_mol_m3[0],
        cr + flux * resistance[1],
        2e-7,
        1e-8,
    );
    close(
        permeability
            * (report.membrane_left_concentration_mol_m3[0]
                - report.membrane_right_concentration_mol_m3[0]),
        flux,
        2e-7,
        1e-14,
    );
    assert_eq!(report.solvent_transfer_m3, 0.0);
    // The original 100 mol remains split across BOTH bulks and all four cells.
    let mut total = [0.0];
    pair.total_solute_moles_into(&mut total).unwrap();
    close(total[0], 100.0, 5e-14, 1e-13);
}

#[test]
fn layered_films_preserve_exponential_polarization_in_both_flow_directions() {
    let widths = [1e-5, 3e-5, 6e-5];
    let diffusivities = [2e-9, 5e-10, 1e-9];
    let resistance: f64 = widths.iter().zip(diffusivities).map(|(dx, d)| dx / d).sum();
    for pressure_right in [101_325.0, 501_325.0] {
        let mut pair = TransientOsmoticPair::new(
            [
                bulk(&[50.0], 1.0, 101_325.0),
                bulk(&[100.0], 1.0, pressure_right),
            ],
            1e-4,
            membrane(2e-11, 1.0, 0.0, 1),
            Some(
                TransientFilm::from_cells(
                    widths.to_vec(),
                    1,
                    diffusivities.to_vec(),
                    vec![50.0; 3],
                )
                .unwrap(),
            ),
            Some(
                TransientFilm::from_cells(
                    widths.into_iter().rev().collect(),
                    1,
                    diffusivities.into_iter().rev().collect(),
                    vec![100.0; 3],
                )
                .unwrap(),
            ),
            1 << 20,
        )
        .unwrap();
        let mut report = pair
            .step(100.0, TransientStepControls::default(), || false)
            .unwrap();
        for _ in 0..19 {
            report = pair
                .step(100.0, TransientStepControls::default(), || false)
                .unwrap();
        }
        let concentration = |side: usize| {
            let b = pair.bulks()[side].solution();
            b.solute_moles()[0] / b.solvent_volume_m3()
        };
        let peclet = report.volume_flux_m_s * resistance;
        close(
            report.membrane_left_concentration_mol_m3[0],
            concentration(0) * peclet.exp(),
            2e-6,
            1e-8,
        );
        close(
            report.membrane_right_concentration_mol_m3[0],
            concentration(1) * (-peclet).exp(),
            2e-6,
            1e-8,
        );
        assert_eq!(report.volume_flux_m_s < 0.0, pressure_right > 200_000.0);
    }
}

#[test]
fn refining_nonuniform_diffusion_layers_resolves_transient_uptake() {
    let uptake = |cells_per_layer: usize| {
        let mut widths = vec![1e-4 / cells_per_layer as f64; cells_per_layer];
        widths.extend(vec![2e-4 / cells_per_layer as f64; cells_per_layer]);
        let mut diffusivities = vec![1e-9; cells_per_layer];
        diffusivities.extend(vec![1e-10; cells_per_layer]);
        let mut pair = TransientOsmoticPair::new(
            [bulk(&[100.0], 1e-3, 0.0), bulk(&[0.0], 1e-3, 0.0)],
            1e-4,
            membrane(0.0, 1.0, 0.0, 1),
            Some(
                TransientFilm::from_cells(widths, 1, diffusivities, vec![0.0; 2 * cells_per_layer])
                    .unwrap(),
            ),
            None,
            1 << 20,
        )
        .unwrap();
        for _ in 0..250 {
            pair.step(0.2, TransientStepControls::default(), || false)
                .unwrap();
        }
        (0..2 * cells_per_layer)
            .map(|cell| pair.film_solute_moles(0, cell).unwrap()[0])
            .sum::<f64>()
    };
    let reference = uptake(32);
    let coarse = (uptake(1) - reference).abs();
    let fine = (uptake(2) - reference).abs();
    assert!(
        fine < 0.7 * coarse,
        "layered transient refinement errors: {coarse:e}, {fine:e}"
    );
}

#[test]
fn pure_water_counts_film_storage_and_exact_pressure_work() {
    let empty_film = || TransientFilm::new(1e-3, 3, vec![], vec![]).unwrap();
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[], 1e-3, 200_000.0), bulk(&[], 1e-3, 100_000.0)],
        1e-3,
        membrane(1e-10, 0.0, 0.0, 0),
        Some(empty_film()),
        Some(empty_film()),
        1 << 20,
    )
    .unwrap();
    let initial_water = pair.total_solvent_volume_m3();
    close(initial_water, 0.002002, 1e-15, 0.0);
    let report = pair
        .step(10.0, TransientStepControls::default(), || false)
        .unwrap();
    close(report.solvent_transfer_m3, 1e-7, 1e-12, 1e-20);
    close(report.pressure_work_into_solution_j, 0.01, 1e-12, 1e-16);
    close(pair.total_solvent_volume_m3(), initial_water, 1e-15, 0.0);
    close(
        pair.film_cell_volume_m3(0, 0).unwrap(),
        1e-6 / 3.0,
        1e-15,
        0.0,
    );
    assert_eq!(report.species_count, 0);
}

#[test]
fn no_film_limit_has_exact_two_compartment_backward_euler_diffusion() {
    let mut pair = TransientOsmoticPair::new(
        [
            bulk(&[0.0], 1e-6, 101_325.0),
            bulk(&[100.0], 1e-6, 101_325.0),
        ],
        1e-4,
        membrane(0.0, 0.0, 2e-6, 1),
        None,
        None,
        1 << 20,
    )
    .unwrap();
    let report = pair
        .step(1000.0, TransientStepControls::default(), || false)
        .unwrap();
    let expected_left = 0.5e-4 * (1.0 - 1.0 / 1.4);
    close(
        pair.bulks()[0].solution().solute_moles()[0],
        expected_left,
        1e-14,
        1e-20,
    );
    assert_eq!(report.solvent_transfer_m3, 0.0);
    assert_eq!(pair.film_cell_count(0), Some(0));
    close(report.solute_balance_residual_mol[0], 0.0, 0.0, 1e-19);
}

#[test]
fn retained_film_inventory_diffuses_into_finite_bulk_without_leaking_membrane() {
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[0.0], 1e-6, 0.0), bulk(&[0.0], 1e-6, 0.0)],
        1e-4,
        membrane(0.0, 1.0, 0.0, 1),
        Some(film(8, 100.0)),
        Some(film(8, 0.0)),
        1 << 20,
    )
    .unwrap();
    let mut before = [0.0];
    pair.total_solute_moles_into(&mut before).unwrap();
    close(before[0], 1e-6, 1e-15, 0.0);
    let report = pair
        .step(1000.0, TransientStepControls::default(), || false)
        .unwrap();
    let mut after = [0.0];
    pair.total_solute_moles_into(&mut after).unwrap();
    close(after[0], before[0], 1e-11, 1e-18);
    assert!(pair.bulks()[0].solution().solute_moles()[0] > 0.9e-6);
    assert_eq!(pair.bulks()[1].solution().solute_moles()[0], 0.0);
    for cell in 0..8 {
        assert_eq!(pair.film_solute_moles(1, cell).unwrap()[0], 0.0);
    }
    assert_eq!(report.solvent_transfer_m3, 0.0);
}

#[test]
fn permeating_species_with_distinct_diffusivities_conserve_all_bulk_and_film_moles() {
    let law = MembraneTransport::new(
        HydraulicPermeance::new(2e-11).unwrap(),
        vec![
            (
                ReflectionCoefficient::new(0.7).unwrap(),
                SolutePermeance::new(2e-7).unwrap(),
            ),
            (
                ReflectionCoefficient::new(0.1).unwrap(),
                SolutePermeance::new(4e-7).unwrap(),
            ),
        ],
    )
    .unwrap();
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[80.0, 0.0], 1e-6, 0.0), bulk(&[0.0, 30.0], 1e-6, 0.0)],
        1e-4,
        law,
        Some(TransientFilm::new(1e-4, 4, vec![1e-9, 5e-10], [80.0, 0.0].repeat(4)).unwrap()),
        Some(TransientFilm::new(2e-4, 3, vec![1e-9, 5e-10], [0.0, 30.0].repeat(3)).unwrap()),
        1 << 20,
    )
    .unwrap();
    let mut before = [0.0; 2];
    pair.total_solute_moles_into(&mut before).unwrap();
    let water = pair.total_solvent_volume_m3();
    let report = pair
        .step(100.0, TransientStepControls::default(), || false)
        .unwrap();
    let mut after = [0.0; 2];
    pair.total_solute_moles_into(&mut after).unwrap();
    for i in 0..2 {
        close(before[i], after[i], 1e-12, 1e-19);
    }
    assert!(pair.bulks()[0].solution().solute_moles()[1] > 0.0);
    assert!(pair.bulks()[1].solution().solute_moles()[0] > 0.0);
    assert_eq!(report.species_count, 2);
    close(pair.total_solvent_volume_m3(), water, 1e-14, 0.0);
    for side in 0..2 {
        for cell in 0..pair.film_cell_count(side).unwrap() {
            assert!(
                pair.film_solute_moles(side, cell)
                    .unwrap()
                    .iter()
                    .all(|n| n.is_finite() && *n >= 0.0)
            );
        }
    }
}

#[test]
fn nonideal_owner_is_refused_by_neutral_ideal_film_transport() {
    use fs_material::osmosis::{AqueousSolution, nonideal::VirialSolutionModel};
    let model = std::sync::Arc::new(
        VirialSolutionModel::new(298.15, 18e-6, 1.0, vec![200.0], vec![0.002]).unwrap(),
    );
    let left = OsmoticCompartment::new(
        AqueousSolution::nonideal(model, 1e-6, vec![50e-6]).unwrap(),
        OsmoticPressureVolume::new(1e-6, 0.0, 0.0).unwrap(),
        1e-9,
    )
    .unwrap();
    assert!(
        TransientOsmoticPair::new(
            [left, bulk(&[100.0], 1e-6, 0.0)],
            1e-4,
            membrane(2e-11, 1.0, 0.0, 1),
            Some(film(4, 50.0)),
            Some(film(4, 100.0)),
            1 << 20,
        )
        .is_err()
    );
}

#[test]
fn transient_profile_relaxes_to_original_steady_film_law_in_both_directions() {
    for pressure_right in [101_325.0, 501_325.0] {
        // Large but finite bulk volumes separate film relaxation from changing
        // bulk composition. Their water and moles still enter every balance.
        let mut pair = retaining_pair(8, 1.0, pressure_right);
        let law = PolarizedMembraneTransport::new(
            membrane(2e-11, 1.0, 0.0, 1),
            vec![Some(
                FilmPermeance::from_diffusion_layer(1e-9, 1e-4).unwrap(),
            )],
            vec![Some(
                FilmPermeance::from_diffusion_layer(1e-9, 1e-4).unwrap(),
            )],
        )
        .unwrap();
        let mut total_before = [0.0];
        pair.total_solute_moles_into(&mut total_before).unwrap();
        let mut report = pair
            .step(100.0, TransientStepControls::default(), || false)
            .unwrap();
        for _ in 0..9 {
            report = pair
                .step(100.0, TransientStepControls::default(), || false)
                .unwrap();
        }
        let bulks = pair.bulks();
        let cl = bulks[0].solution().solute_moles()[0] / bulks[0].solution().solvent_volume_m3();
        let cr = bulks[1].solution().solute_moles()[0] / bulks[1].solution().solvent_volume_m3();
        let (fl, fr) = law
            .coefficients(0, report.volume_flux_m_s)
            .unwrap()
            .membrane_face_concentrations(cl, cr)
            .unwrap();
        close(report.membrane_left_concentration_mol_m3[0], fl, 2e-6, 1e-8);
        close(
            report.membrane_right_concentration_mol_m3[0],
            fr,
            2e-6,
            1e-8,
        );
        let mut total_after = [0.0];
        pair.total_solute_moles_into(&mut total_after).unwrap();
        close(total_before[0], total_after[0], 2e-13, 1e-15);
        if pressure_right > 200_000.0 {
            assert!(report.volume_flux_m_s < 0.0);
            assert!(report.pressure_work_into_solution_j > 0.0);
            assert!(fl < cl && fr > cr);
        } else {
            assert!(report.volume_flux_m_s > 0.0);
            assert!(fl > cl && fr < cr);
        }
    }
}

fn uptake(cells: usize, step: f64, duration: f64) -> f64 {
    let mut pair = retaining_pair(cells, 1e-6, 101_325.0);
    let controls = TransientStepControls {
        absolute_water_tolerance_m3: 1e-22,
        relative_root_tolerance: 1e-11,
        ..TransientStepControls::default()
    };
    for _ in 0..(duration / step).round() as usize {
        pair.step(step, controls, || false).unwrap();
    }
    pair.bulks()[1].solution().solvent_volume_m3() - 1e-6
}

#[test]
fn time_and_space_refinement_resolve_film_startup() {
    let reference_time = uptake(8, 0.03125, 4.0);
    let coarse_time = uptake(8, 0.5, 4.0);
    let fine_time = uptake(8, 0.25, 4.0);
    assert!((fine_time - reference_time).abs() < 0.7 * (coarse_time - reference_time).abs());
    let reference_space = uptake(32, 0.03125, 4.0);
    let coarse_space = uptake(2, 0.03125, 4.0);
    let fine_space = uptake(4, 0.03125, 4.0);
    assert!((fine_space - reference_space).abs() < 0.7 * (coarse_space - reference_space).abs());
}

#[test]
fn steep_water_law_accepts_a_resolved_transfer_bracket_without_hiding_its_residual() {
    // Deliberately extreme numerical stiffness: the source law remains intact,
    // while the root's transfer-location tolerance can outresolve its residual.
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[50.0], 1e-6, 0.0), bulk(&[100.0], 1e-6, 0.0)],
        1e-4,
        membrane(1e6, 1.0, 0.0, 1),
        Some(film(8, 50.0)),
        Some(film(8, 100.0)),
        1 << 20,
    )
    .unwrap();
    let report = pair
        .step(
            10.0,
            TransientStepControls {
                absolute_water_tolerance_m3: 1e-24,
                relative_root_tolerance: 1e-13,
                ..TransientStepControls::default()
            },
            || false,
        )
        .unwrap();
    assert!(report.solvent_transfer_m3 > 0.0);
    assert!(report.water_root_residual_m3.is_finite());
    assert!(report.water_root_residual_m3.abs() > 1e-24 + 1e-13 * report.solvent_transfer_m3);
    close(report.solute_balance_residual_mol[0], 0.0, 0.0, 1e-18);
}

#[test]
fn cancellation_and_root_exhaustion_preserve_all_finite_inventories_and_clock() {
    let mut pair = retaining_pair(8, 1e-6, 101_325.0);
    let before = snapshot(&pair);
    let mut polls = 0;
    let result = pair.step(10.0, TransientStepControls::default(), || {
        polls += 1;
        polls == 10
    });
    assert!(matches!(
        result,
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::Cancelled
        ))
    ));
    assert_eq!(snapshot(&pair), before);
    let result = pair.step(
        10.0,
        TransientStepControls {
            maximum_root_evaluations: 2,
            ..TransientStepControls::default()
        },
        || false,
    );
    assert!(matches!(
        result,
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::NonlinearConvergence
        ))
    ));
    assert_eq!(snapshot(&pair), before);
    pair.step(10.0, TransientStepControls::default(), || false)
        .unwrap();
    assert_eq!(pair.time_s(), 10.0);
}
