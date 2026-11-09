//! G1/G3/G4 open species/enthalpy diffusion against independent physical limits.
//! Coefficients and partial enthalpies are authored numerical fixtures.
use super::*;

fn accuracy<const N: usize>(h: f64, energy: f64) -> BoundarySpeciesEnthalpyAccuracy<N> {
    BoundarySpeciesEnthalpyAccuracy {
        absolute_boundary_transfer_tolerance_kg: [1.0; N],
        absolute_specific_enthalpy_tolerance_j_kg: h,
        absolute_boundary_enthalpy_tolerance_j: energy,
    }
}

#[test]
fn g1_reservoir_enthalpy_matches_single_cell_film_relaxation_and_caloric_inventory() {
    let domain = VoxelDomain::graded([vec![0.7], vec![1.2], vec![2.1]]).unwrap();
    let mass = 4.0;
    let state = carrier(&domain, &[mass]);
    let before = state.clone();
    let initial = [[2.4, 1.2, 0.4]];
    let initial_h = [90_000.0];
    let partial = [100_000.0, 200_000.0, -300_000.0];
    let bath = [0.1, 0.3, 0.6];
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        for film in [0.0, 3.0, 1e6] {
            let face = 2 * axis;
            let boundary = SpeciesDiffusionBoundary {
                cell: 0,
                face,
                reservoir_mass_fractions: bath,
                film_resistance_m2_s_kg: film,
            };
            let density = mass / domain.volume(0);
            let g = domain.face_area(0, axis)
                / (0.5 * domain.widths(0)[axis] / (density * 0.06) + film);
            for duration in [1e-4, 0.7, 1e6] {
                let calls = std::cell::Cell::new(0);
                let result = state
                    .diffuse_species_with_enthalpy_and_boundaries(
                        &domain,
                        &initial,
                        &initial_h,
                        &[0.06],
                        &[boundary],
                        duration,
                        |_, _| panic!("one cell has no interior face"),
                        |index| {
                            assert_eq!(index, 0);
                            calls.set(calls.get() + 1);
                            Ok(partial)
                        },
                        common_controls(1),
                        1e-6,
                        &gate,
                    )
                    .unwrap();
                let relaxation = -(-g * duration / mass).exp_m1();
                let expected_h = 90_000.0 - 200_000.0 * relaxation;
                let expected_energy = mass * (expected_h - initial_h[0]);
                near(result.step.specific_enthalpy_j_kg[0], expected_h, 2e-9);
                near(
                    result.thermal_boundary.enthalpy_into_domain_j,
                    expected_energy,
                    8e-9,
                );
                near(
                    result.step.diffusive_enthalpy_change_j[0],
                    expected_energy,
                    8e-9,
                );
                let transfer = &result.thermal_boundary.faces[0];
                near(transfer.enthalpy_into_cell_j, expected_energy, 8e-9);
                near(transfer.positive_enthalpy_into_cell_j, 0.0, 1e-9);
                near(
                    transfer.negative_enthalpy_into_cell_j,
                    -expected_energy,
                    8e-9,
                );
                near(
                    result.step.report.stored_enthalpy_change_j,
                    expected_energy,
                    8e-9,
                );
                assert!(result.step.report.absolute_enthalpy_defect_j <= 1e-6);
                let caloric = result.step.species.species_mass_kg[0]
                    .iter()
                    .zip(partial)
                    .map(|(&m, h)| m * h)
                    .sum::<f64>()
                    / mass;
                near(result.step.specific_enthalpy_j_kg[0], caloric, 2e-9);
                assert_eq!(calls.get(), 1);
                assert_eq!(result.step.species.report.pair_updates, 2);
            }
        }
    }
    assert_eq!(state, before);
    assert_eq!(initial_h, [90_000.0]);
    assert_eq!(initial, [[2.4, 1.2, 0.4]]);
}

#[test]
fn g3_boundary_enthalpy_recovers_closed_owner_and_common_formation_datum() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0, 2.0]);
    let initial = [[0.75, 0.25], [0.5, 1.5]];
    let h = [265.0, 755.0];
    let gate = CancelGate::new_clock_free();
    let closed = state
        .diffuse_species_with_enthalpy(
            &domain,
            &initial,
            &h,
            &[0.2; 2],
            0.7,
            |_, _| Ok([20.0, 1000.0]),
            common_controls(3),
            1e-8,
            &gate,
        )
        .unwrap();
    let absent = state
        .diffuse_species_with_enthalpy_and_boundaries(
            &domain,
            &initial,
            &h,
            &[0.2; 2],
            &[],
            0.7,
            |_, _| Ok([20.0, 1000.0]),
            |_| panic!("no boundary source"),
            common_controls(3),
            1e-8,
            &gate,
        )
        .unwrap();
    assert_eq!(absent.step.species, closed.species);
    assert_eq!(
        absent.step.specific_enthalpy_j_kg,
        closed.specific_enthalpy_j_kg
    );
    assert_eq!(
        absent.step.diffusive_enthalpy_change_j,
        closed.diffusive_enthalpy_change_j
    );
    assert_eq!(absent.step.report, closed.report);
    assert_eq!(absent.thermal_boundary.enthalpy_into_domain_j, 0.0);

    let boundaries = [SpeciesDiffusionBoundary {
        cell: 0,
        face: 0,
        reservoir_mass_fractions: [0.2, 0.8],
        film_resistance_m2_s_kg: 0.5,
    }];
    let run = |datum: f64| {
        state
            .diffuse_species_with_enthalpy_and_boundaries(
                &domain,
                &initial,
                &h.map(|v| v + datum),
                &[0.2; 2],
                &boundaries,
                0.7,
                |_, _| Ok([20.0 + datum, 1000.0 + datum]),
                |_| Ok([20.0 + datum, 1000.0 + datum]),
                common_controls(5),
                0.01,
                &gate,
            )
            .unwrap()
    };
    let base = run(0.0);
    let datum = 2.0_f64.powi(40);
    let shifted = run(datum);
    assert_eq!(base.step.species, shifted.step.species);
    assert_eq!(
        base.step.diffusive_enthalpy_change_j,
        shifted.step.diffusive_enthalpy_change_j
    );
    assert_eq!(base.thermal_boundary.faces, shifted.thermal_boundary.faces);
    assert_eq!(
        base.thermal_boundary.enthalpy_into_domain_j,
        shifted.thermal_boundary.enthalpy_into_domain_j
    );
    for (a, b) in base
        .step
        .specific_enthalpy_j_kg
        .iter()
        .zip(&shifted.step.specific_enthalpy_j_kg)
    {
        near(b - datum, *a, 0.001);
    }
}

#[test]
fn g1_adaptive_opposing_reservoirs_resolve_energy_exchange_despite_zero_net_energy() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let initial = [[0.5, 0.5]];
    let boundaries = [0, 1].map(|face| SpeciesDiffusionBoundary {
        cell: 0,
        face,
        reservoir_mass_fractions: if face == 0 { [0.0, 1.0] } else { [1.0, 0.0] },
        film_resistance_m2_s_kg: 0.0,
    });
    let run = |energy| {
        state
            .diffuse_species_with_enthalpy_and_boundaries_adaptive(
                &domain,
                &initial,
                &[150_000.0],
                &[0.1],
                &boundaries,
                2.0,
                |_, _| panic!("one cell has no interior face"),
                |_| Ok([-200_000.0, 500_000.0]),
                adaptive_controls([1.0; 2]),
                accuracy(1e9, energy),
                1e-6,
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let loose = run(1e5);
    let tight = run(1e-3);
    assert!(tight.refinements > loose.refinements);
    assert!(tight.maximum_error_ratio <= 1.0);
    assert!(tight.estimated_boundary_enthalpy_error_j <= 1e-3);
    near(tight.step.step.specific_enthalpy_j_kg[0], 150_000.0, 0.005);
    near(
        tight.step.thermal_boundary.faces[0].enthalpy_into_cell_j,
        140_000.0,
        0.005,
    );
    near(
        tight.step.thermal_boundary.faces[1].enthalpy_into_cell_j,
        -140_000.0,
        0.005,
    );
    near(
        tight.step.thermal_boundary.enthalpy_into_domain_j,
        0.0,
        0.005,
    );
    near(
        tight.step.thermal_boundary.faces[0].negative_enthalpy_into_cell_j,
        0.0,
        1e-8,
    );
    near(
        tight.step.thermal_boundary.faces[1].positive_enthalpy_into_cell_j,
        0.0,
        1e-8,
    );
    assert_eq!(
        tight.total_pair_updates,
        4 * (2 * tight.step.step.species.report.substeps - 1),
    );
}

#[test]
fn g1_adaptive_caloric_tolerance_resolves_independent_neumann_modes() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 3]);
    let initial = [[0.8, 0.2], [0.2, 0.8], [0.4, 0.6]];
    let initial_h = initial.map(|m| 1000.0 + 8192.0 * m[0]);
    let result = state
        .diffuse_species_with_enthalpy_and_boundaries_adaptive(
            &domain,
            &initial,
            &initial_h,
            &[0.1; 3],
            &[],
            2.0,
            |_, _| Ok([9192.0, 1000.0]),
            |_| panic!("no reservoir"),
            adaptive_controls([1.0; 2]),
            accuracy(1e-4, 1.0),
            1e-7,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // Three-cell native Neumann eigenvalues are D and 3D. The coefficient
    // vectors below come from the independent mean/odd/even decomposition.
    let mean = 7.0 / 15.0;
    let odd = 0.2 * (-0.2_f64).exp();
    let even = (4.0 / 15.0) * (-0.6_f64).exp();
    let expected = [
        mean + odd + 0.5 * even,
        mean - even,
        mean - odd + 0.5 * even,
    ];
    for (cell, y) in expected.into_iter().enumerate() {
        near(
            result.step.step.specific_enthalpy_j_kg[cell],
            1000.0 + 8192.0 * y,
            8e-4,
        );
    }
    assert!(result.refinements > 1);
    assert!(result.estimated_specific_enthalpy_error_j_kg <= 1e-4);
    assert_eq!(result.step.thermal_boundary.enthalpy_into_domain_j, 0.0);
}

#[test]
fn g1_adaptive_caloric_error_survives_a_large_common_stored_datum() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 3]);
    let epsilon = 2.0_f64.powi(-20);
    let initial = [
        [0.5 + epsilon, 0.5 - epsilon],
        [0.5 - epsilon, 0.5 + epsilon],
        [0.5, 0.5],
    ];
    let datum = 2.0_f64.powi(80);
    let contrast = 2.0_f64.powi(40);
    let h = [datum + 0.5 * contrast; 3];
    // Partial-enthalpy contrasts remain representable. Transported changes
    // are smaller than the endpoint h ULP; the intentionally wide physical
    // storage-defect allowance must not relax the independent time tolerance.
    let run = |tolerance| {
        state
            .diffuse_species_with_enthalpy_and_boundaries_adaptive(
                &domain,
                &initial,
                &h,
                &[0.1; 3],
                &[],
                2.0,
                |_, _| Ok([datum + contrast, datum]),
                |_| panic!("no reservoir"),
                adaptive_controls([1.0; 2]),
                accuracy(tolerance, 1.0),
                1e9,
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let loose = run(1e9);
    let tight = run(1.0);
    assert_eq!(loose.step.step.specific_enthalpy_j_kg, h);
    assert_eq!(tight.step.step.specific_enthalpy_j_kg, h);
    assert_eq!(loose.refinements, 1);
    assert!(tight.refinements > loose.refinements);
    assert!(tight.estimated_specific_enthalpy_error_j_kg > 0.0);
    assert!(tight.estimated_specific_enthalpy_error_j_kg <= 1.0);
    assert!(tight.step.step.report.absolute_enthalpy_defect_j > 0.0);
}

#[test]
fn g0_g4_boundary_caloric_sources_barriers_budgets_and_cancellation_are_atomic() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let before = state.clone();
    let initial = [[0.2, 0.8], [0.8, 0.2]];
    let h = [804.0, 216.0];
    let boundaries = [(0, 0), (1, 1)].map(|(cell, face)| SpeciesDiffusionBoundary {
        cell,
        face,
        reservoir_mass_fractions: [0.5, 0.5],
        film_resistance_m2_s_kg: 0.0,
    });
    let gate = CancelGate::new_clock_free();
    let calls = std::cell::Cell::new(0);
    let result = state
        .diffuse_species_with_enthalpy_and_boundaries(
            &domain,
            &initial,
            &h,
            &[0.0, 0.1],
            &boundaries,
            1.0,
            |_, _| panic!("zero D blocks the interior"),
            |index| {
                assert_eq!(index, 1);
                calls.set(calls.get() + 1);
                Ok([20.0, 1000.0])
            },
            common_controls(1),
            1e-9,
            &gate,
        )
        .unwrap();
    assert_eq!(calls.get(), 1);
    assert_eq!(result.step.specific_enthalpy_j_kg[0], h[0]);
    assert_eq!(result.thermal_boundary.faces[0].enthalpy_into_cell_j, 0.0);
    assert!(result.thermal_boundary.faces[1].enthalpy_into_cell_j > 0.0);
    let broken = state.diffuse_species_with_enthalpy_and_boundaries(
        &domain,
        &initial,
        &h,
        &[0.0, 0.1],
        &boundaries,
        1.0,
        |_, _| unreachable!(),
        |_| Ok([0.0, f64::NAN]),
        common_controls(1),
        1e-9,
        &gate,
    );
    assert!(broken.is_err());
    for (pairs, bytes) in [(1, 1 << 22), (2, 1)] {
        let mut control = common_controls(1);
        control.maximum_pair_updates = pairs;
        control.maximum_storage_bytes = bytes;
        assert!(
            state
                .diffuse_species_with_enthalpy_and_boundaries(
                    &domain,
                    &initial,
                    &h,
                    &[0.0, 0.1],
                    &boundaries,
                    1.0,
                    |_, _| unreachable!(),
                    |_| Ok([20.0, 1000.0]),
                    control,
                    1e-9,
                    &gate,
                )
                .is_err()
        );
    }
    let mut control = adaptive_controls([1.0; 2]);
    control.solver.maximum_pair_updates = 5;
    assert!(
        state
            .diffuse_species_with_enthalpy_and_boundaries_adaptive(
                &domain,
                &initial,
                &h,
                &[0.0, 0.1],
                &boundaries,
                1.0,
                |_, _| unreachable!(),
                |_| Ok([20.0, 1000.0]),
                control,
                accuracy(1.0, 1.0),
                1e-9,
                &gate,
            )
            .is_err()
    );
    gate.request();
    assert_eq!(
        state
            .diffuse_species_with_enthalpy_and_boundaries(
                &domain,
                &initial,
                &h,
                &[0.0, 0.1],
                &boundaries,
                1.0,
                |_, _| unreachable!(),
                |_| Ok([20.0, 1000.0]),
                common_controls(1),
                1e-9,
                &gate,
            )
            .unwrap_err(),
        ChtError::Cancelled
    );
    assert_eq!(state, before);
    assert_eq!(initial, [[0.2, 0.8], [0.8, 0.2]]);
    assert_eq!(h, [804.0, 216.0]);
}
