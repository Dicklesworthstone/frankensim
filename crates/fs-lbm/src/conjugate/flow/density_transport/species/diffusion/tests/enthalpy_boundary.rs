//! Original reservoir maps with physical partial-enthalpy receipts and audits.
use super::*;
use std::cell::Cell;

fn reservoir(cell: usize, face: usize, fractions: [f64; 2]) -> SpeciesDiffusionBoundary<2> {
    SpeciesDiffusionBoundary {
        cell,
        face,
        reservoir_mass_fractions: fractions,
        film_resistance_m2_s_kg: 0.0,
    }
}

#[test]
fn reservoir_enthalpy_matches_analytic_film_exchange_and_original_species_maps() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[2.0]);
    let species = [[1.6, 0.4]];
    let h = [900.0];
    let mut boundary = reservoir(0, 1, [0.25, 0.75]);
    boundary.film_resistance_m2_s_kg = 0.75;
    let boundaries = [boundary];
    let gate = CancelGate::new_clock_free();
    let mass = state
        .diffuse_species_with_boundaries(
            &domain,
            &species,
            &[0.4],
            &boundaries,
            0.6,
            controls(7),
            &gate,
        )
        .unwrap();
    let calls = Cell::new(0);
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries(
            &domain,
            &species,
            &h,
            &[0.4],
            &boundaries,
            0.6,
            |_, _| panic!("a single cell has no interior face"),
            |index| {
                assert_eq!(index, 0);
                calls.set(calls.get() + 1);
                Ok([1000.0, 500.0])
            },
            controls(7),
            1e-9,
            &gate,
        )
        .unwrap();
    assert_eq!(calls.get(), 1);
    assert_eq!(got.step.species, mass.step);
    assert_eq!(got.boundary, mass.boundary);
    // rho=2, half-cell resistance=.5/(rho*D), G=A/(R_half+R_film).
    // The fixed reservoir solves dY/dt=(G/M)*(Y_res-Y) exactly.
    let conductance = 1.0 / (0.5 / (2.0 * 0.4) + 0.75);
    let relaxation = -(-conductance * 0.6 / 2.0_f64).exp_m1();
    let dm_a = 2.0 * (0.25 - 0.8) * relaxation;
    let heat = 500.0 * dm_a;
    near(got.step.species.species_mass_kg[0][0], 1.6 + dm_a, 2e-14);
    near(
        got.step.specific_enthalpy_j_kg[0],
        900.0 + heat / 2.0,
        2e-12,
    );
    near(got.thermal_boundary.enthalpy_into_domain_j, heat, 2e-12);
    near(got.step.report.stored_enthalpy_change_j, heat, 2e-12);
    near(got.step.report.transported_enthalpy_change_j, heat, 2e-12);
    assert_eq!(got.thermal_boundary.positive_enthalpy_into_domain_j, 0.0);
    near(
        got.thermal_boundary.negative_enthalpy_into_domain_j,
        -heat,
        2e-12,
    );
    assert_eq!(state.mass_kg, [2.0]);
    assert_eq!(species, [[1.6, 0.4]]);
    assert_eq!(h, [900.0]);
}

#[test]
fn reservoir_and_interior_enthalpy_share_one_original_mass_candidate() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0, 2.0]);
    let species = [[0.8, 0.2], [0.2, 1.8]];
    let h = [768.0, -128.0];
    let d = [0.3, 0.1];
    let boundaries = [reservoir(0, 0, [0.25, 0.75]), reservoir(1, 1, [0.75, 0.25])];
    let gate = CancelGate::new_clock_free();
    let interior_calls = Cell::new(0);
    let boundary_calls = Cell::new(0);
    let mass = state
        .diffuse_species_with_boundaries(
            &domain,
            &species,
            &d,
            &boundaries,
            1.0,
            controls(4),
            &gate,
        )
        .unwrap();
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries(
            &domain,
            &species,
            &h,
            &d,
            &boundaries,
            1.0,
            |cells, axis| {
                assert_eq!(cells, [0, 1]);
                assert_eq!(axis, 0);
                interior_calls.set(interior_calls.get() + 1);
                Ok([1024.0, -256.0])
            },
            |_| {
                boundary_calls.set(boundary_calls.get() + 1);
                Ok([1024.0, -256.0])
            },
            controls(4),
            1e-9,
            &gate,
        )
        .unwrap();
    assert_eq!(interior_calls.get(), 1);
    assert_eq!(boundary_calls.get(), 2);
    assert_eq!(got.step.species, mass.step);
    assert_eq!(got.boundary, mass.boundary);
    for c in 0..2 {
        let inventory = got.step.species.species_mass_kg[c];
        let expected = (1024.0 * inventory[0] - 256.0 * inventory[1]) / state.mass_kg[c];
        near(got.step.specific_enthalpy_j_kg[c], expected, 2e-10);
    }
    for (thermal, mass) in got.thermal_boundary.faces.iter().zip(&got.boundary.faces) {
        near(
            thermal.enthalpy_into_cell_j,
            1280.0 * mass.mass_into_cell_kg[0],
            2e-10,
        );
    }
    // No declared reservoirs reduces exactly to the existing closed entry,
    // including its independent local/global thermal and species reports.
    let closed = state
        .diffuse_species_with_enthalpy(
            &domain,
            &species,
            &h,
            &d,
            1.0,
            |_, _| Ok([1024.0, -256.0]),
            controls(4),
            1e-9,
            &gate,
        )
        .unwrap();
    let no_boundaries = state
        .diffuse_species_with_enthalpy_and_boundaries(
            &domain,
            &species,
            &h,
            &d,
            &[],
            1.0,
            |_, _| Ok([1024.0, -256.0]),
            |_| panic!("an absent reservoir must not query its source"),
            controls(4),
            1e-9,
            &gate,
        )
        .unwrap();
    assert_eq!(no_boundaries.step, closed);
    assert!(no_boundaries.thermal_boundary.faces.is_empty());
    assert_eq!(no_boundaries.thermal_boundary.enthalpy_into_domain_j, 0.0);
}

#[test]
fn reservoir_heat_reversals_keep_gross_signs_and_common_datum_invariance() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let boundaries = [
        reservoir(0, 0, [0.0, 1.0]),
        reservoir(0, 1, [0.5, 0.5]),
        reservoir(0, 2, [1.0, 0.0]),
    ];
    let gate = CancelGate::new_clock_free();
    let run = |datum: f64| {
        state
            .diffuse_species_with_enthalpy_and_boundaries(
                &domain,
                &[[0.5, 0.5]],
                &[datum + 4.0],
                &[1e30],
                &boundaries,
                1.0,
                |_, _| unreachable!(),
                |_| Ok([datum + 8.0, datum]),
                controls(1),
                0.0,
                &gate,
            )
            .unwrap()
    };
    let got = run(0.0);
    // Forward reservoirs give Y_A=0,.5,1; reverse gives 1,.5,0.
    // The middle face therefore carries +4 J and -4 J despite zero net heat.
    assert_eq!(got.step.species.species_mass_kg, [[0.0, 1.0]]);
    assert_eq!(got.step.specific_enthalpy_j_kg, [0.0]);
    let middle = &got.thermal_boundary.faces[1];
    assert_eq!(middle.enthalpy_into_cell_j, 0.0);
    assert_eq!(middle.positive_enthalpy_into_cell_j, 4.0);
    assert_eq!(middle.negative_enthalpy_into_cell_j, 4.0);
    assert_eq!(got.thermal_boundary.enthalpy_into_domain_j, -4.0);
    assert_eq!(got.thermal_boundary.positive_enthalpy_into_domain_j, 8.0);
    assert_eq!(got.thermal_boundary.negative_enthalpy_into_domain_j, 12.0);
    assert_eq!(got.step.report.absolute_enthalpy_defect_j, 0.0);
    let datum = 2.0_f64.powi(50);
    let shifted = run(datum);
    assert_eq!(shifted.step.specific_enthalpy_j_kg, [datum]);
    assert_eq!(shifted.thermal_boundary, got.thermal_boundary);
    assert_eq!(shifted.step.report, got.step.report);
    assert_eq!(shifted.boundary, got.boundary);
}

#[test]
fn reservoir_formation_cancellation_retains_external_and_species_product_tails() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let gate = CancelGate::new_clock_free();
    let formation = 2.0_f64.powi(54);
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries(
            &domain,
            &[[1.0, 0.0]],
            &[-formation],
            &[1e30],
            &[reservoir(0, 0, [0.5, 0.5])],
            1.0,
            |_, _| unreachable!(),
            |_| Ok([-formation, formation + 4.0]),
            controls(1),
            0.0,
            &gate,
        )
        .unwrap();
    // Actual heat is (F, +2). Both public energy scalars round to F, but the
    // exact represented endpoint and full external account agree at h=2.
    assert_eq!(got.step.specific_enthalpy_j_kg, [2.0]);
    assert_eq!(got.thermal_boundary.enthalpy_into_domain_j, formation);
    assert_eq!(got.step.diffusive_enthalpy_change_j, [formation]);
    assert_eq!(got.step.report.absolute_enthalpy_defect_j, 0.0);

    // Opposing ternary chemical products also retain the original difference
    // low part: (-.25)*(2F) + .25*(2F+4) is exactly one joule.
    let boundary = SpeciesDiffusionBoundary {
        cell: 0,
        face: 0,
        reservoir_mass_fractions: [0.5, 0.125, 0.375],
        film_resistance_m2_s_kg: 0.0,
    };
    let ternary = SpeciesDiffusionControls::<3> {
        substeps: 1,
        maximum_pair_updates: 10,
        maximum_storage_bytes: 1 << 20,
        maximum_mass_defect_kg: [0.0; 3],
        maximum_partition_defect_kg: 0.0,
    };
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries(
            &domain,
            &[[0.5, 0.375, 0.125]],
            &[0.5],
            &[1e30],
            &[boundary],
            1.0,
            |_, _| unreachable!(),
            |_| Ok([-formation, formation, formation + 4.0]),
            ternary,
            0.0,
            &gate,
        )
        .unwrap();
    assert_eq!(got.step.specific_enthalpy_j_kg, [1.5]);
    assert_eq!(got.thermal_boundary.enthalpy_into_domain_j, 1.0);
    assert_eq!(got.step.report.absolute_enthalpy_defect_j, 0.0);
}

#[test]
fn reservoir_enthalpy_charges_unrepresentable_storage_in_absolute_joules() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let gate = CancelGate::new_clock_free();
    let h = [2.0_f64.powi(54)];
    let run = |budget| {
        state.diffuse_species_with_enthalpy_and_boundaries(
            &domain,
            &[[0.5, 0.5]],
            &h,
            &[1e30],
            &[reservoir(0, 0, [1.0, 0.0])],
            1.0,
            |_, _| unreachable!(),
            |_| Ok([2.0, 0.0]),
            controls(1),
            budget,
            &gate,
        )
    };
    assert!(run(1.0).is_err());
    let got = run(2.0).unwrap();
    assert_eq!(got.step.specific_enthalpy_j_kg, h);
    assert_eq!(got.thermal_boundary.enthalpy_into_domain_j, 1.0);
    assert_eq!(got.step.report.stored_enthalpy_change_j, 0.0);
    // One joule local mismatch plus one joule domain storage/external mismatch.
    assert_eq!(got.step.report.absolute_enthalpy_defect_j, 2.0);
}

#[test]
fn reservoir_enthalpy_only_queries_active_faces_and_refuses_atomically() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let species = [[0.5, 0.5]; 2];
    let h = [1.0; 2];
    let boundaries = [reservoir(0, 0, [0.5, 0.5]), reservoir(1, 1, [0.5, 0.5])];
    let calls = Cell::new(0);
    let gate = CancelGate::new_clock_free();
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries(
            &domain,
            &species,
            &h,
            &[0.0, 0.2],
            &boundaries,
            1.0,
            |_, _| panic!("zero D blocks the interior face"),
            |index| {
                assert_eq!(index, 1);
                calls.set(calls.get() + 1);
                Ok([2.0, 0.0])
            },
            controls(2),
            0.0,
            &gate,
        )
        .unwrap();
    assert_eq!(calls.get(), 1);
    assert_eq!(got.step.species.species_mass_kg, species);
    assert_eq!(got.thermal_boundary.faces.len(), 2);
    assert_eq!(got.thermal_boundary.faces[0].enthalpy_into_cell_j, 0.0);
    for storage in [false, true] {
        calls.set(0);
        let mut limited = controls(1);
        if storage {
            limited.maximum_storage_bytes = 0;
        } else {
            limited.maximum_pair_updates = 0;
        }
        assert!(
            state
                .diffuse_species_with_enthalpy_and_boundaries(
                    &domain,
                    &species,
                    &h,
                    &[0.0, 0.2],
                    &boundaries,
                    1.0,
                    |_, _| unreachable!(),
                    |_| {
                        calls.set(calls.get() + 1);
                        Ok([2.0, 0.0])
                    },
                    limited,
                    0.0,
                    &gate,
                )
                .is_err()
        );
        assert_eq!(
            calls.get(),
            0,
            "resources must be admitted before original source calls"
        );
    }
    let error = refuse("original reservoir caloric source refused");
    assert_eq!(
        state
            .diffuse_species_with_enthalpy_and_boundaries(
                &domain,
                &species,
                &h,
                &[0.0, 0.2],
                &boundaries,
                1.0,
                |_, _| unreachable!(),
                |_| Err(error.clone()),
                controls(1),
                0.0,
                &gate,
            )
            .unwrap_err(),
        error
    );
    assert!(
        state
            .diffuse_species_with_enthalpy_and_boundaries(
                &domain,
                &species,
                &h,
                &[0.0, 0.2],
                &boundaries,
                1.0,
                |_, _| unreachable!(),
                |_| {
                    gate.request();
                    Ok([2.0, 0.0])
                },
                controls(1),
                0.0,
                &gate,
            )
            .is_err()
    );
    assert_eq!(state.mass_kg, [1.0; 2]);
    assert_eq!(species, [[0.5, 0.5]; 2]);
    assert_eq!(h, [1.0; 2]);
}
