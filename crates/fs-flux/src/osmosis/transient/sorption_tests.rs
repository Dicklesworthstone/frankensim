use super::*;

#[test]
fn linear_sorption_admission_and_r_one_preserve_the_original_path() {
    let definition = film(2, 50.0);
    for ratios in [
        vec![],
        vec![2.0],
        vec![0.5, 1.0],
        vec![1.0, f64::NAN],
        vec![f64::INFINITY, 1.0],
    ] {
        assert!(definition.clone().with_linear_sorption(ratios).is_err());
    }
    let prepare = |left| {
        TransientOsmoticPair::new(
            [
                bulk(&[50.0], 1e-6, 101_325.0),
                bulk(&[100.0], 1e-6, 101_325.0),
            ],
            1e-4,
            membrane(2e-11, 1.0, 0.0, 1),
            Some(left),
            Some(film(2, 100.0)),
            1 << 20,
        )
        .unwrap()
    };
    let mut original = prepare(definition.clone());
    let mut equivalent = prepare(
        definition
            .clone()
            .with_linear_sorption(vec![1.0; 2])
            .unwrap(),
    );
    assert_eq!(equivalent.workspace_bytes(), original.workspace_bytes());
    for _ in 0..3 {
        original
            .step(0.1, TransientStepControls::default(), || false)
            .unwrap();
        equivalent
            .step(0.1, TransientStepControls::default(), || false)
            .unwrap();
    }
    assert_eq!(snapshot(&original), snapshot(&equivalent));
    assert_eq!(equivalent.film_sorbed_solute_moles(0, 0, 0), Some(0.0));
    assert_eq!(
        equivalent.film_mobile_solute_moles(0, 0, 0),
        Some(equivalent.film_solute_moles(0, 0).unwrap()[0])
    );
    assert_eq!(equivalent.film_mobile_concentration_mol_m3(2, 0, 0), None);
    assert_eq!(equivalent.film_mobile_solute_moles(0, 2, 0), None);
    assert_eq!(equivalent.film_sorbed_solute_moles(0, 0, 1), None);
    let sorbing = prepare(definition.with_linear_sorption(vec![2.0, 10.0]).unwrap());
    assert!(sorbing.workspace_bytes() > original.workspace_bytes());

    let overflow = TransientFilm::new(1.0, 1, vec![1e-9], vec![1.0])
        .unwrap()
        .with_linear_sorption(vec![f64::MAX])
        .unwrap();
    assert!(
        TransientOsmoticPair::new(
            [bulk(&[0.0], 1e-6, 0.0), bulk(&[0.0], 1e-6, 0.0)],
            10.0,
            membrane(0.0, 0.0, 0.0, 1),
            Some(overflow),
            None,
            1 << 20,
        )
        .is_err()
    );
}

#[test]
fn heterogeneous_total_partition_preserves_mobile_equilibrium_and_osmotic_force() {
    let concentrations = [25.0, 5.0];
    let ratios = [
        vec![2.0, 7.0, 100.0, 1.0, 3.0, 4.0],
        vec![1.0, 30.0, 5.0, 2.0],
    ];
    let left = TransientFilm::from_cells(
        vec![1e-5, 3e-5, 6e-5],
        2,
        vec![1e-9, 2e-9, 3e-10, 5e-10, 8e-11, 1e-9],
        concentrations.repeat(3),
    )
    .unwrap()
    .with_linear_sorption(ratios[0].clone())
    .unwrap();
    let right = TransientFilm::from_cells(
        vec![8e-5, 2e-5],
        2,
        vec![1e-9, 2e-9, 4e-10, 3e-9],
        concentrations.repeat(2),
    )
    .unwrap()
    .with_linear_sorption(ratios[1].clone())
    .unwrap();
    let mut pair = TransientOsmoticPair::new(
        [
            bulk(&concentrations, 1e-6, 0.0),
            bulk(&concentrations, 1e-6, 0.0),
        ],
        1e-4,
        membrane(2e-11, 1.0, 1e-6, 2),
        Some(left),
        Some(right),
        1 << 20,
    )
    .unwrap();
    let before = [pair.total_amount(0), pair.total_amount(1)];
    let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
    let report = pair
        .advance(
            10.0,
            TransientAdvanceControls {
                maximum_step_s: 10.0,
                ..Default::default()
            },
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(
        report.solvent_transfer_m3, 0.0,
        "bound particles must not create an osmotic gradient"
    );
    close(pair.total_solvent_volume_m3(), 2.02e-6, 1e-15, 0.0);
    for species in 0..2 {
        close(pair.total_amount(species), before[species], 1e-12, 1e-19);
        close(
            report.endpoint.membrane_left_concentration_mol_m3[species],
            concentrations[species],
            1e-12,
            1e-12,
        );
        close(
            report.endpoint.membrane_right_concentration_mol_m3[species],
            concentrations[species],
            1e-12,
            1e-12,
        );
        for side in 0..2 {
            for cell in 0..pair.film_cell_count(side).unwrap() {
                let total = pair.film_solute_moles(side, cell).unwrap()[species];
                let mobile = pair.film_mobile_solute_moles(side, cell, species).unwrap();
                let bound = pair.film_sorbed_solute_moles(side, cell, species).unwrap();
                close(
                    pair.film_mobile_concentration_mol_m3(side, cell, species)
                        .unwrap(),
                    concentrations[species],
                    1e-12,
                    1e-12,
                );
                close(
                    total / pair.film_cell_volume_m3(side, cell).unwrap(),
                    ratios[side][2 * cell + species] * concentrations[species],
                    1e-12,
                    1e-12,
                );
                close(total, mobile + bound, 1e-15, 1e-24);
                assert!(bound >= 0.0);
            }
        }
    }
}

#[test]
fn finite_bath_sorption_matches_retarded_diffusion_and_dissipates_chemical_potential() {
    let area = 1e-4;
    let width = 1e-3;
    let diffusivity = 1e-9;
    let bulk_volume = 1e-6;
    let duration = 1000.0;
    let mut mobile_endpoints = Vec::new();
    for ratio in [1.0, 10.0] {
        let capacity = ratio * area * width;
        let conductance = area * diffusivity / (0.5 * width);
        let equilibrium = 50.0 * bulk_volume / (bulk_volume + capacity);
        let rate: f64 = conductance * (1.0 / bulk_volume + 1.0 / capacity);
        let expected = equilibrium * (1.0_f64 - (-rate * duration).exp());
        let left = TransientFilm::new(width, 1, vec![diffusivity], vec![0.0])
            .unwrap()
            .with_linear_sorption(vec![ratio])
            .unwrap();
        let mut pair = TransientOsmoticPair::new(
            [
                bulk(&[50.0], bulk_volume, 0.0),
                bulk(&[0.0], bulk_volume, 0.0),
            ],
            area,
            membrane(0.0, 1.0, 0.0, 1),
            Some(left),
            None,
            1 << 20,
        )
        .unwrap();
        // Reduced ideal sorption free energy / RT: dF/dn = ln(n/(R*V*c_ref)).
        // The arbitrary fixed reference c_ref=1 mol/m³ cancels in differences.
        let potential = |state: &TransientOsmoticPair| {
            (0..state.rows)
                .map(|row| {
                    let n = state.amount(row, 0);
                    if n == 0.0 {
                        0.0
                    } else {
                        n * ((n / state.solute_capacity(row, 0, 0.0)).ln() - 1.0)
                    }
                })
                .sum::<f64>()
        };
        let before_potential = potential(&pair);
        let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
        let report = pair
            .advance(
                duration,
                TransientAdvanceControls {
                    maximum_step_s: duration,
                    relative_tolerance: 1e-5,
                    ..Default::default()
                },
                &mut workspace,
                || false,
            )
            .unwrap();
        let mobile = pair.film_mobile_concentration_mol_m3(0, 0, 0).unwrap();
        close(mobile, expected, 1e-3, 1e-6);
        close(
            pair.film_solute_moles(0, 0).unwrap()[0],
            capacity * mobile,
            1e-14,
            1e-20,
        );
        close(pair.total_amount(0), 50.0 * bulk_volume, 1e-12, 1e-19);
        assert!(potential(&pair) < before_potential);
        assert!(report.maximum_scaled_error <= 1.0);
        assert!(report.solute_balance_residual_mol[0].abs() < 1e-16);
        mobile_endpoints.push(mobile);
    }
    assert!(mobile_endpoints[1] < mobile_endpoints[0] * 0.25);
}

#[test]
fn strong_advection_transports_mobile_solute_without_advecting_the_bound_pool() {
    for side in 0..2 {
        let film = TransientFilm::new(1e-4, 1, vec![1e-15], vec![100.0])
            .unwrap()
            .with_linear_sorption(vec![10.0])
            .unwrap();
        let mut films = [None, None];
        films[side] = Some(film);
        let [left, right] = films;
        let pressure = if side == 0 {
            [201_325.0, 101_325.0]
        } else {
            [101_325.0, 201_325.0]
        };
        let mut pair = TransientOsmoticPair::new(
            [
                bulk(&[0.0], 1e-6, pressure[0]),
                bulk(&[0.0], 1e-6, pressure[1]),
            ],
            1e-4,
            membrane(1e-11, 0.0, 0.0, 1),
            left,
            right,
            1 << 20,
        )
        .unwrap();
        // Pe=50,000: the original exponential law reaches its upwind limit.
        // Total inventory obeys dn/dt=-A*|Jv|*n/(R*V), not -A*|Jv|*n/V.
        let initial = pair.total_amount(0);
        let report = pair
            .step(100.0, TransientStepControls::default(), || false)
            .unwrap();
        let expected = initial / 1.1;
        close(report.volume_flux_m_s.abs(), 1e-6, 1e-12, 1e-18);
        close(
            pair.film_solute_moles(side, 0).unwrap()[0],
            expected,
            1e-12,
            1e-19,
        );
        close(
            pair.film_mobile_solute_moles(side, 0, 0).unwrap(),
            expected / 10.0,
            1e-12,
            1e-19,
        );
        close(
            pair.film_sorbed_solute_moles(side, 0, 0).unwrap(),
            0.9 * expected,
            1e-12,
            1e-19,
        );
        close(
            pair.bulks()[1 - side].solution().solute_moles()[0],
            initial - expected,
            1e-12,
            1e-19,
        );
        close(pair.total_amount(0), initial, 1e-12, 1e-19);
    }
}
