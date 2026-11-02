//! Passive zero-storage film contacts, including both exterior ends.
use super::*;

#[test]
fn zero_resistance_is_the_original_path_and_admission_is_explicit() {
    for bad in [
        vec![0.0; 2],
        vec![0.0, -1.0, 0.0],
        vec![f64::NAN; 3],
        vec![f64::INFINITY; 3],
    ] {
        assert!(film(2, 20.0).with_interface_resistances(bad).is_err());
    }
    let make = |resistance: Option<Vec<f64>>, budget| {
        let f = film(2, 20.0);
        let f = match resistance {
            Some(r) => f.with_interface_resistances(r).unwrap(),
            None => f,
        };
        TransientOsmoticPair::new(
            [
                bulk(&[50.0], 1e-6, 101_325.0),
                bulk(&[100.0], 1e-6, 101_325.0),
            ],
            1e-4,
            membrane(2e-11, 0.5, 1e-5, 1),
            Some(f),
            Some(film(2, 100.0)),
            budget,
        )
    };
    let mut original = make(None, 1 << 20).unwrap();
    let mut zero = make(Some(vec![0.0; 3]), 1 << 20).unwrap();
    assert_eq!(original.workspace_bytes(), zero.workspace_bytes());
    for _ in 0..5 {
        original
            .step(0.03, TransientStepControls::default(), || false)
            .unwrap();
        zero.step(0.03, TransientStepControls::default(), || false)
            .unwrap();
        assert_eq!(snapshot(&original), snapshot(&zero));
    }
    let finite = make(Some(vec![1.0, 0.0, 2.0]), 1 << 20).unwrap();
    assert_eq!(
        finite.workspace_bytes(),
        original.workspace_bytes() + 3 * size_of::<f64>()
    );
    assert!(matches!(
        make(Some(vec![1.0, 0.0, 2.0]), finite.workspace_bytes() - 1),
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::Workspace { .. }
        ))
    ));
    // Pure solvent needs no contact species coefficients.
    assert!(
        TransientFilm::new(1.0, 2, vec![], vec![])
            .unwrap()
            .with_interface_resistances(vec![])
            .is_ok()
    );
}

#[test]
fn finite_bath_uptake_uses_each_exterior_species_resistance_and_mobile_concentration() {
    let dt = 3.0;
    let area = 2e-4;
    let width = 4e-5;
    let volume = 1e-7;
    let diffusivity = [1e-9, 4e-10];
    let resistance = [2e5, 8e5];
    let ratios = [2.0, 5.0];
    for side in 0..2 {
        let mut contacts = vec![0.0; 4];
        let outer = if side == 0 { 0 } else { 1 };
        contacts[2 * outer..2 * outer + 2].copy_from_slice(&resistance);
        // This opposite contact touches the sealed membrane and must add no
        // fictitious leak, despite a different resistance for every species.
        contacts[2 * (1 - outer)..2 * (1 - outer) + 2].copy_from_slice(&[5e4, 1e4]);
        let f = TransientFilm::new(width, 1, diffusivity.to_vec(), vec![0.0; 2])
            .unwrap()
            .with_linear_sorption(ratios.to_vec())
            .unwrap()
            .with_interface_resistances(contacts)
            .unwrap();
        let mut definitions = [None, None];
        definitions[side] = Some(f);
        let [left, right] = definitions;
        let mut pair = TransientOsmoticPair::new(
            [
                bulk(&[100.0, 40.0], volume, 0.0),
                bulk(&[100.0, 40.0], volume, 0.0),
            ],
            area,
            membrane(0.0, 1.0, 0.0, 2),
            left,
            right,
            1 << 20,
        )
        .unwrap();
        let report = pair
            .step(dt, TransientStepControls::default(), || false)
            .unwrap();
        for species in 0..2 {
            let initial = [100.0, 40.0][species];
            let capacity = ratios[species] * area * width;
            let conductance = area / (0.5 * width / diffusivity[species] + resistance[species]);
            let difference = initial / (1.0 + dt * conductance * (1.0 / volume + 1.0 / capacity));
            let expected_bulk = (initial * volume + capacity * difference) / (volume + capacity);
            let expected_cell = expected_bulk - difference;
            close(
                pair.bulks()[side].solution().solute_moles()[species] / volume,
                expected_bulk,
                2e-14,
                1e-13,
            );
            close(
                pair.film_mobile_concentration_mol_m3(side, 0, species)
                    .unwrap(),
                expected_cell,
                3e-13,
                1e-13,
            );
            close(
                pair.bulks()[1 - side].solution().solute_moles()[species] / volume,
                initial,
                1e-15,
                1e-13,
            );
            assert!(report.solute_balance_residual_mol[species].abs() < 1e-18);
        }
    }
}

#[test]
fn contact_coefficients_retain_large_signed_peclet_and_extreme_resistance_limits() {
    let definition =
        TransientFilm::from_cells(vec![2e-5, 6e-5], 1, vec![2e-9, 5e-10], vec![10.0, 20.0])
            .unwrap()
            .with_interface_resistances(vec![0.0, 3e5, 0.0])
            .unwrap();
    let resistance = 1e-5 / 2e-9 + 3e-5 / 5e-10 + 3e5;
    let operator = film_operator(&definition, 0, Some(1), 1).unwrap();
    for peclet in [-1000.0_f64, -10.0, -0.001, 0.0, 0.001, 10.0, 1000.0] {
        let speed = peclet / resistance;
        let (left, right) = operator.solute_stream_coefficients(0, speed).unwrap();
        let p = peclet.abs();
        let down = if p == 0.0 {
            1.0 / resistance
        } else {
            speed.abs() / p.exp_m1()
        };
        let up = down + speed.abs();
        let expected = if peclet >= 0.0 {
            (up, down)
        } else {
            (down, up)
        };
        close(left, expected.0, 2e-13, 1e-300);
        close(right, expected.1, 2e-13, 1e-300);
        assert!(left >= 0.0 && right >= 0.0);
        close(left - right, speed, 2e-12, 1e-20);
    }
    // Product k*r overflows, while the physical series conductance is finite.
    let extreme = TransientFilm::new(2.0, 1, vec![1e200], vec![1.0])
        .unwrap()
        .with_interface_resistances(vec![1e200, 0.0])
        .unwrap();
    let operator = film_operator(&extreme, 0, None, 0).unwrap();
    let (a, b) = operator.solute_stream_coefficients(0, 0.0).unwrap();
    close(a, 1e-200, 1e-15, 0.0);
    assert_eq!(a, b);
}

#[test]
fn whole_film_profiles_include_every_contact_and_the_membrane_faces_in_both_directions() {
    let widths = [2e-5, 3e-5];
    let d = [2e-9, 5e-10];
    let contacts = [[3e4, 8e4, 2e4], [5e4, 4e4, 6e4]];
    let resistance = contacts.map(|c| widths[0] / d[0] + widths[1] / d[1] + c.iter().sum::<f64>());
    let membrane_permeance = 2e-5;
    for speed in [-5e-6_f64, 0.0, 5e-6] {
        let mut pair = TransientOsmoticPair::new(
            [bulk(&[80.0], 1.0, speed / 1e-10), bulk(&[20.0], 1.0, 0.0)],
            1e-4,
            membrane(1e-10, 0.0, membrane_permeance, 1),
            Some(
                TransientFilm::from_cells(widths.to_vec(), 1, d.to_vec(), vec![80.0; 2])
                    .unwrap()
                    .with_interface_resistances(contacts[0].to_vec())
                    .unwrap(),
            ),
            Some(
                TransientFilm::from_cells(widths.to_vec(), 1, d.to_vec(), vec![20.0; 2])
                    .unwrap()
                    .with_interface_resistances(contacts[1].to_vec())
                    .unwrap(),
            ),
            1 << 20,
        )
        .unwrap();
        let mut report = pair
            .step(100.0, TransientStepControls::default(), || false)
            .unwrap();
        for _ in 0..29 {
            report = pair
                .step(100.0, TransientStepControls::default(), || false)
                .unwrap();
        }
        let concentration = |side: usize| {
            let b = pair.bulks()[side].solution();
            b.solute_moles()[0] / b.solvent_volume_m3()
        };
        let cl = concentration(0);
        let cr = concentration(1);
        let total = resistance[0] + 1.0 / membrane_permeance + resistance[1];
        let flux = if speed == 0.0 {
            (cl - cr) / total
        } else {
            speed * (cl * (speed * total).exp() - cr) / (speed * total).exp_m1()
        };
        let at = |r: f64| {
            if speed == 0.0 {
                cl - flux * r
            } else {
                cl * (speed * r).exp() - flux / speed * (speed * r).exp_m1()
            }
        };
        close(report.volume_flux_m_s, speed, 1e-12, 1e-16);
        close(
            report.membrane_left_concentration_mol_m3[0],
            at(resistance[0]),
            5e-7,
            1e-7,
        );
        close(
            report.membrane_right_concentration_mol_m3[0],
            at(resistance[0] + 1.0 / membrane_permeance),
            5e-7,
            1e-7,
        );
        // Original cell storage remains at its physical center, before/after
        // the zero-storage face jump; contact resistance is never cell width.
        let mut cumulative = contacts[0][0];
        for cell in 0..2 {
            cumulative += 0.5 * widths[cell] / d[cell];
            close(
                pair.film_mobile_concentration_mol_m3(0, cell, 0).unwrap(),
                at(cumulative),
                5e-7,
                1e-7,
            );
            cumulative += 0.5 * widths[cell] / d[cell] + contacts[0][cell + 1];
        }
        assert!(report.solute_balance_residual_mol[0].abs() < 1e-12);
    }
}

fn sorbing_contact_pair() -> TransientOsmoticPair {
    let f = TransientFilm::from_cells(vec![2e-5, 3e-5], 1, vec![1e-9, 3e-10], vec![0.0; 2])
        .unwrap()
        .with_interface_resistances(vec![4e5, 8e5, 2e5])
        .unwrap()
        .with_linear_sorption(vec![2.0, 3.0])
        .unwrap()
        .with_langmuir_sorption(vec![
            Some(
                FilmLangmuirSite::new(1e-6, 1e-3, 0.01, 0.0).unwrap()
            );
            2
        ])
        .unwrap();
    TransientOsmoticPair::new(
        [bulk(&[100.0], 1e-6, 0.0), bulk(&[0.0], 1e-6, 0.0)],
        1e-4,
        membrane(0.0, 0.0, 1e-6, 1),
        Some(f),
        None,
        1 << 20,
    )
    .unwrap()
}

#[test]
fn adaptive_contact_limited_uptake_evolves_and_conserves_both_sorption_pools() {
    let mut reference = sorbing_contact_pair();
    for _ in 0..512 {
        reference
            .step(0.25 / 512.0, TransientStepControls::default(), || false)
            .unwrap();
    }
    let mut pair = sorbing_contact_pair();
    let initial = pair.total_amount(0);
    let mut work = pair.adaptive_workspace(1 << 20).unwrap();
    let report = pair
        .advance(
            0.25,
            TransientAdvanceControls {
                maximum_step_s: 0.25,
                relative_tolerance: 1e-5,
                ..Default::default()
            },
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.accepted_steps > 1);
    assert_eq!(pair.time_s(), 0.25);
    close(pair.total_amount(0), initial, 3e-13, 1e-18);
    for cell in 0..2 {
        let mobile = pair.film_mobile_concentration_mol_m3(0, cell, 0).unwrap();
        let trapped = pair.film_trapped_solute_moles(0, cell, 0).unwrap();
        assert!(mobile > 0.0 && trapped > 0.0);
        close(
            mobile,
            reference
                .film_mobile_concentration_mol_m3(0, cell, 0)
                .unwrap(),
            2e-3,
            1e-7,
        );
        close(
            trapped,
            reference.film_trapped_solute_moles(0, cell, 0).unwrap(),
            2e-3,
            1e-15,
        );
    }
}
