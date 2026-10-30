use super::*;

fn properties(density: f64, mixing: ConductivityMixing) -> CharProperties<'static> {
    CharProperties {
        density_kg_m3: density,
        conductivity_w_m_k: 0.05,
        conductivity_mixing: mixing,
        source: "manufactured incompressible char; no wood material claim",
    }
}

fn totals(slab: &RecedingSlab<'_>) -> [f64; 2] {
    slab.state()
        .inventories()
        .values()
        .fold([0.0; 2], |mut sum, cell| {
            sum[0] += cell.virgin_mass_kg;
            sum[1] += cell.char_mass_kg;
            sum
        })
}

#[test]
fn retained_char_sets_independent_volume_and_both_conductivity_limits() {
    let gas = |t| CALORIC.specific_enthalpy_j_kg(t);
    for mixing in [ConductivityMixing::Parallel, ConductivityMixing::Series] {
        let initial = fixture_with_components(
            &CALORIC,
            &gas,
            &[600.0; 8],
            &[1.0; 8],
            2.0,
            Some((&CALORIC, 0.2, properties(200.0, mixing))),
        );
        let mut work = initial.state().implicit_workspace().unwrap();
        let (next, audit) = initial
            .advance_implicit_to(
                0.5,
                PyrolysisControls {
                    maximum_dt_s: 0.1,
                    ..Default::default()
                },
                ImplicitSolve::default(),
                SlabSolve::default(),
                &mut work,
                || Ok(()),
            )
            .unwrap();
        let virgin = initial.mass_kg() / 1.2_f64.powi(5);
        let char = 0.2 * (initial.mass_kg() - virgin);
        let volume = virgin / DENSITY + char / 200.0;
        let observed = totals(&next);
        approx(observed[0], virgin, 1e-16);
        approx(observed[1], char, 1e-16);
        approx(next.mass_kg(), virgin + char, 1e-16);
        approx(next.thickness_m(), volume / AREA, 1e-16);
        approx(
            audit.pyrolysis.volatile_mass_kg,
            0.8 * (initial.mass_kg() - virgin),
            1e-16,
        );
        approx(
            next.energy_j() + audit.pyrolysis.volatile_enthalpy_j,
            initial.energy_j(),
            1e-8,
        );
        assert_eq!(next.back_m().to_bits(), initial.back_m().to_bits());
        let phi = char / (200.0 * volume);
        let expected_k = match mixing {
            ConductivityMixing::Parallel => (1.0 - phi) * CONDUCTIVITY + phi * 0.05,
            ConductivityMixing::Series => 1.0 / ((1.0 - phi) / CONDUCTIVITY + phi / 0.05),
        };
        for &id in &next.ids {
            approx(next.cell_conductivity_w_m_k(id).unwrap(), expected_k, 1e-14);
            approx(next.state().temperature_k(id).unwrap(), 600.0, 1e-9);
        }
        approx(
            next.front_half_cell_conductance_w_k().unwrap(),
            2.0 * expected_k * AREA / (next.thickness_m() / 8.0),
            1e-11,
        );
        assert!(audit.absolute_remap_component_mass_residual_kg < 1e-15);
        assert!(audit.remap_char_mass_residual_kg.abs() < 1e-16);
        // A virgin-density-only geometry would lose most of the retained volume.
        assert!((next.thickness_m() - next.mass_kg() / (DENSITY * AREA)).abs() > 1e-3);
    }
}

fn set_composition(slab: &mut RecedingSlab<'_>, id: CellId, char_fraction: f64, mass_scale: f64) {
    let cell = slab.solid.thermal.cells.get_mut(&id).unwrap();
    let inventory = slab.solid.inventories.get_mut(&id).unwrap();
    cell.mass_kg *= mass_scale;
    inventory.char_mass_kg = cell.mass_kg * char_fraction;
    inventory.virgin_mass_kg = cell.mass_kg - inventory.char_mass_kg;
    let h = slab.solid.materials[&cell.material_id]
        .enthalpies(inventory.temperature_k)
        .unwrap();
    cell.specific_h_j_kg = h[0] + char_fraction * (h[1] - h[0]);
}

#[test]
fn mixed_nonlinear_caloric_overlap_preserves_each_component_and_actual_energy() {
    let virgin = Caloric {
        cp: 2000.0,
        offset: -4e6,
        quadratic: 0.4,
    };
    let char = Caloric {
        cp: 1100.0,
        offset: -2e6,
        quadratic: 1.2,
    };
    let gas = |t| Ok(3e6 + 2500.0 * (t - 300.0));
    let mut slab = fixture_with_components(
        &virgin,
        &gas,
        &[700.0, 350.0],
        &[1.0, 1.0],
        0.0,
        Some((&char, 0.2, properties(200.0, ConductivityMixing::Series))),
    );
    set_composition(&mut slab, CellId(0), 0.2, 0.5);
    set_composition(&mut slab, CellId(1), 0.8, 1.0);
    let old_h = slab.state().cells()[&CellId(0)].specific_h_j_kg;
    let cold_h = slab.state().cells()[&CellId(1)].specific_h_j_kg;
    let expected_h = (2.0 * old_h + cold_h) / 3.0;
    let before = totals(&slab);
    let energy = slab.energy_j();
    let report = slab.rezone(SlabSolve::default(), &mut || Ok(())).unwrap();
    let after = totals(&slab);
    approx(before[0], after[0], 1e-16);
    approx(before[1], after[1], 1e-16);
    approx(slab.energy_j(), energy, 1e-8);
    let first = slab.state().inventories()[&CellId(0)];
    approx(
        first.char_mass_kg / slab.state().cells()[&CellId(0)].mass_kg,
        0.4,
        1e-15,
    );
    approx(
        slab.state().cells()[&CellId(0)].specific_h_j_kg,
        expected_h,
        1e-7,
    );
    let recovered = 0.6 * virgin.specific_enthalpy_j_kg(first.temperature_k).unwrap()
        + 0.4 * char.specific_enthalpy_j_kg(first.temperature_k).unwrap();
    approx(recovered, expected_h, 1e-7);
    assert!((first.temperature_k - (2.0 * 700.0 + 350.0) / 3.0).abs() > 10.0);
    approx(
        slab.thickness_m(),
        (after[0] / DENSITY + after[1] / 200.0) / AREA,
        1e-16,
    );
    for pair in slab.ids.windows(2) {
        let left = &slab.state().cells()[&pair[0]];
        let right = &slab.state().cells()[&pair[1]];
        let expected = series_conductance_w_k(
            AREA,
            [
                0.5 * (left.geometry.end_m[0] - left.geometry.start_m[0]),
                0.5 * (right.geometry.end_m[0] - right.geometry.start_m[0]),
            ],
            [
                slab.cell_conductivity_w_m_k(pair[0]).unwrap(),
                slab.cell_conductivity_w_m_k(pair[1]).unwrap(),
            ],
        )
        .unwrap();
        approx(
            slab.solid.thermal.links[&1].conductance_w_k,
            expected,
            1e-14,
        );
    }
    assert!(report.energy_j.abs() < 1e-8 && report.char_mass_kg.abs() < 1e-16);
}

#[test]
fn char_reconstruction_refines_without_phase_mass_or_energy_datum_loss() {
    let run = |n: usize, shifts: [f64; 2]| {
        let virgin = Caloric {
            cp: 2000.0,
            offset: -4e6 + shifts[0],
            quadratic: 0.4,
        };
        let char = Caloric {
            cp: 1100.0,
            offset: -2e6 + shifts[1],
            quadratic: 1.2,
        };
        let gas_shift = (shifts[0] - 0.2 * shifts[1]) / 0.8;
        let gas = |t| Ok(3e6 + gas_shift + 2500.0 * (t - 300.0));
        let weights: Vec<_> = (0..n).map(|i| if i % 2 == 0 { 1.0 } else { 2.0 }).collect();
        let mut slab = fixture_with_components(
            &virgin,
            &gas,
            &vec![600.0; n],
            &weights,
            0.0,
            Some((&char, 0.2, properties(200.0, ConductivityMixing::Parallel))),
        );
        let total: f64 = weights.iter().sum();
        let mut a = 0.0;
        for (i, w) in weights.iter().enumerate() {
            let b = a + w / total;
            let fraction = 0.1 + 0.5 * (a * a + a * b + b * b) / 3.0;
            set_composition(&mut slab, CellId(i as u64), fraction, 1.0);
            a = b;
        }
        let before = totals(&slab);
        let energy = slab.energy_j();
        slab.rezone(SlabSolve::default(), &mut || Ok(())).unwrap();
        let after = totals(&slab);
        approx(after[0], before[0], 1e-15);
        approx(after[1], before[1], 1e-15);
        approx(slab.energy_j(), energy, 1e-7);
        let mut error = 0.0;
        let mut temperatures = Vec::new();
        for (i, &id) in slab.ids.iter().enumerate() {
            let a = i as f64 / n as f64;
            let b = (i + 1) as f64 / n as f64;
            let expected = 0.1 + 0.5 * (a * a + a * b + b * b) / 3.0;
            let inventory = slab.solid.inventories[&id];
            let actual = inventory.char_mass_kg / slab.solid.cells()[&id].mass_kg;
            error += (actual - expected).abs() / n as f64;
            temperatures.push(inventory.temperature_k);
            approx(inventory.temperature_k, 600.0, 2e-8);
        }
        (error, temperatures)
    };
    let errors = [8, 16, 32].map(|n| run(n, [0.0, 0.0]).0);
    assert!(
        errors[0] / errors[1] > 3.0 && errors[1] / errors[2] > 3.0,
        "{errors:?}"
    );
    let base = run(16, [0.0, 0.0]).1;
    let shifted = run(16, [1e7, -2e7]).1;
    for (a, b) in base.into_iter().zip(shifted) {
        approx(a, b, 2e-8);
    }
}

#[test]
fn swelling_retries_on_component_volume_even_when_mass_loss_is_small() {
    let gas = |t| CALORIC.specific_enthalpy_j_kg(t);
    let slab = fixture_with_components(
        &CALORIC,
        &gas,
        &[600.0; 3],
        &[1.0; 3],
        5.0,
        Some((
            &CALORIC,
            0.8,
            properties(100.0, ConductivityMixing::Parallel),
        )),
    );
    let controls = PyrolysisControls {
        maximum_dt_s: 0.05,
        ..Default::default()
    };
    let geometry = SlabSolve {
        maximum_cell_volume_change_fraction: 0.1,
        ..Default::default()
    };
    let mut work = slab.state().implicit_workspace().unwrap();
    let refused = slab.advance_implicit_to(
        0.05,
        controls,
        ImplicitSolve::default(),
        SlabSolve {
            maximum_halvings: 0,
            ..geometry
        },
        &mut work,
        || Ok(()),
    );
    assert!(matches!(refused, Err(SlabError::GeometrySteps)));
    let (next, audit) = slab
        .advance_implicit_to(
            0.05,
            controls,
            ImplicitSolve::default(),
            geometry,
            &mut work,
            || Ok(()),
        )
        .unwrap();
    assert!(audit.rejected_geometry_attempts > 0);
    assert!(audit.maximum_cell_volume_change_fraction <= 0.1);
    assert!(audit.maximum_cell_mass_loss_fraction < 0.02);
    assert!(next.recession_m() < 0.0 && next.front_m() < slab.front_m());
    assert_eq!(next.back_m().to_bits(), slab.back_m().to_bits());
    approx(
        next.mass_kg() + audit.pyrolysis.volatile_mass_kg,
        slab.mass_kg(),
        1e-15,
    );
    approx(
        next.energy_j() + audit.pyrolysis.volatile_enthalpy_j,
        slab.energy_j(),
        1e-8,
    );
    let masses = totals(&next);
    approx(
        next.thickness_m(),
        (masses[0] / DENSITY + masses[1] / 100.0) / AREA,
        1e-16,
    );
}

#[test]
fn char_source_refusals_and_late_cancellation_publish_no_partial_state() {
    let gas = |t| CALORIC.specific_enthalpy_j_kg(t);
    let slab = fixture_with_components(
        &CALORIC,
        &gas,
        &[650.0, 500.0, 300.0],
        &[1.0; 3],
        0.4,
        Some((&CALORIC, 0.2, properties(200.0, ConductivityMixing::Series))),
    );
    for bad in [
        CharProperties {
            density_kg_m3: 0.0,
            ..slab.char_properties().unwrap()
        },
        CharProperties {
            conductivity_w_m_k: f64::NAN,
            ..slab.char_properties().unwrap()
        },
        CharProperties {
            source: " ",
            ..slab.char_properties().unwrap()
        },
    ] {
        assert!(
            RecedingSlab::with_char(slab.solid.clone(), slab.geometry(), bad, || Ok(())).is_err()
        );
    }
    let mut work = slab.state().implicit_workspace().unwrap();
    let mut calls = 0usize;
    let reference = slab
        .advance_implicit_to(
            0.02,
            PyrolysisControls::default(),
            ImplicitSolve::default(),
            SlabSolve::default(),
            &mut work,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    let cells = slab.state().cells().clone();
    let inventories = slab.state().inventories().clone();
    let mut count = 0usize;
    let cancelled = slab.advance_implicit_to(
        0.02,
        PyrolysisControls::default(),
        ImplicitSolve::default(),
        SlabSolve::default(),
        &mut work,
        || {
            count += 1;
            if count >= calls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(cancelled, Err(SlabError::Cancelled)));
    assert_eq!(slab.state().cells(), &cells);
    assert_eq!(slab.state().inventories(), &inventories);
    let replay = slab
        .advance_implicit_to(
            0.02,
            PyrolysisControls::default(),
            ImplicitSolve::default(),
            SlabSolve::default(),
            &mut work,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(reference.0.state().cells(), replay.0.state().cells());
    assert_eq!(
        reference.0.state().inventories(),
        replay.0.state().inventories()
    );
    assert_eq!(reference.1, replay.1);
}
