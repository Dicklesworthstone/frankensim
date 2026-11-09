use super::super::cells::partition;
use super::*;
use crate::radiation::spectral::planck::BlackbodyBand;

mod irradiation;

fn cell_state<'a>(
    law: &'a Caloric,
    response: &SlabCellResponse,
    energies: &[f64],
    area: f64,
) -> DepositionSimulation<'a> {
    let events = response
        .cells()
        .iter()
        .enumerate()
        .map(|(i, cell)| {
            let volume = area * (cell.interval_m[1] - cell.interval_m[0]);
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(i as u64),
                    segment_id: cell.layer as u64,
                    body_id: 7,
                    material_id: MaterialId(7),
                    mass_kg: 300.0 * volume,
                    specific_h_j_kg: energies[i],
                    geometry: CellGeometry {
                        start_m: [cell.interval_m[0], 0.0, 0.0],
                        end_m: [cell.interval_m[1], 0.0, 0.0],
                        volume_m3: volume,
                        surface_area_m2: 6.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(7), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 8,
            max_substeps: 100_000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn g1_nonuniform_cells_exchange_native_heat_inside_one_original_material() {
    let original = asymmetric::between(&[3000.0], &[1.5], [1.0, 1.3]);
    let response = partition(&original, 1.0 / 3.0);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = cell_state(&law, &response, &[900_000.0, 300_000.0], 0.01);
    let indices = original.exterior_refractive_indices();
    let irradiation = [
        indices[0].powi(2) * exitance(900.0),
        indices[1].powi(2) * exitance(300.0),
    ];
    let source = GraySlabThermal::new_cells_full_spectrum(
        &response,
        0.01,
        irradiation,
        &bindings(&law, 2),
        &old,
        || Ok(()),
    )
    .unwrap();
    let transport = source.endpoint_transport(&old, || Ok(())).unwrap();
    assert!(transport.heat.bulk.sources_w.is_empty());
    assert_eq!(transport.transfers.len(), 1);
    let power = 0.01 * 2.25 * (exitance(900.0) - exitance(300.0));
    close(transport.transfers[0].power_w, power, 2e-10);
    let (next, _, audit) = old
        .advance_observed(&source, 1e-4, config(), || Ok(()))
        .unwrap();
    for (i, sign) in [(0, -1.0), (1, 1.0)] {
        let id = CellId(i);
        close(
            next.cells()[&id].specific_h_j_kg,
            old.cells()[&id].specific_h_j_kg + sign * 1e-4 * power / old.cells()[&id].mass_kg,
            2e-10,
        );
        assert_eq!(next.cells()[&id].geometry, old.cells()[&id].geometry);
    }
    assert_eq!(audit.external_heat_j, 0.0);
    assert!(next.links().is_empty());
    assert_eq!(source.slab().layer_count(), 1);
    assert_eq!(source.cell_response().unwrap().cell_count(), 2);
    assert_eq!(source.slab().layers().as_ptr(), original.layers().as_ptr());
}

#[test]
fn g1_spectral_original_layer_cells_drive_native_implicit_energy_without_new_interfaces() {
    let originals = [
        asymmetric::between(&[3000.0], &[1.5], [1.2, 1.7]),
        asymmetric::between(&[3000.0], &[2.0], [1.5, 2.2]),
    ];
    let responses = originals.each_ref().map(|s| partition(s, 1.0 / 3.0));
    let windows = [
        BlackbodyBand::new(1e-6, 5e-6).unwrap(),
        BlackbodyBand::new(8e-6, 20e-6).unwrap(),
    ];
    let power: f64 = (0..2)
        .map(|b| {
            0.01 * originals[b].refractive_index().unwrap().powi(2)
                * windows[b].difference(550.0, 350.0).unwrap()
        })
        .sum();
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let known = cell_state(&law, &responses[0], &[550_000.0, 350_000.0], 0.01);
    let old = cell_state(
        &law,
        &responses[0],
        &[
            550_000.0 + power / known.cells()[&CellId(0)].mass_kg,
            350_000.0 - power / known.cells()[&CellId(1)].mass_kg,
        ],
        0.01,
    );
    let bands = (0..2)
        .map(|b| SlabThermalCellBand {
            window: windows[b],
            response: responses[b].clone(),
            irradiation_w_m2: core::array::from_fn(|side| {
                originals[b].exterior_refractive_indices()[side].powi(2)
                    * windows[b]
                        .evaluate([550.0, 350.0][side])
                        .unwrap()
                        .power_w_m2
            }),
        })
        .collect();
    let source =
        SpectralSlabThermal::new_cells(bands, 0.01, &bindings(&law, 2), &old, || Ok(())).unwrap();
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_transport(
            &[-power, power],
            |candidate| source.endpoint_transport(candidate, || Ok(())),
            || Ok(()),
        )
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-8);
    let (next, (), audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
    close(next.temperature_k(CellId(0)).unwrap(), 550.0, 1e-12);
    close(next.temperature_k(CellId(1)).unwrap(), 350.0, 1e-12);
    close(audit.external_heat_j, 0.0, 1e-10);
    assert_eq!(source.bands()[0].optics.layer_count(), 1);
    assert_eq!(
        source.cell_response(0).unwrap().cells(),
        source.cell_response(1).unwrap().cells()
    );
    assert_eq!(
        source
            .endpoint_transport(&known, || Ok(()))
            .unwrap()
            .transfers
            .len(),
        1
    );
}

#[test]
fn g3_gray_and_spectral_native_isothermal_cell_sums_recover_original_layer_heat() {
    let original = asymmetric::between(&[0.8], &[1.5], [1.2, 1.7]);
    let whole = original
        .cell_response(
            vec![SlabCell {
                layer: 0,
                interval_m: [0.0, 1.0],
            }],
            SlabQuadrature::default(),
            || false,
        )
        .unwrap();
    let response = partition(&original, 0.4);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let coarse = cell_state(&law, &whole, &[400_000.0], 0.01);
    let fine = cell_state(&law, &response, &[400_000.0; 2], 0.01);
    let index = original.exterior_refractive_indices();
    let irradiation = [
        index[0].powi(2) * exitance(600.0),
        index[1].powi(2) * exitance(300.0),
    ];
    let a = GraySlabThermal::new_full_spectrum(
        &original,
        0.01,
        irradiation,
        &bindings(&law, 1),
        &coarse,
        || Ok(()),
    )
    .unwrap();
    let b = GraySlabThermal::new_cells_full_spectrum(
        &response,
        0.01,
        irradiation,
        &bindings(&law, 2),
        &fine,
        || Ok(()),
    )
    .unwrap();
    let pa = a.endpoint_transport(&coarse, || Ok(())).unwrap();
    let pb = b.endpoint_transport(&fine, || Ok(())).unwrap();
    close(
        pa.heat.bulk.sources_w.values().sum(),
        pb.heat.bulk.sources_w.values().sum(),
        2e-8,
    );
    assert!(pb.transfers.is_empty());
    let (_, _, aa) = coarse
        .advance_observed(&a, 1e-4, config(), || Ok(()))
        .unwrap();
    let (_, _, ab) = fine
        .advance_observed(&b, 1e-4, config(), || Ok(()))
        .unwrap();
    close(aa.external_heat_j, ab.external_heat_j, 2e-12);
    let window = BlackbodyBand::new(3e-6, 14e-6).unwrap();
    let irradiation = [
        index[0].powi(2) * window.evaluate(600.0).unwrap().power_w_m2,
        index[1].powi(2) * window.evaluate(300.0).unwrap().power_w_m2,
    ];
    let a = SpectralSlabThermal::new(
        vec![SlabThermalBand {
            window,
            optics: original.clone(),
            irradiation_w_m2: irradiation,
        }],
        0.01,
        &bindings(&law, 1),
        &coarse,
        || Ok(()),
    )
    .unwrap();
    let b = SpectralSlabThermal::new_cells(
        vec![SlabThermalCellBand {
            window,
            response: response.clone(),
            irradiation_w_m2: irradiation,
        }],
        0.01,
        &bindings(&law, 2),
        &fine,
        || Ok(()),
    )
    .unwrap();
    let pa = a.endpoint_transport(&coarse, || Ok(())).unwrap();
    let pb = b.endpoint_transport(&fine, || Ok(())).unwrap();
    close(
        pa.heat.bulk.sources_w.values().sum(),
        pb.heat.bulk.sources_w.values().sum(),
        2e-8,
    );
    assert!(pb.transfers.is_empty());
    let (_, _, aa) = coarse
        .advance_observed(&a, 1e-4, config(), || Ok(()))
        .unwrap();
    let (_, _, ab) = fine
        .advance_observed(&b, 1e-4, config(), || Ok(()))
        .unwrap();
    close(aa.external_heat_j, ab.external_heat_j, 2e-12);
}

#[test]
fn g2_native_depth_cells_keep_original_geometry_support_and_common_spectral_partition() {
    let original = asymmetric::between(&[1.0], &[1.5], [1.2, 1.7]);
    let response = partition(&original, 0.5);
    let changed = partition(&original, 0.4);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = cell_state(&law, &response, &[400_000.0; 2], 0.01);
    assert!(
        GraySlabThermal::new_cells_full_spectrum(
            &changed,
            0.01,
            [0.0; 2],
            &bindings(&law, 2),
            &old,
            || Ok(())
        )
        .is_err()
    );
    let unsupported = cell_state(&law, &response, &[400_000.0, 3_500_000.0], 0.01);
    assert!(matches!(
        GraySlabThermal::new_cells_full_spectrum(
            &response,
            0.01,
            [0.0; 2],
            &bindings(&law, 2),
            &unsupported,
            || Ok(())
        ),
        Err(DepositionError::MaterialDomain)
    ));
    let bands = vec![
        SlabThermalCellBand {
            window: BlackbodyBand::new(1e-6, 4e-6).unwrap(),
            response,
            irradiation_w_m2: [0.0; 2],
        },
        SlabThermalCellBand {
            window: BlackbodyBand::new(4e-6, 8e-6).unwrap(),
            response: changed,
            irradiation_w_m2: [0.0; 2],
        },
    ];
    assert!(
        SpectralSlabThermal::new_cells(bands, 0.01, &bindings(&law, 2), &old, || Ok(())).is_err()
    );
    assert_eq!(old.time_s(), 0.0);
}
