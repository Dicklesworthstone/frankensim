//! G1/G3/G4: original-material shielding, daughter heat and binding rollback.
//! Mass kernels and the two-normal-direction quadrature are manufactured;
//! the evaluated-data, decay, photon and thermal owners are production paths.

use super::*;
use crate::decay::thermal::DecayCellSource;
use crate::decay::{DecayBudget, DecayInventory};
use crate::endf::decay::model::emissions::{DecayEmissionLimits, PreparedDecayEmission};
use crate::endf::decay::spectra::DecayRadiation;
use crate::radiation::decay::{DecayRadiationClosure, PreparedDecayRadiation};
use crate::radiation::transport::material::{PhotonMassCollision, PhotonMassMaterial};
use crate::radiation::transport::refinement::PhotonSlabRefinement;
use crate::radiation::transport::{PhotonOrdinate, PhotonSlabLimits, PhotonSlabResponse};
use fs_la::decay::DecayLimits;

const EV_J: f64 = 1.602_176_634e-19;
const INITIAL: f64 = 1e14;

fn edges() -> Vec<f64> {
    vec![0.0, 1.5e6 * EV_J, 3e6 * EV_J]
}

fn partition() -> Arc<PreparedDecayRadiation> {
    let assembly = crate::radiation::decay::tests::gamma_chain(1e6, 2e6, 1.0, 1.0);
    Arc::new(
        PreparedDecayRadiation::new(
            Arc::new(
                PreparedDecayEmission::new(
                    assembly,
                    DecayRadiation::Gamma,
                    edges(),
                    DecayEmissionLimits::default(),
                    || Ok(()),
                )
                .unwrap(),
            ),
            DecayRadiationClosure {
                absolute_j_per_decay: 0.0,
                relative: 1e-12,
            },
            DecayEmissionLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn material(id: u64, coefficients: [f64; 2], scattering: bool) -> Arc<PhotonMassMaterial> {
    Arc::new(
        PhotonMassMaterial::new(
            format!("manufactured original material {id}"),
            id,
            DecayRadiation::Gamma,
            edges(),
            coefficients
                .into_iter()
                .enumerate()
                .map(|(group, coefficient)| {
                    let mut redistribution = vec![0.0; 2];
                    if scattering {
                        redistribution[group] = coefficient;
                    }
                    PhotonMassCollision {
                        total_m2_kg: coefficient,
                        heat_m2_kg: if scattering { 0.0 } else { coefficient },
                        redistribution_m2_kg: redistribution,
                    }
                })
                .collect(),
            PhotonSlabLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn layers(scattering_shield: bool) -> Vec<PhotonSlabLayer> {
    vec![
        PhotonSlabLayer {
            cell_id: CellId(1),
            material: material(1, [0.0, 0.0], false),
        },
        PhotonSlabLayer {
            cell_id: CellId(2),
            material: material(2, [0.5, 0.125], scattering_shield),
        },
        PhotonSlabLayer {
            cell_id: CellId(3),
            material: material(3, [0.25, 0.75], false),
        },
    ]
}

fn cells(shield_mass: f64) -> Vec<DepositedCell> {
    let mut cells = vec![
        cell(1, 1, 2.0, 1200.0, 0.0),
        cell(2, 2, shield_mass, 600.0, 0.0),
        cell(3, 3, 4.0, 1200.0, 0.0),
    ];
    let mut front = 0.0;
    for (cell, depth) in cells.iter_mut().zip([0.5, 1.0, 2.0]) {
        cell.geometry.start_m = [front, 0.0, 0.0];
        front += depth;
        cell.geometry.end_m = [front, 0.0, 0.0];
        cell.geometry.volume_m3 = depth;
        cell.geometry.surface_area_m2 = 2.0 + 4.0 * depth;
    }
    cells
}

fn thermal_grid(cells: Vec<DepositedCell>) -> DepositionSimulation<'static> {
    thermal(
        cells,
        &[
            (MaterialId(1), &HOT),
            (MaterialId(2), &COOL),
            (MaterialId(3), &HOT),
        ],
        vec![],
    )
}

fn plain(
    thermal: &DepositionSimulation<'_>,
    partition: &PreparedDecayRadiation,
) -> NuclearHeating<DecayCellSource> {
    NuclearHeating::new(
        thermal,
        vec![DecayCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            inventory: DecayInventory::new(Arc::clone(partition.model()), vec![INITIAL, 0.0, 0.0])
                .unwrap(),
        }],
        DecayBudget {
            propagation: DecayLimits::default(),
            nuclei_tolerance: 1e5,
            energy_tolerance_j: 1e-8,
        },
    )
    .unwrap()
}

fn bind(
    thermal: &DepositionSimulation<'_>,
    partition: Arc<PreparedDecayRadiation>,
    layers: Vec<PhotonSlabLayer>,
    geometry: PhotonSlabGeometry,
    checkpoint: impl FnMut() -> Result<(), FissionError>,
) -> Result<(NuclearHeating<DecayCellSource>, PhotonSlabResponse), FissionThermalError> {
    plain(thermal, &partition).with_photon_slab(
        thermal,
        "manufactured three-material slab".into(),
        geometry,
        layers,
        vec![DecayRadiationEmitter {
            cell_id: CellId(1),
            partition,
        }],
        vec![PhotonOrdinate {
            mu: 1.0,
            weight: 1.0,
        }],
        PhotonSlabLimits::default(),
        checkpoint,
    )
}

fn geometry() -> PhotonSlabGeometry {
    PhotonSlabGeometry {
        area_m2: 1.0,
        position_tolerance_m: 0.0,
    }
}

#[test]
fn g1_evaluated_parent_daughter_photons_heat_original_shield_and_receiver() {
    let original = thermal_grid(cells(3.0));
    let (source, transport) =
        bind(&original, partition(), layers(false), geometry(), || Ok(())).unwrap();
    let (heated, evolved, audit) = source
        .advance_to(&original, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let parent = INITIAL * 0.5 * 1e6 * EV_J;
    let daughter = INITIAL * (1.0 - 0.5 * (1.0 + std::f64::consts::LN_2)) * 2e6 * EV_J;
    // Transparent source: half the energy enters each normal direction. The
    // right-going half traverses original mass/area optical depths of each layer.
    let shield = 0.5 * (parent * (1.0 - (-1.5_f64).exp()) + daughter * (1.0 - (-0.375_f64).exp()));
    let receiver = 0.5
        * (parent * (-1.5_f64).exp() * (1.0 - (-1.0_f64).exp())
            + daughter * (-0.375_f64).exp() * (1.0 - (-3.0_f64).exp()));
    let escape = 0.5 * (parent + daughter)
        + 0.5 * (parent * (-2.5_f64).exp() + daughter * (-3.375_f64).exp());
    close(heated.temperature_k(CellId(1)).unwrap(), 300.0, 1e-13);
    close(
        heated.temperature_k(CellId(2)).unwrap(),
        300.0 + shield / 6.0,
        2e-12,
    );
    close(
        heated.temperature_k(CellId(3)).unwrap(),
        300.0 + receiver / 16.0,
        2e-12,
    );
    let account = evolved.audit();
    close(account.released.radiation_j, parent + daughter, 2e-12);
    close(account.radiation_deposited_j, shield + receiver, 2e-12);
    close(account.radiation_escaped_j, escape, 2e-12);
    close(account.thermal_credit_j, shield + receiver, 2e-12);
    close(audit.external_heat_j, shield + receiver, 2e-12);
    assert_eq!(account.released.deposited_j, 0.0);
    assert_eq!(account.radiation_ungrouped_j, 0.0);
    assert_eq!(
        transport.response.receiver_ids().collect::<Vec<_>>(),
        [1, 2, 3]
    );
    assert!(transport.max_first_flight_defect < 2e-14);
    assert!(transport.max_response_defect < 2e-14);
    for id in [CellId(1), CellId(2), CellId(3)] {
        assert_eq!(heated.cells()[&id].mass_kg, original.cells()[&id].mass_kg);
        assert_eq!(heated.cells()[&id].geometry, original.cells()[&id].geometry);
    }
    // The same prepared response accompanies the original subdivided history.
    let (middle, middle_source, _) = source
        .advance_to(&original, 0.4, controls(1.0), || Ok(()))
        .unwrap();
    let (split, _, _) = middle_source
        .advance_to(&middle, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    for id in [CellId(1), CellId(2), CellId(3)] {
        close(
            split.temperature_k(id).unwrap(),
            heated.temperature_k(id).unwrap(),
            2e-12,
        );
    }
}

#[test]
fn g3_original_mass_changes_shielding_and_zero_heat_layers_remain_bound() {
    let original = thermal_grid(cells(3.0));
    let denser = thermal_grid(cells(6.0));
    let partition = partition();
    let (_, thin) = bind(
        &original,
        Arc::clone(&partition),
        layers(false),
        geometry(),
        || Ok(()),
    )
    .unwrap();
    let (_, thick) = bind(
        &denser,
        Arc::clone(&partition),
        layers(false),
        geometry(),
        || Ok(()),
    )
    .unwrap();
    let unit = [1.0, 0.0];
    let emission = [(1, unit.as_slice())];
    let thin = thin.response.apply(&emission, 1e-12, || Ok(())).unwrap();
    let thick = thick.response.apply(&emission, 1e-12, || Ok(())).unwrap();
    close(
        thin.receiver_energy_j[&2],
        0.5 * (1.0 - (-1.5_f64).exp()),
        2e-14,
    );
    close(
        thick.receiver_energy_j[&2],
        0.5 * (1.0 - (-3.0_f64).exp()),
        2e-14,
    );
    assert!(thick.receiver_energy_j[&3] < thin.receiver_energy_j[&3]);
    let (source, scattering) =
        bind(&original, partition, layers(true), geometry(), || Ok(())).unwrap();
    let scattering = scattering
        .response
        .apply(&emission, 1e-12, || Ok(()))
        .unwrap();
    assert_eq!(scattering.receiver_energy_j[&2], 0.0);
    // This cell receives no heat but its original mass changes scattering.
    let failure = source
        .advance_to(&denser, 1.0, controls(1.0), || Ok(()))
        .err()
        .unwrap();
    assert_eq!(failure, FissionThermalError::BindingChanged(CellId(2)));
    assert_eq!(source.time_s(), 0.0);
    assert_eq!(
        source.source(CellId(1)).unwrap().inventory.populations()[0],
        INITIAL
    );
}

#[test]
fn g4_geometry_material_and_exact_callback_refusals_preserve_original_history() {
    let original = thermal_grid(cells(3.0));
    let partition = partition();
    let mut wrong_material = layers(false);
    wrong_material[1].material = material(99, [0.5, 0.125], false);
    assert!(matches!(
        bind(
            &original,
            Arc::clone(&partition),
            wrong_material,
            geometry(),
            || Ok(())
        ),
        Err(FissionThermalError::MaterialMismatch {
            cell_id: CellId(2),
            ..
        })
    ));
    assert!(
        bind(
            &original,
            Arc::clone(&partition),
            layers(false),
            PhotonSlabGeometry {
                area_m2: 2.0,
                position_tolerance_m: 0.0,
            },
            || Ok(())
        )
        .is_err()
    );
    let mut gapped_cells = cells(3.0);
    gapped_cells[1].geometry.start_m[0] += 0.125;
    gapped_cells[1].geometry.end_m[0] += 0.125;
    let gapped = thermal_grid(gapped_cells);
    assert!(
        bind(
            &gapped,
            Arc::clone(&partition),
            layers(false),
            geometry(),
            || Ok(())
        )
        .is_err()
    );
    let mut calls = 0;
    bind(
        &original,
        Arc::clone(&partition),
        layers(false),
        geometry(),
        || {
            calls += 1;
            Ok(())
        },
    )
    .unwrap();
    let stop = calls / 2;
    let mut calls = 0;
    let expected = FissionError::InvalidInput("original photon preparation interruption");
    let failure = bind(&original, partition, layers(false), geometry(), || {
        calls += 1;
        if calls == stop {
            Err(expected.clone())
        } else {
            Ok(())
        }
    })
    .err()
    .unwrap();
    assert_eq!(failure, FissionThermalError::Fission(expected));
    assert_eq!(original.time_s(), 0.0);
    for id in [CellId(1), CellId(2), CellId(3)] {
        assert_eq!(original.temperature_k(id).unwrap(), 300.0);
    }
}

#[test]
fn g3_refined_transport_credits_original_cells_once_and_replays_decay_history() {
    let original = thermal_grid(cells(3.0));
    let partition = partition();
    let source = plain(&original, &partition);
    let (coupled, refined) = source
        .clone()
        .with_refined_photon_slab(
            &original,
            "manufactured spatially refined original materials".into(),
            geometry(),
            layers(true),
            vec![DecayRadiationEmitter {
                cell_id: CellId(1),
                partition,
            }],
            vec![PhotonOrdinate {
                mu: 1.0,
                weight: 1.0,
            }],
            PhotonSlabLimits::default(),
            PhotonSlabRefinement {
                absolute_tolerance: 5e-3,
                relative_tolerance: 1e-3,
                max_levels: 4,
            },
            || Ok(()),
        )
        .unwrap();
    assert!(refined.levels >= 1);
    assert_eq!(refined.subdivisions_per_layer, 1 << refined.levels);
    assert!(refined.max_absolute_change <= 6e-3);
    assert_eq!(
        refined.transport.response.source_ids().collect::<Vec<_>>(),
        [1]
    );
    assert_eq!(
        refined
            .transport
            .response
            .receiver_ids()
            .collect::<Vec<_>>(),
        [1, 2, 3]
    );
    let parent = INITIAL * 0.5 * 1e6 * EV_J;
    let daughter = INITIAL * (1.0 - 0.5 * (1.0 + std::f64::consts::LN_2)) * 2e6 * EV_J;
    let photons = [parent, daughter];
    let expected = refined
        .transport
        .response
        .apply(&[(1, &photons)], 1e-10, || Ok(()))
        .unwrap();
    let (heated, evolved, audit) = coupled
        .advance_to(&original, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    assert_eq!(heated.cells().len(), 3);
    assert_eq!(evolved.sources().count(), 1);
    assert_eq!(expected.receiver_energy_j[&1], 0.0);
    assert_eq!(expected.receiver_energy_j[&2], 0.0);
    assert!(expected.receiver_energy_j[&3] > 0.0);
    close(heated.temperature_k(CellId(1)).unwrap(), 300.0, 1e-13);
    close(heated.temperature_k(CellId(2)).unwrap(), 300.0, 1e-13);
    close(
        heated.temperature_k(CellId(3)).unwrap(),
        300.0 + expected.receiver_energy_j[&3] / 16.0,
        2e-12,
    );
    close(
        evolved.audit().thermal_credit_j,
        expected.deposited_j,
        2e-12,
    );
    close(
        evolved.audit().radiation_escaped_j,
        expected.escaped_j,
        2e-12,
    );
    close(audit.external_heat_j, expected.deposited_j, 2e-12);
    // Photon-only subdivisions cannot alter the original nuclear history.
    let (_, plain_evolved, _) = source
        .advance_to(&original, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    assert_eq!(
        evolved.source(CellId(1)).unwrap().inventory.populations(),
        plain_evolved
            .source(CellId(1))
            .unwrap()
            .inventory
            .populations(),
    );
    let (middle, middle_source, _) = coupled
        .advance_to(&original, 0.4, controls(1.0), || Ok(()))
        .unwrap();
    let (split, split_source, _) = middle_source
        .advance_to(&middle, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    for id in [CellId(1), CellId(2), CellId(3)] {
        close(
            split.temperature_k(id).unwrap(),
            heated.temperature_k(id).unwrap(),
            2e-12,
        );
        assert_eq!(heated.cells()[&id].mass_kg, original.cells()[&id].mass_kg);
        assert_eq!(heated.cells()[&id].geometry, original.cells()[&id].geometry);
    }
    for (whole, split) in evolved
        .source(CellId(1))
        .unwrap()
        .inventory
        .populations()
        .iter()
        .zip(
            split_source
                .source(CellId(1))
                .unwrap()
                .inventory
                .populations(),
        )
    {
        close(*whole, *split, 32.0);
    }
}

#[test]
fn g4_transport_dimensions_refuse_before_macroscopic_coefficient_expansion() {
    let original = thermal_grid(cells(3.0));
    let partition = partition();
    let source = plain(&original, &partition);
    let prepare = |layers, limits| {
        source.clone().with_photon_slab(
            &original,
            "bounded original-material preparation".into(),
            geometry(),
            layers,
            vec![DecayRadiationEmitter {
                cell_id: CellId(1),
                partition: Arc::clone(&partition),
            }],
            vec![PhotonOrdinate {
                mu: 1.0,
                weight: 1.0,
            }],
            limits,
            || Ok(()),
        )
    };
    let mut excessive_coefficients = layers(false);
    // These admitted mass coefficients overflow when multiplied by this cell's
    // original density. The explicit zero storage budget must refuse first.
    excessive_coefficients[0].material = material(1, [f64::MAX, f64::MAX], false);
    assert!(matches!(
        prepare(
            excessive_coefficients,
            PhotonSlabLimits {
                max_values: 0,
                ..PhotonSlabLimits::default()
            }
        ),
        Err(FissionThermalError::Radiation(RadiationError::Budget(_)))
    ));
    // Raising the caller's cell cap cannot bypass the core's bounded LU size;
    // dimension refusal also precedes examining these repeated physical IDs.
    let excessive_layers = vec![layers(false)[0].clone(); 129];
    assert!(matches!(
        prepare(
            excessive_layers,
            PhotonSlabLimits {
                max_cells: 129,
                response: crate::radiation::grouped::GroupedRadiationResponseLimits {
                    max_receivers: 129,
                    ..PhotonSlabLimits::default().response
                },
                ..PhotonSlabLimits::default()
            }
        ),
        Err(FissionThermalError::Radiation(RadiationError::Budget(_)))
    ));
    assert_eq!(source.time_s(), 0.0);
    assert_eq!(
        source.source(CellId(1)).unwrap().inventory.populations()[0],
        INITIAL
    );
}

#[test]
fn g4_coordinate_tolerance_cannot_admit_a_reversed_diagonal_layer() {
    let n = 1.0 / fs_math::det::sqrt(3.0);
    let mut reversed: Vec<_> = cells(3.0).into_iter().take(2).collect();
    reversed[0].geometry.start_m = [0.0; 3];
    reversed[0].geometry.end_m = [n; 3];
    reversed[0].geometry.volume_m3 = 1.0;
    reversed[1].geometry.start_m = [n + 0.48; 3];
    reversed[1].geometry.end_m = [2.0 * n - 0.48; 3];
    reversed[1].geometry.volume_m3 = 1.0;
    let original = thermal_grid(reversed);
    // Every componentwise endpoint/seam discrepancy is below 0.49 m and the
    // tolerance is below half the declared 1 m depth. Nevertheless the second
    // original segment points backward; a positive-depth transport substitute
    // would silently reverse its physical geometry.
    let result = bind(
        &original,
        partition(),
        layers(false).into_iter().take(2).collect(),
        PhotonSlabGeometry {
            area_m2: 1.0,
            position_tolerance_m: 0.49,
        },
        || Ok(()),
    );
    assert!(matches!(result, Err(FissionThermalError::Invalid(_))));
    assert_eq!(original.time_s(), 0.0);
}
