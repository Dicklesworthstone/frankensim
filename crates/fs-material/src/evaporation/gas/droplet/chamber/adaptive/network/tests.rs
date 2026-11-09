//! Full-network accuracy and accepted physical accounts on original material owners.
use super::*;
use super::super::{GasLiquidChamberSolve, SphericalDropletState, ThermochemicalGasModel};
use crate::evaporation::gas::{
    GasLiquidBinding,
    droplet::chamber::{
        GasDropletChamberReport, GasDropletSurface, advance_droplets_with_gas,
    },
    tests::{Phase, binding, catalog, controls, gas, liquid, temperature},
};
use crate::thermochemical_gas::{
    energy::{
        GasEnergySolve,
        parcel::transfer::{
            nozzle::{GasNozzle, network::snapshot::SnapshotNozzleLink},
            species::network::{
                SpeciesNetworkSolve,
                transport::{
                    self, GasTransportPaths, GasTransportSolve, GasTransportStep,
                },
            },
        },
    },
    tests::{point, species_with_energy_reference, transport as transport_material},
};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, mixture::FrozenIdealGasMixtureModelV1,
};

const REGION: u64 = 7;
const DONOR: u64 = 8;
const IDLE: u64 = 9;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}

fn solve() -> GasLiquidChamberSolve {
    GasLiquidChamberSolve {
        surface: controls(),
        maximum_surfaces: 4,
    }
}

fn adaptive() -> AdaptiveGasDropletControls {
    let mut result = AdaptiveGasDropletControls::new(solve());
    result.mass_conservation_tolerance_kg = 1e-12;
    result.energy_conservation_tolerance_j = 1e-6;
    result
}

fn parcel(
    source: &ThermochemicalGasModel,
    t: f64,
    p: f64,
    vapor_fraction: f64,
    volume: f64,
) -> ClosedGasParcel {
    let composition = gas(source, t, vapor_fraction * 100000.0, volume);
    ClosedGasParcel::new(
        composition.model().clone(),
        point(t, p).with("local-lot", 11.0).unwrap(),
        volume,
    )
    .unwrap()
}

fn pure_carrier(t: f64, p: f64, volume: f64) -> ClosedGasParcel {
    // B is absent from the source catalog, rather than present at zero fraction.
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport_material("A", 1.7e-5, 0.025)],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(source, point(t, p), volume).unwrap()
}

fn row<'a>(
    id: u64,
    binding: &'a GasLiquidBinding,
    mass: f64,
    t: f64,
    sherwood_number: f64,
) -> GasDropletSurface<'a> {
    GasDropletSurface {
        id,
        law: &Phase,
        binding,
        droplet: SphericalDropletState {
            mass_kg: mass,
            specific_enthalpy_j_kg: liquid(t).specific_enthalpy_j_kg,
        },
        temperature_from_enthalpy: &temperature,
        binary_diffusivity_m2_s: 2e-5,
        sherwood_number,
    }
}

fn transport_solve() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    }
}

#[allow(clippy::type_complexity)]
fn original_interval(
    before: &BTreeMap<u64, ClosedGasParcel>,
    rows: &[GasDropletSurface<'_>],
    source: &ThermochemicalGasModel,
    dt: f64,
) -> (
    BTreeMap<u64, ClosedGasParcel>,
    Vec<GasDropletSurfaceStep>,
    GasTransportStep,
    GasDropletChamberReport,
) {
    let bounds: BTreeMap<u64, GasEnergySolve> =
        before.keys().map(|&id| (id, controls().gas.gas)).collect();
    let nozzles = [SnapshotNozzleLink {
        id: 31,
        parcels: [DONOR, REGION],
        nozzle: GasNozzle {
            area_m2: 2e-8,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: controls().gas.gas,
    }];
    let (mut after, transport) = transport::step_candidate_with_paths(
        before,
        GasTransportPaths {
            nozzles: &nozzles,
            ..GasTransportPaths::default()
        },
        dt,
        &bounds,
        transport_solve(),
        || false,
    )
    .unwrap()
    .into_parts();
    let (droplets, coupled_gas, _, report) =
        advance_droplets_with_gas(rows, &after[&REGION], source, dt, solve(), || false)
            .unwrap()
            .into_parts();
    after.insert(REGION, coupled_gas);
    transport
        .admit_endpoint_driving_forces(&after, || false)
        .unwrap();
    (after, droplets, transport, report)
}

#[test]
fn empty_cloud_still_measures_remote_temperature_and_newly_arriving_species() {
    let source = catalog(0.0);
    let before = BTreeMap::from([
        (REGION, gas(&source, 320.0, 200.0, 1e-4)),
        (IDLE, pure_carrier(400.0, 100000.0, 1e-5)),
    ]);
    let coarse = before.clone();
    let mut fine = before.clone();
    // Hold mass and composition fixed while changing only a remote temperature.
    fine.insert(IDLE, pure_carrier(410.0, 100000.0 * 410.0 / 400.0, 1e-5));
    let mut accuracy = adaptive();
    accuracy.relative_tolerance = 0.0;
    accuracy.temperature_absolute_tolerance_k = 0.25;
    accuracy.mass_absolute_tolerance_kg = 1.0;
    let discrepancy = time_error(
        &[],
        &before,
        (&[], &coarse),
        (&[], &fine),
        accuracy,
        1.0,
        GasDropletNetworkLimits::default(),
        &mut || false,
    )
    .unwrap();
    near(discrepancy.gas_temperature_difference_k, 10.0, 0.0);
    near(discrepancy.normalized, 40.0, 0.0);
    near(discrepancy.maximum_temperature_difference_k, 0.0, 0.0);
    near(discrepancy.maximum_liquid_mass_difference_kg, 0.0, 0.0);
    assert_eq!(fine[&REGION], before[&REGION]);

    let vapor = SpeciesId::new("B").unwrap();
    assert!(
        !before[&IDLE]
            .species_inventory_kg()
            .unwrap()
            .contains_key(&vapor)
    );
    // Preserve the old carrier partial pressure while adding previously absent B.
    fine.insert(IDLE, parcel(&source, 400.0, 100000.0 / 0.9, 0.1, 1e-5));
    accuracy.mass_absolute_tolerance_kg = 1e-12;
    let arrival = fine[&IDLE].species_mass_kg(&vapor).unwrap();
    assert!(arrival > 0.0);
    let discrepancy = time_error(
        &[],
        &before,
        (&[], &coarse),
        (&[], &fine),
        accuracy,
        1.0,
        GasDropletNetworkLimits::default(),
        &mut || false,
    )
    .unwrap();
    near(discrepancy.gas_temperature_difference_k, 0.0, 0.0);
    near(discrepancy.maximum_species_mass_difference_kg, arrival, 1e-20);
    assert!(discrepancy.normalized >= arrival / accuracy.mass_absolute_tolerance_kg);
    assert_eq!(coarse, before);
}

#[test]
fn two_original_intervals_close_every_node_species_and_liquid_without_coarse_receipts() {
    for reference in [0.0, -1e7] {
        let source = catalog(reference);
        let binding = binding(reference);
        let before = BTreeMap::from([
            (REGION, parcel(&source, 320.0, 100000.0, 0.016, 1e-4)),
            (DONOR, parcel(&source, 350.0, 130000.0, 0.03, 1e-3)),
            (IDLE, parcel(&source, 300.0, 90000.0, 0.01, 2e-4)),
        ]);
        let saved = before.clone();
        let inputs = [
            row(2, &binding, 1e-6, 360.0, 5.0),
            row(9, &binding, 8e-6, 280.0, 2.0),
        ];
        let seeds: Vec<_> = inputs.iter().copied().map(GasDropletSeed::from).collect();
        let dt = 0.01;
        // This real coarse proposal is deliberately never recorded as accepted.
        let (_, _, _, coarse) = original_interval(&before, &inputs, &source, 2.0 * dt);
        let (middle, first, transport_first, report_first) =
            original_interval(&before, &inputs, &source, dt);
        let next_inputs: Vec<_> = inputs
            .iter()
            .zip(&first)
            .map(|(input, step)| GasDropletSurface {
                droplet: step.droplet,
                ..*input
            })
            .collect();
        let (after, final_rows, transport_last, report_last) =
            original_interval(&middle, &next_inputs, &source, dt);
        let limits = GasDropletNetworkLimits::default();
        let accuracy = adaptive();
        let audit = GasDropletNetworkWindow::new(
            &seeds, REGION, &before, accuracy, limits, &mut || false,
        )
        .unwrap()
        .record(&transport_first, &report_first, &first, &mut || false)
        .unwrap()
        .record(&transport_last, &report_last, &final_rows, &mut || false)
        .unwrap()
        .finish(&after, &final_rows, &mut || false)
        .unwrap();
        assert_eq!(audit.accepted_intervals, 2);
        near(audit.duration_s, 2.0 * dt, 0.0);
        let transferred =
            report_first.exchange.transferred_mass_kg + report_last.exchange.transferred_mass_kg;
        let evaporated =
            report_first.exchange.evaporated_mass_kg + report_last.exchange.evaporated_mass_kg;
        let condensed =
            report_first.exchange.condensed_mass_kg + report_last.exchange.condensed_mass_kg;
        near(audit.exchange.transferred_mass_kg, transferred, 1e-23);
        near(audit.exchange.evaporated_mass_kg, evaporated, 1e-23);
        near(audit.exchange.condensed_mass_kg, condensed, 1e-23);
        assert!(evaporated > 0.0 && condensed > 0.0);
        assert!(coarse.exchange.evaporated_mass_kg > 1e-23);
        assert!(
            (audit.exchange.evaporated_mass_kg - evaporated - coarse.exchange.evaporated_mass_kg)
                .abs()
                > 1e-23
        );
        assert_eq!(
            audit.parcels.iter().map(|p| p.id).collect::<Vec<_>>(),
            [REGION, DONOR, IDLE]
        );
        assert_eq!(audit.exchange.droplets.len(), inputs.len());
        let mut liquid_mass_change = 0.0;
        let mut energy_change = 0.0;
        for (((input, first), last), account) in inputs
            .iter()
            .zip(&first)
            .zip(&final_rows)
            .zip(&audit.exchange.droplets)
        {
            assert_eq!(last.geometry[0], first.geometry[1]);
            assert_eq!(account.id, input.id);
            assert_eq!(account.mass_kg, [input.droplet.mass_kg, last.droplet.mass_kg]);
            let change = last.droplet.mass_kg - input.droplet.mass_kg;
            liquid_mass_change += change;
            near(change + account.transferred_mass_kg, 0.0, 1e-18);
            near(
                account.transferred_mass_kg,
                first.exchange.transferred_mass_kg + last.exchange.transferred_mass_kg,
                1e-23,
            );
            energy_change += last.droplet.mass_kg * last.droplet.specific_enthalpy_j_kg
                - input.droplet.mass_kg * input.droplet.specific_enthalpy_j_kg
                + reference * change;
        }
        let mut gas_mass_change = 0.0;
        for node in &audit.parcels {
            let old = &before[&node.id];
            let new = &after[&node.id];
            assert_eq!(node.mass_kg, [old.mass_kg(), new.mass_kg()]);
            assert_eq!(old.volume_m3(), new.volume_m3());
            let change = new.mass_kg() - old.mass_kg();
            gas_mass_change += change;
            near(
                change,
                node.net_transported_mass_kg + node.phase_transferred_mass_kg,
                accuracy.mass_conservation_tolerance_kg,
            );
            near(node.mass_residual_kg, 0.0, accuracy.mass_conservation_tolerance_kg);
            near(
                node.maximum_species_mass_residual_kg,
                0.0,
                accuracy.mass_conservation_tolerance_kg,
            );
            if node.id == REGION {
                near(node.phase_transferred_mass_kg, transferred, 1e-23);
                assert!(node.net_transported_mass_kg > 0.0);
            } else {
                near(node.phase_transferred_mass_kg, 0.0, 0.0);
            }
            for species in &node.species {
                near(
                    species.mass_kg[0],
                    old.species_mass_kg(&species.species).unwrap(),
                    1e-18,
                );
                near(
                    species.mass_kg[1],
                    new.species_mass_kg(&species.species).unwrap(),
                    1e-18,
                );
                near(
                    species.mass_kg[1] - species.mass_kg[0],
                    species.transported_mass_kg + species.phase_transferred_mass_kg,
                    accuracy.mass_conservation_tolerance_kg,
                );
                near(species.mass_residual_kg, 0.0, accuracy.mass_conservation_tolerance_kg);
            }
            let physical_energy_change = new.mass_kg() * new.state().specific_internal_energy_j_kg()
                - old.mass_kg() * old.state().specific_internal_energy_j_kg();
            energy_change += physical_energy_change;
            near(
                node.gas_internal_energy_change_j,
                physical_energy_change,
                accuracy.energy_conservation_tolerance_j,
            );
        }
        let vapor = SpeciesId::new("B").unwrap();
        for species in &audit.exchange.species {
            let initial: f64 = before
                .values()
                .map(|gas| gas.species_mass_kg(&species.species).unwrap())
                .sum();
            let final_mass: f64 = after
                .values()
                .map(|gas| gas.species_mass_kg(&species.species).unwrap())
                .sum();
            near(species.gas_mass_kg[0], initial, 1e-18);
            near(species.gas_mass_kg[1], final_mass, 1e-18);
            let liquid_change = if species.species == vapor {
                liquid_mass_change
            } else {
                0.0
            };
            near(
                final_mass - initial + liquid_change,
                0.0,
                accuracy.mass_conservation_tolerance_kg,
            );
            near(species.closed_mass_residual_kg, 0.0, accuracy.mass_conservation_tolerance_kg);
        }
        near(liquid_mass_change + gas_mass_change, 0.0, accuracy.mass_conservation_tolerance_kg);
        near(energy_change, 0.0, accuracy.energy_conservation_tolerance_j);
        near(audit.exchange.mass_residual_kg, 0.0, accuracy.mass_conservation_tolerance_kg);
        near(audit.exchange.energy_residual_j, 0.0, accuracy.energy_conservation_tolerance_j);
        assert!(
            audit.exchange.sum_absolute_mass_residual_kg
                <= accuracy.mass_conservation_tolerance_kg
        );
        assert!(
            audit.exchange.sum_absolute_energy_residual_j
                <= accuracy.energy_conservation_tolerance_j
        );
        assert!(audit.parcels[1].net_transported_mass_kg < 0.0);
        near(audit.parcels[2].net_transported_mass_kg, 0.0, 0.0);
        near(audit.parcels[2].gas_internal_energy_change_j, 0.0, 0.0);
        assert_eq!(after[&IDLE], before[&IDLE]);
        assert_eq!(before, saved);
    }
}

#[test]
fn full_network_admission_bounds_new_species_and_preserves_original_identity() {
    let source = catalog(0.0);
    let before = BTreeMap::from([
        (REGION, gas(&source, 320.0, 200.0, 1e-4)),
        (IDLE, pure_carrier(400.0, 100000.0, 1e-5)),
    ]);
    let saved = before.clone();
    let limits = GasDropletNetworkLimits {
        maximum_parcels: 2,
        maximum_species_entries: 3,
    };
    admit_work(&before, limits, &mut || false).unwrap();
    for too_small in [
        GasDropletNetworkLimits { maximum_parcels: 1, ..limits },
        GasDropletNetworkLimits { maximum_species_entries: 2, ..limits },
    ] {
        assert!(admit_work(&before, too_small, &mut || false).is_err());
    }
    assert!(matches!(
        admit_work(&before, limits, &mut || true),
        Err(GasEvaporationError::Cancelled)
    ));
    let mut arrival = before.clone();
    arrival.insert(IDLE, parcel(&source, 400.0, 100000.0 / 0.9, 0.1, 1e-5));
    assert!(time_error(
        &[], &before, (&[], &before), (&[], &arrival), adaptive(), 1.0, limits, &mut || false,
    ).is_err());
    let mut renamed = before.clone();
    let original = renamed.remove(&IDLE).unwrap();
    renamed.insert(IDLE + 1, original);
    assert!(time_error(
        &[], &before, (&[], &before), (&[], &renamed), adaptive(), 1.0, limits, &mut || false,
    ).is_err());
    let mut resized = before.clone();
    resized.insert(IDLE, pure_carrier(400.0, 100000.0, 2e-5));
    assert!(time_error(
        &[], &before, (&[], &before), (&[], &resized), adaptive(), 1.0, limits, &mut || false,
    ).is_err());
    assert!(GasDropletNetworkWindow::new(
        &[], DONOR, &before, adaptive(), limits, &mut || false,
    ).is_err());
    let unused = GasDropletNetworkWindow::new(
        &[], REGION, &before, adaptive(), limits, &mut || false,
    ).unwrap();
    assert!(unused.finish(&before, &[], &mut || false).is_err());
    assert_eq!(before, saved);
}
