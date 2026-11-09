//! Actual nozzle transport feeds one sourced droplet cloud and one closed budget.
#[path = "adaptive/tests.rs"]
mod adaptive;

use super::*;
use crate::graph::thermal::medium::mass_transfer::{
    evaporation::droplet::{
        chamber::{
            self as cloud, SourceDropletFilm, SourceDropletFilmSample, SourceGasDropletHeatSurface,
            SourceGasDropletSurface,
            transport::{
                DropletGasTransportError, GasDropletTransportSolve,
                advance_droplets_with_gas_transport, advance_droplets_with_gas_transport_and_heat,
            },
        },
        forced::SourceForcedSphericalDropletFilm,
    },
    sphere::ForcedSphereMassTransfer,
};
use fs_material::{
    evaporation::gas::chamber::{GasLiquidChamberSolve, heat::GasLiquidChamberHeatSolve},
    fluid_transport::FluidTransport,
    thermochemical_gas::energy::parcel::transfer::{
        nozzle::{
            GasNozzle,
            network::snapshot::{SnapshotNozzleError, SnapshotNozzleLink},
        },
        species::network::{
            SpeciesNetworkSolve,
            transport::{
                self as plumbing, GasTransportError, GasTransportPaths, GasTransportSolve,
            },
        },
    },
};

const REGION: u64 = 7;
const DONOR: u64 = 8;
const IDLE: u64 = 9;
const DURATION: f64 = 0.02;

fn parcel(
    catalog: &ThermochemicalGasModel,
    temperature: f64,
    pressure: f64,
    vapor: f64,
    volume: f64,
) -> ClosedGasParcel {
    let original = gas(catalog, temperature, pressure, vapor);
    ClosedGasParcel::new(original.model().clone(), original.point().clone(), volume).unwrap()
}
fn interval_controls() -> GasDropletTransportSolve {
    GasDropletTransportSolve {
        transport: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.2,
                mass_tolerance_kg: 1e-12,
                energy_tolerance_j: 1e-6,
            },
            entropy_tolerance_j_k: 1e-10,
        },
        phase: GasLiquidChamberSolve {
            surface: phase_controls().phase,
            maximum_surfaces: 4,
        },
        mass_tolerance_kg: 1e-12,
        energy_tolerance_j: 1e-6,
    }
}
fn endpoints(parcels: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    parcels
        .keys()
        .map(|&id| (id, phase_controls().phase.gas.gas))
        .collect()
}
fn nozzle(area_m2: f64) -> SnapshotNozzleLink {
    SnapshotNozzleLink {
        id: 31,
        parcels: [DONOR, REGION],
        nozzle: GasNozzle {
            area_m2,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: phase_controls().phase.gas.gas,
    }
}
fn joint_heat() -> GasLiquidChamberHeatSolve {
    GasLiquidChamberHeatSolve {
        maximum_iterations: 128,
    }
}

#[test]
fn original_nozzle_refreshes_cloud_sources_and_closes_all_inventories_once() {
    for with_heat in [false, true] {
        let catalog = catalog();
        let caloric = LiquidCaloric::new(&catalog);
        let binding = liquid_binding();
        let parcels = BTreeMap::from([
            (REGION, parcel(&catalog, 400.0, 100000.0, 0.003, 1e-5)),
            (DONOR, parcel(&catalog, 440.0, 130000.0, 0.04, 1e-4)),
            (IDLE, parcel(&catalog, 410.0, 90000.0, 0.01, 1e-5)),
        ]);
        let preserved = parcels.clone();
        let card = diffusion();
        let context = point(400.0, 100000.0);
        let prepared = SourceSphericalDropletFilm::new(
            &source(&card, &context).source,
            &binding.species,
            &parcels[&REGION],
            || false,
        )
        .unwrap();
        let forced = SourceForcedSphericalDropletFilm::try_new(
            &prepared,
            ForcedSphereMassTransfer::Frossling1938,
            1.0,
            0.05,
        )
        .unwrap();
        let temperature = |h| inverse(&caloric, h);
        let rows =
            [(11, 390.0, 1e-6), (29, 410.0, 8e-6)].map(|(id, t, mass)| SourceGasDropletSurface {
                id,
                law: &LiquidPhase,
                binding: &binding,
                droplet: sphere(&caloric, t, mass),
                temperature_from_enthalpy: &temperature,
                film: if id == 11 {
                    SourceDropletFilm::Stagnant(&prepared)
                } else {
                    SourceDropletFilm::Forced(forced)
                },
            });
        let initial_droplets = rows.map(|row| row.droplet);
        let warm_rows = rows.map(|surface| SourceGasDropletHeatSurface {
            surface,
            heat_transfer_coefficient_w_m2_k: 80.0,
            heat_controls: heat_controls(&caloric),
        });
        let controls = interval_controls();
        let bounds = endpoints(&parcels);
        let links = [nozzle(1e-8)];
        let paths = GasTransportPaths {
            nozzles: &links,
            ..Default::default()
        };
        let moved = plumbing::step_candidate_with_paths(
            &parcels,
            paths,
            DURATION,
            &bounds,
            controls.transport,
            || false,
        )
        .unwrap();
        assert!(moved.step().nozzles[0].transferred_mass_kg > 0.0);
        let direct = |gas: &ClosedGasParcel| {
            if with_heat {
                cloud::advance_droplets_with_gas_and_heat(
                    &warm_rows,
                    gas,
                    &catalog,
                    DURATION,
                    controls.phase,
                    joint_heat(),
                    || false,
                )
            } else {
                cloud::advance_droplets_with_gas(
                    &rows,
                    gas,
                    &catalog,
                    DURATION,
                    controls.phase,
                    || false,
                )
            }
            .unwrap()
        };
        let after_transport = direct(&moved.parcels()[&REGION]);
        let isolated = direct(&parcels[&REGION]);
        let next = if with_heat {
            advance_droplets_with_gas_transport_and_heat(
                &warm_rows,
                REGION,
                &parcels,
                &catalog,
                paths,
                &bounds,
                DURATION,
                controls,
                joint_heat(),
                || false,
            )
        } else {
            advance_droplets_with_gas_transport(
                &rows,
                REGION,
                &parcels,
                &catalog,
                paths,
                &bounds,
                DURATION,
                controls,
                || false,
            )
        }
        .unwrap();
        assert_eq!(&next.report().transport, moved.step());
        assert_eq!(next.droplets(), after_transport.original().droplets());
        assert_eq!(&next.parcels()[&REGION], after_transport.original().gas());
        assert_eq!(&next.report().coupled, after_transport.original().report());
        assert_eq!(next.films(), after_transport.films());
        assert_ne!(
            next.droplets()[0].phase.requested_transferred_mass_kg,
            isolated.original().droplets()[0]
                .phase
                .requested_transferred_mass_kg,
        );
        assert_eq!(next.parcels()[&DONOR], moved.parcels()[&DONOR]);
        assert_eq!(next.parcels()[&IDLE], parcels[&IDLE]);
        assert!(std::ptr::eq(next.catalog(), &catalog));
        assert_eq!(next.report().duration_s, DURATION);
        assert_eq!(next.report().transport.duration_s, DURATION);
        assert_eq!(next.report().coupled.exchange.duration_s, DURATION);
        assert_eq!(next.report().gas_region, REGION);
        assert_eq!(next.report().coupled.heat.is_some(), with_heat);

        let actual_gas = &moved.parcels()[&REGION];
        let bulk = actual_gas.state().state();
        assert_ne!(
            bulk.temperature,
            parcels[&REGION].state().state().temperature
        );
        let SourceDropletFilmSample::Stagnant(stagnant) = &next.films()[0].films[0] else {
            panic!("the explicitly stagnant first droplet changed its transport model");
        };
        assert_eq!(stagnant.temperature_k, bulk.temperature);
        assert_eq!(stagnant.total_pressure_pa, bulk.pressure);
        let d = 1e-5 + (bulk.temperature - 200.0) * (2e-5 / 700.0);
        near(stagnant.diffusion.diffusivity_m2_s, d, 1e-18, 1e-13);
        near(
            stagnant.mass_transfer_coefficient_m_s,
            d / stagnant.geometry.radius_m(),
            1e-15,
            1e-13,
        );
        let SourceDropletFilmSample::Forced(forced_sample) = &next.films()[1].films[0] else {
            panic!("the explicitly forced second droplet changed its transport model");
        };
        let expected_t = 0.5 * temperature(rows[1].droplet.specific_enthalpy_j_kg).unwrap()
            + 0.5 * bulk.temperature;
        let actual_material = actual_gas
            .model()
            .evaluate(&point(expected_t, bulk.pressure))
            .unwrap();
        assert_eq!(forced_sample.film_temperature_k, expected_t);
        assert_eq!(forced_sample.bulk_temperature_k, bulk.temperature);
        assert_eq!(forced_sample.total_pressure_pa, bulk.pressure);
        assert_eq!(
            forced_sample.transport,
            actual_material.transport_state().unwrap()
        );
        near(
            forced_sample.diffusion.diffusivity_m2_s,
            1e-5 + (expected_t - 200.0) * (2e-5 / 700.0),
            1e-18,
            1e-13,
        );
        assert_ne!(
            next.films()[1].films[0].mass_transfer_coefficient_m_s(),
            isolated.films()[1].films[0].mass_transfer_coefficient_m_s(),
        );

        // Independent original-to-final inventories include every finite gas
        // node and every liquid, with carried enthalpy counted only once.
        let mut liquid_mass = 0.0;
        let mut liquid_energy = 0.0;
        for (old, new) in rows.iter().zip(next.droplets()) {
            assert_eq!(new.id, old.id);
            liquid_mass += new.droplet.mass_kg - old.droplet.mass_kg;
            liquid_energy += new.droplet.mass_kg
                * (new.droplet.specific_enthalpy_j_kg + old.binding.liquid_enthalpy_offset_j_kg)
                - old.droplet.mass_kg
                    * (old.droplet.specific_enthalpy_j_kg
                        + old.binding.liquid_enthalpy_offset_j_kg);
        }
        let mut gas_mass = 0.0;
        let mut gas_energy = 0.0;
        let mut vapor_mass = 0.0;
        let mut carrier_mass = 0.0;
        let carrier = SpeciesId::new("A").unwrap();
        for (id, old) in &parcels {
            let new = &next.parcels()[id];
            gas_mass += new.mass_kg() - old.mass_kg();
            gas_energy += new.mass_kg() * new.state().specific_internal_energy_j_kg()
                - old.mass_kg() * old.state().specific_internal_energy_j_kg();
            vapor_mass += new.species_mass_kg(&binding.species).unwrap()
                - old.species_mass_kg(&binding.species).unwrap();
            carrier_mass +=
                new.species_mass_kg(&carrier).unwrap() - old.species_mass_kg(&carrier).unwrap();
            assert_eq!(new.volume_m3(), old.volume_m3());
            assert_eq!(new.point().axes().get("local-lot"), Some(&9.0));
        }
        near(liquid_mass + gas_mass, 0.0, controls.mass_tolerance_kg, 0.0);
        near(
            liquid_mass + vapor_mass,
            0.0,
            controls.mass_tolerance_kg,
            0.0,
        );
        near(carrier_mass, 0.0, controls.mass_tolerance_kg, 0.0);
        near(
            liquid_energy + gas_energy,
            0.0,
            controls.energy_tolerance_j,
            0.0,
        );
        assert!(next.report().absolute_mass_defect_kg <= controls.mass_tolerance_kg);
        assert!(next.report().absolute_energy_defect_j <= controls.energy_tolerance_j);
        assert_eq!(rows.map(|row| row.droplet), initial_droplets);
        assert_eq!(parcels, preserved);
    }
}

#[test]
fn empty_paths_recover_the_original_shared_phase_and_heat_and_cancel_without_publication() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let parcels = BTreeMap::from([
        (REGION, parcel(&catalog, 400.0, 100000.0, 0.003, 1e-5)),
        (IDLE, parcel(&catalog, 410.0, 90000.0, 0.01, 1e-5)),
    ]);
    let preserved = parcels.clone();
    let card = diffusion();
    let context = point(400.0, 100000.0);
    let prepared = SourceSphericalDropletFilm::new(
        &source(&card, &context).source,
        &binding.species,
        &parcels[&REGION],
        || false,
    )
    .unwrap();
    let temperature = |h| inverse(&caloric, h);
    let rows =
        [(11, 390.0, 1e-6), (29, 410.0, 8e-6)].map(|(id, t, mass)| SourceGasDropletSurface {
            id,
            law: &LiquidPhase,
            binding: &binding,
            droplet: sphere(&caloric, t, mass),
            temperature_from_enthalpy: &temperature,
            film: SourceDropletFilm::Stagnant(&prepared),
        });
    let initial_droplets = rows.map(|row| row.droplet);
    let warm_rows = rows.map(|surface| SourceGasDropletHeatSurface {
        surface,
        heat_transfer_coefficient_w_m2_k: 80.0,
        heat_controls: heat_controls(&caloric),
    });
    let controls = interval_controls();
    let bounds = endpoints(&parcels);
    for with_heat in [false, true] {
        let direct = if with_heat {
            cloud::advance_droplets_with_gas_and_heat(
                &warm_rows,
                &parcels[&REGION],
                &catalog,
                DURATION,
                controls.phase,
                joint_heat(),
                || false,
            )
        } else {
            cloud::advance_droplets_with_gas(
                &rows,
                &parcels[&REGION],
                &catalog,
                DURATION,
                controls.phase,
                || false,
            )
        }
        .unwrap();
        let next = if with_heat {
            advance_droplets_with_gas_transport_and_heat(
                &warm_rows,
                REGION,
                &parcels,
                &catalog,
                GasTransportPaths::default(),
                &bounds,
                DURATION,
                controls,
                joint_heat(),
                || false,
            )
        } else {
            advance_droplets_with_gas_transport(
                &rows,
                REGION,
                &parcels,
                &catalog,
                GasTransportPaths::default(),
                &bounds,
                DURATION,
                controls,
                || false,
            )
        }
        .unwrap();
        assert_eq!(next.droplets(), direct.original().droplets());
        assert_eq!(&next.parcels()[&REGION], direct.original().gas());
        assert_eq!(next.films(), direct.films());
        assert_eq!(&next.report().coupled, direct.original().report());
        assert_eq!(next.parcels()[&IDLE], parcels[&IDLE]);
        assert!(next.report().transport.nozzles.is_empty());
        assert_eq!(
            next.report().transport.exchange.total_energy_residual_j,
            0.0
        );
    }
    assert!(matches!(
        advance_droplets_with_gas_transport(
            &rows,
            REGION,
            &parcels,
            &catalog,
            GasTransportPaths::default(),
            &bounds,
            DURATION,
            controls,
            || true,
        ),
        Err(DropletGasTransportError::Cancelled)
    ));
    assert!(matches!(
        advance_droplets_with_gas_transport_and_heat(
            &warm_rows,
            REGION,
            &parcels,
            &catalog,
            GasTransportPaths::default(),
            &bounds,
            DURATION,
            controls,
            joint_heat(),
            || true,
        ),
        Err(DropletGasTransportError::Cancelled)
    ));
    assert_eq!(parcels, preserved);
    assert_eq!(rows.map(|row| row.droplet), initial_droplets);
}

#[test]
fn post_cloud_pressure_reversal_refuses_the_entire_network_then_a_closed_nozzle_retries() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let parcels = BTreeMap::from([
        (REGION, parcel(&catalog, 395.0, 100000.0, 0.003, 1e-5)),
        (DONOR, parcel(&catalog, 395.0, 100000.01, 0.003, 1e-4)),
        (IDLE, parcel(&catalog, 410.0, 90000.0, 0.01, 1e-5)),
    ]);
    let preserved = parcels.clone();
    let card = diffusion();
    let context = point(400.0, 100000.0);
    let prepared = SourceSphericalDropletFilm::new(
        &source(&card, &context).source,
        &binding.species,
        &parcels[&REGION],
        || false,
    )
    .unwrap();
    let temperature = |h| inverse(&caloric, h);
    let rows =
        [(11, 400.0, 1e-6), (29, 410.0, 2e-6)].map(|(id, t, mass)| SourceGasDropletSurface {
            id,
            law: &LiquidPhase,
            binding: &binding,
            droplet: sphere(&caloric, t, mass),
            temperature_from_enthalpy: &temperature,
            film: SourceDropletFilm::Stagnant(&prepared),
        });
    let initial_droplets = rows.map(|row| row.droplet);
    let warm_rows = rows.map(|surface| SourceGasDropletHeatSurface {
        surface,
        heat_transfer_coefficient_w_m2_k: 80.0,
        heat_controls: heat_controls(&caloric),
    });
    let controls = interval_controls();
    let bounds = endpoints(&parcels);
    let links = [nozzle(1e-12)];
    let paths = GasTransportPaths {
        nozzles: &links,
        ..Default::default()
    };
    let moved = plumbing::step_candidate_with_paths(
        &parcels,
        paths,
        DURATION,
        &bounds,
        controls.transport,
        || false,
    )
    .unwrap();
    assert!(moved.step().nozzles[0].transferred_mass_kg > 0.0);
    let direct = cloud::advance_droplets_with_gas_and_heat(
        &warm_rows,
        &moved.parcels()[&REGION],
        &catalog,
        DURATION,
        controls.phase,
        joint_heat(),
        || false,
    )
    .unwrap();
    assert!(
        direct.original().gas().state().state().pressure
            > moved.parcels()[&DONOR].state().state().pressure
    );
    let run = |paths| {
        advance_droplets_with_gas_transport_and_heat(
            &warm_rows,
            REGION,
            &parcels,
            &catalog,
            paths,
            &bounds,
            DURATION,
            controls,
            joint_heat(),
            || false,
        )
    };
    assert!(matches!(
        run(paths),
        Err(DropletGasTransportError::Transport(
            GasTransportError::Nozzle(SnapshotNozzleError::EquilibriumCrossed { id: 31 })
        ))
    ));
    assert_eq!(parcels, preserved);
    assert_eq!(rows.map(|row| row.droplet), initial_droplets);
    let closed = [nozzle(0.0)];
    let retry = run(GasTransportPaths {
        nozzles: &closed,
        ..Default::default()
    })
    .unwrap();
    let expected = cloud::advance_droplets_with_gas_and_heat(
        &warm_rows,
        &parcels[&REGION],
        &catalog,
        DURATION,
        controls.phase,
        joint_heat(),
        || false,
    )
    .unwrap();
    assert_eq!(retry.droplets(), expected.original().droplets());
    assert_eq!(&retry.parcels()[&REGION], expected.original().gas());
    assert_eq!(retry.films(), expected.films());
    assert_eq!(retry.parcels()[&DONOR], parcels[&DONOR]);
    assert_eq!(retry.parcels()[&IDLE], parcels[&IDLE]);
    assert_eq!(parcels, preserved);
    assert_eq!(rows.map(|row| row.droplet), initial_droplets);
}
