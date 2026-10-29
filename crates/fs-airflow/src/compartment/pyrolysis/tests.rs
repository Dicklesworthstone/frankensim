use super::*;
use crate::compartment::{Ambient, NetworkControl, NetworkError, Room, ZoneNetwork};
use fs_gas::multicomponent::compartment::ZoneSolve;
use fs_qty::{Acceleration, Area, Length, MassFlowRate, Pressure};
use fs_thermal::enthalpy::pyrolysis::{
    Arrhenius, CondensedCaloric, PyrolysisControls, PyrolysisMaterial, PyrolysisSimulation,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use std::cell::Cell;

#[path = "../../../examples/ventilated_compartment_fire/nasa9.rs"]
#[allow(dead_code)] // The reused source bank also exposes combustion-only constants.
mod nasa9;

fn room(gas: MixtureGas, air: Vec<f64>) -> TwoZone {
    TwoZone::new(
        gas,
        Area::new(1.0),
        Length::new(1.0),
        Length::new(0.0),
        Pressure::new(101325.0),
        [Temperature::new(300.0); 2],
        [air.clone(), air],
        || false,
    )
    .unwrap()
}
fn methane(gas: &MixtureGas) -> VolatileMixture {
    VolatileMixture::new(gas, &[1.0, 0.0, 0.0, 0.0, 0.0], || false).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

#[test]
fn original_nasa_mixture_enthalpy_and_reference_rates_match_existing_inflow() {
    let (gas, _, air) = nasa9::source().unwrap();
    let fractions = [0.5, 0.0, 0.0, 0.5, 0.0];
    let model = VolatileMixture::new(&gas, &fractions, || false).unwrap();
    let room = room(gas.clone(), air);
    for t in [300.0, 700.0, 1200.0] {
        let h = model
            .specific_enthalpy(Temperature::new(t), || false)
            .unwrap()
            .value();
        let original: f64 = gas
            .source()
            .models()
            .iter()
            .zip(fractions)
            .map(|(species, y)| {
                y * species
                    .evaluate_properties(Temperature::new(t))
                    .unwrap()
                    .h()
                    .value()
                    / species.molar_mass().value()
            })
            .sum();
        close(h, original, 2e-9);
        close(model.specific_enthalpy_j_kg(t).unwrap(), h, 0.0);
        let release = VolatileRelease {
            mass_kg: 0.002,
            enthalpy_j: 0.002 * h,
            consumed_virgin_mass_kg: 0.0025,
        };
        let converted = model
            .release_rates(
                &room,
                release,
                Time::new(0.25),
                VolatileTransferControl::default(),
                || false,
            )
            .unwrap();
        let inflow = room
            .inflow(MassFlowRate::new(0.008), &fractions, Temperature::new(t))
            .unwrap();
        assert_eq!(converted.rates.species_kg_s, inflow.species_kg_s);
        close(
            converted.rates.referenced_energy_w,
            inflow.referenced_energy_w,
            2e-10,
        );
        close(converted.enthalpy_residual_j, 0.0, 1e-12);
        close(converted.mass_residual_kg, 0.0, 1e-18);
    }
    assert_eq!(model.gas(), &gas);
    assert_eq!(model.mass_fractions(), fractions);
    assert_eq!(model.calorics().components().len(), 2);
}

#[test]
fn variable_temperature_release_retains_integrated_enthalpy_and_actual_species() {
    let (gas, _, air) = nasa9::source().unwrap();
    let mixture = methane(&gas);
    let mut room = room(gas, air);
    let before = room.clone();
    let release = VolatileRelease {
        mass_kg: 0.003,
        enthalpy_j: 0.001 * mixture.specific_enthalpy_j_kg(400.0).unwrap()
            + 0.002 * mixture.specific_enthalpy_j_kg(800.0).unwrap(),
        consumed_virgin_mass_kg: 0.003,
    };
    let injection = mixture
        .release_rates(
            &room,
            release,
            Time::new(0.4),
            VolatileTransferControl::default(),
            || false,
        )
        .unwrap();
    let mut rates = [LayerRates::default(); 2];
    injection.add_to(Layer::Upper, &mut rates).unwrap();
    room.advance(Time::new(0.4), rates, ZoneSolve::default(), || false)
        .unwrap();
    close(
        room.layer(Layer::Upper).species_masses()[0]
            - before.layer(Layer::Upper).species_masses()[0],
        release.mass_kg,
        1e-15,
    );
    for i in 1..5 {
        assert_eq!(
            room.layer(Layer::Upper).species_masses()[i],
            before.layer(Layer::Upper).species_masses()[i]
        );
    }
    close(
        room.physical_internal_energy().value() - before.physical_internal_energy().value(),
        release.enthalpy_j,
        2e-7,
    );
    assert!(
        (release.enthalpy_j - release.mass_kg * mixture.specific_enthalpy_j_kg(800.0).unwrap())
            .abs()
            > 500.0
    );
    assert_eq!(injection.release, release);
    assert_eq!(room.time(), Time::new(0.4));
}

#[test]
fn wrong_original_reference_composition_domain_and_empty_release_refuse() {
    let (gas, _, air) = nasa9::source().unwrap();
    let mixture = methane(&gas);
    let original = room(gas.clone(), air.clone());
    let changed = MixtureGas::new(
        gas.source().clone(),
        fs_gas::multicomponent::CaloricDomain {
            minimum_temperature: 250.0,
            ..gas.domain()
        },
        || false,
    )
    .unwrap();
    let changed = room(changed, air);
    assert_eq!(
        mixture.release_rates(
            &changed,
            VolatileRelease::default(),
            Time::new(1.0),
            VolatileTransferControl::default(),
            || false
        ),
        Err(VolatileGasError::SourceMismatch)
    );
    assert!(VolatileMixture::new(&gas, &[1.0], || false).is_err());
    assert!(VolatileMixture::new(&gas, &[0.9, 0.0, 0.0, 0.0, 0.0], || false).is_err());
    assert!(VolatileMixture::new(&gas, &[1.0, f64::NAN, 0.0, 0.0, 0.0], || false).is_err());
    assert!(
        mixture
            .specific_enthalpy(Temperature::new(199.0), || false)
            .is_err()
    );
    let empty = mixture
        .release_rates(
            &original,
            VolatileRelease::default(),
            Time::new(1.0),
            VolatileTransferControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(empty.rates, LayerRates::default());
    for release in [
        VolatileRelease {
            mass_kg: 0.0,
            enthalpy_j: 1.0,
            consumed_virgin_mass_kg: 0.0,
        },
        VolatileRelease {
            mass_kg: 1.0,
            enthalpy_j: 1e12,
            consumed_virgin_mass_kg: 1.0,
        },
        VolatileRelease {
            mass_kg: -1.0,
            enthalpy_j: 0.0,
            consumed_virgin_mass_kg: 1.0,
        },
        VolatileRelease {
            mass_kg: 1.0,
            enthalpy_j: 0.0,
            consumed_virgin_mass_kg: 0.5,
        },
    ] {
        assert!(
            mixture
                .release_rates(
                    &original,
                    release,
                    Time::new(1.0),
                    VolatileTransferControl::default(),
                    || false
                )
                .is_err()
        );
    }
    assert!(
        mixture
            .release_rates(
                &original,
                VolatileRelease::default(),
                Time::new(0.0),
                VolatileTransferControl::default(),
                || false
            )
            .is_err()
    );
}

#[test]
fn cancellation_during_species_conversion_leaves_gas_and_caller_rates_unchanged() {
    let (gas, _, air) = nasa9::source().unwrap();
    let mixture = methane(&gas);
    let room = room(gas, air);
    let before = room.clone();
    let release = VolatileRelease {
        mass_kg: 0.001,
        enthalpy_j: 0.001 * mixture.specific_enthalpy_j_kg(500.0).unwrap(),
        consumed_virgin_mass_kg: 0.001,
    };
    let calls = Cell::new(0);
    let mut rates = [LayerRates::default(); 2];
    let attempted = mixture
        .release_rates(
            &room,
            release,
            Time::new(0.1),
            VolatileTransferControl::default(),
            || {
                calls.set(calls.get() + 1);
                calls.get() == 4
            },
        )
        .and_then(|injection| injection.add_to(Layer::Upper, &mut rates));
    assert_eq!(attempted, Err(VolatileGasError::Cancelled));
    assert_eq!(rates, [LayerRates::default(); 2]);
    assert_eq!(room, before);
}

// Manufactured solid on the SAME datum as original methane. Its Cp and latent
// offset are explicit test inputs; no PMMA or other real-material claim is made.
struct Solid {
    h300: f64,
}
impl ThermalLaw for Solid {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.h300) / 500.0;
        if (200.0..=1500.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 500.0
    }
}
impl CondensedCaloric for Solid {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, DepositionError> {
        if (200.0..=1500.0).contains(&t) {
            Ok(self.h300 + 500.0 * (t - 300.0))
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}
fn solid<'a>(law: &'a Solid, volatile: &'a VolatileMixture) -> PyrolysisSimulation<'a> {
    let thermal = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(1),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: 0.02,
                specific_h_j_kg: law.specific_enthalpy_j_kg(700.0).unwrap(),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.01, 0.0, 0.0],
                    volume_m3: 1e-5,
                    surface_area_m2: 0.002,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 500.0,
                },
            }),
        }],
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.02,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    PyrolysisSimulation::new(
        thermal,
        &[PyrolysisMaterial {
            material_id: MaterialId(1),
            virgin: law,
            char_caloric: None,
            volatile,
            char_mass_yield: 0.0,
            kinetics: Arrhenius {
                prefactor_s_inv: 0.05,
                activation_energy_j_mol: 0.0,
            },
            temperature_range_k: [300.0, 1200.0],
            source: "manufactured solid with original NASA methane volatile",
        }],
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn original_reactive_solid_and_gas_close_mass_and_physical_formation_energy() {
    let (gas, _, air) = nasa9::source().unwrap();
    let volatile = methane(&gas);
    let law = Solid {
        h300: volatile.specific_enthalpy_j_kg(300.0).unwrap() - 1e6,
    };
    let solid = solid(&law, &volatile);
    let mut gas = room(gas, air);
    let mass0 = solid.mass_kg() + gas.layer(Layer::Upper).mass_kg();
    let energy0 = solid.energy_j() + gas.physical_internal_energy().value();
    let (candidate, audit) = solid
        .advance_to(0.2, PyrolysisControls::default(), || Ok(()))
        .unwrap();
    let release = audit.releases[&CellId(1)];
    let injection = volatile
        .release_rates(
            &gas,
            release,
            Time::new(0.2),
            VolatileTransferControl::default(),
            || false,
        )
        .unwrap();
    let mut rates = [LayerRates::default(); 2];
    injection.add_to(Layer::Upper, &mut rates).unwrap();
    gas.advance(Time::new(0.2), rates, ZoneSolve::default(), || false)
        .unwrap();
    close(
        candidate.mass_kg() + gas.layer(Layer::Upper).mass_kg(),
        mass0,
        2e-15,
    );
    close(
        candidate.energy_j() + gas.physical_internal_energy().value() - energy0,
        audit.thermal.external_heat_j,
        2e-7,
    );
    close(
        solid.mass_kg() - candidate.mass_kg(),
        release.mass_kg,
        2e-17,
    );
    assert!(release.mass_kg > 0.0);
    let endpoint_proxy = release.mass_kg
        * volatile
            .specific_enthalpy_j_kg(candidate.temperature_k(CellId(1)).unwrap())
            .unwrap();
    assert!((endpoint_proxy - release.enthalpy_j).abs() > 0.1);
    assert_eq!(solid.time_s(), 0.0);
    assert_eq!(candidate.time_s(), gas.time().value());
}

fn network(gas: MixtureGas, air: Vec<f64>) -> ZoneNetwork {
    ZoneNetwork::new(
        vec![Room {
            gas: room(gas, air.clone()),
            floor_elevation: Length::new(0.0),
        }],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(300.0),
            mass_fractions: air,
        },
        vec![],
        vec![],
        Acceleration::new(9.81),
        || false,
    )
    .unwrap()
}

#[test]
fn network_retry_recomputes_solid_candidate_and_includes_gas_source_in_pressure() {
    let (gas, _, air) = nasa9::source().unwrap();
    let volatile = methane(&gas);
    let law = Solid {
        h300: volatile.specific_enthalpy_j_kg(300.0).unwrap() - 1e6,
    };
    let solid = solid(&law, &volatile);
    let mut network = network(gas, air);
    let gas0 = network.rooms()[0].gas.clone();
    let attempts = Cell::new(0);
    let (report, (candidate, audit)) = network
        .step_with_sources(
            Time::new(0.2),
            NetworkControl::default(),
            |rooms, dt, rates, cancel| {
                attempts.set(attempts.get() + 1);
                assert_eq!(solid.time_s(), rooms[0].gas.time().value());
                let (candidate, audit) = solid
                    .advance_to(
                        solid.time_s() + dt.value(),
                        PyrolysisControls::default(),
                        || {
                            if cancel() {
                                Err(DepositionError::Cancelled)
                            } else {
                                Ok(())
                            }
                        },
                    )
                    .map_err(VolatileGasError::from)?;
                volatile
                    .release_rates(
                        &rooms[0].gas,
                        audit.releases[&CellId(1)],
                        dt,
                        VolatileTransferControl::default(),
                        &mut *cancel,
                    )?
                    .add_to(Layer::Upper, &mut rates[0])?;
                if attempts.get() == 1 {
                    return Err(NetworkError::PressureIterations);
                }
                Ok((candidate, audit))
            },
            || false,
        )
        .unwrap();
    assert_eq!(attempts.get(), 2);
    assert_eq!(report.elapsed, Time::new(0.1));
    assert!(report.maximum_pressure_residual_pa <= 1e-5);
    assert_eq!(report.maximum_split_pressure_change_pa, 0.0);
    let gas1 = &network.rooms()[0].gas;
    assert!(gas1.pressure().value() > gas0.pressure().value());
    close(
        gas1.layer(Layer::Upper).species_masses()[0],
        audit.volatile_mass_kg,
        1e-17,
    );
    close(
        candidate.energy_j() + gas1.physical_internal_energy().value()
            - solid.energy_j()
            - gas0.physical_internal_energy().value(),
        audit.thermal.external_heat_j,
        2e-7,
    );
    assert_eq!(candidate.time_s(), gas1.time().value());
    assert_eq!(solid.time_s(), 0.0);
}

#[test]
fn cancellation_after_solid_preparation_publishes_neither_owner() {
    let (gas, _, air) = nasa9::source().unwrap();
    let volatile = methane(&gas);
    let law = Solid {
        h300: volatile.specific_enthalpy_j_kg(300.0).unwrap() - 1e6,
    };
    let solid = solid(&law, &volatile);
    let mass0 = solid.mass_kg();
    let energy0 = solid.energy_j();
    let mut network = network(gas, air);
    let gas0 = network.rooms()[0].gas.clone();
    let cancel_requested = Cell::new(false);
    let result = network.step_with_sources(
        Time::new(0.1),
        NetworkControl::default(),
        |rooms, dt, rates, cancel| {
            let (candidate, audit) = solid
                .advance_to(dt.value(), PyrolysisControls::default(), || Ok(()))
                .map_err(VolatileGasError::from)?;
            volatile
                .release_rates(
                    &rooms[0].gas,
                    audit.releases[&CellId(1)],
                    dt,
                    VolatileTransferControl::default(),
                    cancel,
                )?
                .add_to(Layer::Upper, &mut rates[0])?;
            cancel_requested.set(true);
            Ok(candidate)
        },
        || cancel_requested.get(),
    );
    assert!(matches!(result, Err(NetworkError::Cancelled)));
    assert_eq!(network.rooms()[0].gas, gas0);
    assert_eq!(solid.time_s(), 0.0);
    assert_eq!(solid.mass_kg(), mass0);
    assert_eq!(solid.energy_j(), energy0);
}
