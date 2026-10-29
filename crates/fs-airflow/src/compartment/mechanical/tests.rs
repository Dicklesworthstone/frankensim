use super::*;
use crate::compartment::{Ambient, NetworkControl};
use crate::{FanCurve, FanPoint, LossElement, LossResistance, ToleranceBasis};
use fs_gas::multicomponent::MixtureGas;
use fs_gas::multicomponent::compartment::{LayerRates, TwoZone};
use fs_qty::{Acceleration, Temperature, Time};
use std::cell::Cell;

#[path = "../../../examples/ventilated_compartment_fire/nasa9.rs"]
#[allow(dead_code)]
mod nasa9;

fn source() -> SourceProvenance {
    SourceProvenance::new(
        "authored analytic fan fixture; not a product rating",
        "analytic-fan-v1",
    )
}
fn bank(count: usize, arrangement: FanArrangement, speed: f64) -> FanBank {
    FanBank::new(
        FanCurve::new(
            "linear",
            vec![
                FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(100.0)),
                FanPoint::new(VolumetricFlowRate::new(0.4), Pressure::new(60.0)),
                FanPoint::new(VolumetricFlowRate::new(1.0), Pressure::new(0.0)),
            ],
            source(),
            0.0,
            ToleranceBasis::Analytic,
            VolumetricFlowRate::new(0.0),
            (0.5, 2.0),
        )
        .unwrap(),
        count,
        arrangement,
        speed,
    )
    .unwrap()
}
fn losses() -> LossNetwork {
    // Original series/parallel equivalent: 50 + (200 || 200) = 100.
    let element = |name, r| {
        LossNetwork::Element(
            LossElement::new(
                name,
                LossResistance::new(r),
                0.0,
                source(),
                ToleranceBasis::Analytic,
            )
            .unwrap(),
        )
    };
    LossNetwork::series(vec![
        element("series", 50.0),
        LossNetwork::parallel(vec![element("left", 200.0), element("right", 200.0)]).unwrap(),
    ])
    .unwrap()
}
fn model() -> FanModel {
    make_model(bank(1, FanArrangement::Series, 1.0), 0.25)
}
fn make_model(bank: FanBank, rejected_loss_fraction: f64) -> FanModel {
    FanModel::new(
        bank,
        losses(),
        Density::new(1.0),
        Area::new(1.0),
        FanEnergy {
            total_efficiency: 0.5,
            rejected_loss_fraction,
            source: source(),
        },
        FanDomain {
            density_ratio: [0.2, 5.0],
            maximum_pressure_fraction: 0.02,
        },
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} versus {b:.16e}, tolerance {tolerance}"
    );
}

#[test]
fn analytic_total_pressure_root_density_speed_and_bank_affinity() {
    for (count, arrangement, speed) in [
        (1, FanArrangement::Series, 1.0),
        (1, FanArrangement::Series, 0.5),
        (2, FanArrangement::Series, 1.0),
        (2, FanArrangement::Parallel, 1.5),
    ] {
        let model = make_model(bank(count, arrangement, speed), 0.25);
        let (series, parallel) = if arrangement == FanArrangement::Series {
            (count as f64, 1.0)
        } else {
            (1.0, count as f64)
        };
        let a = 100.0 * series * speed * speed;
        let b = 100.0 * series * speed / parallel;
        let velocity_r = 0.5 * series / (parallel * parallel);
        let effective_r = 100.0 - velocity_r;
        let head = 5.0;
        // Independent positive root of R_eff Q^2 + b Q - (a-head) = 0.
        let expected_q = 2.0 * (a - head) / (b + (b * b + 4.0 * effective_r * (a - head)).sqrt());
        for density in [0.25, 1.0, 3.0] {
            let point = model
                .operating_point(
                    Pressure::new(100_000.0),
                    Pressure::new(100_000.0 + density * head),
                    Density::new(density),
                    FanControl::default(),
                    || false,
                )
                .unwrap();
            close(point.volume_flow.value(), expected_q, 1e-12);
            close(point.mass_rate.value(), density * expected_q, 3e-12);
            let expected_rise =
                density * (a - b * expected_q + velocity_r * expected_q * expected_q);
            close(point.total_pressure_rise.value(), expected_rise, 4e-10);
            close(
                point.duct_pressure_drop.value(),
                density * 100.0 * expected_q * expected_q,
                4e-10,
            );
            close(
                point.shaft_power.value(),
                expected_q * expected_rise / 0.5,
                2e-9,
            );
            close(
                point.gas_power.value() + point.rejected_power.value(),
                point.shaft_power.value(),
                1e-12,
            );
            close(
                point.rejected_power.value(),
                0.25 * 0.5 * point.shaft_power.value(),
                1e-12,
            );
            assert!(point.pressure_residual_pa.abs() <= 1e-7);
        }
    }
}

#[test]
fn pressure_load_reduces_flow_and_missing_forward_domain_refuses() {
    let m = model();
    let evaluate = |head| {
        m.operating_point(
            Pressure::new(100_000.0),
            Pressure::new(100_000.0 + head),
            Density::new(1.0),
            FanControl::default(),
            || false,
        )
    };
    let aiding = evaluate(-30.0).unwrap();
    let neutral = evaluate(0.0).unwrap();
    let adverse = evaluate(30.0).unwrap();
    assert!(aiding.mass_rate.value() > neutral.mass_rate.value());
    assert!(neutral.mass_rate.value() > adverse.mass_rate.value());
    assert!(matches!(
        evaluate(200.0),
        Err(FanError::NoOperatingPoint { .. })
    ));
    assert!(matches!(
        evaluate(-200.0),
        Err(FanError::NoOperatingPoint { .. })
    ));
    assert!(matches!(evaluate(100.0), Err(FanError::Invalid(_)))); // no invented shutoff power
    assert!(matches!(
        m.operating_point(
            Pressure::new(100_000.0),
            Pressure::new(100_000.0),
            Density::new(0.1),
            FanControl::default(),
            || false
        ),
        Err(FanError::Invalid(_))
    ));
    assert!(matches!(
        m.operating_point(
            Pressure::new(1000.0),
            Pressure::new(1000.0),
            Density::new(1.0),
            FanControl::default(),
            || false
        ),
        Err(FanError::Invalid(_))
    ));
    assert_eq!(
        m.operating_point(
            Pressure::new(100_000.0),
            Pressure::new(100_000.0),
            Density::new(1.0),
            FanControl {
                maximum_iterations: 1,
                ..FanControl::default()
            },
            || false
        ),
        Err(FanError::IterationLimit)
    );
    let mut polls = 0;
    assert_eq!(
        m.operating_point(
            Pressure::new(100_000.0),
            Pressure::new(100_000.0),
            Density::new(1.0),
            FanControl::default(),
            || {
                polls += 1;
                polls > 5
            }
        ),
        Err(FanError::Cancelled)
    );
}

#[test]
fn static_total_correction_cannot_hide_nonmonotone_or_unsourced_models() {
    let good = model();
    let mut bad_energy = good.energy().clone();
    bad_energy.source.identifier.clear();
    assert!(matches!(
        FanModel::new(
            good.bank().clone(),
            good.losses().clone(),
            Density::new(1.0),
            Area::new(1.0),
            bad_energy,
            good.domain(),
            || false
        ),
        Err(FanError::Invalid(_))
    ));
    // At tiny outlet area the increasing velocity pressure would overwhelm
    // this duct loss and decreasing static curve, admitting ambiguous roots.
    assert!(matches!(
        FanModel::new(
            good.bank().clone(),
            good.losses().clone(),
            Density::new(1.0),
            Area::new(0.01),
            good.energy().clone(),
            good.domain(),
            || false
        ),
        Err(FanError::Invalid(_))
    ));
    assert_eq!(
        FanModel::new(
            good.bank().clone(),
            good.losses().clone(),
            Density::new(1.0),
            Area::new(1.0),
            good.energy().clone(),
            good.domain(),
            || true
        ),
        Err(FanError::Cancelled)
    );
}

fn room(gas: MixtureGas, air: Vec<f64>, upper_temperature: f64, interface: f64) -> Room {
    Room {
        gas: TwoZone::new(
            gas,
            Area::new(10.0),
            Length::new(3.0),
            Length::new(interface),
            Pressure::new(101325.0),
            [Temperature::new(300.0), Temperature::new(upper_temperature)],
            [air.clone(), air],
            || false,
        )
        .unwrap(),
        floor_elevation: Length::new(0.0),
    }
}
fn network(two: bool, interface: f64) -> ZoneNetwork {
    let (gas, _, air) = nasa9::source().unwrap();
    let mut rooms = vec![room(gas.clone(), air.clone(), 600.0, interface)];
    if two {
        rooms.push(room(gas, air.clone(), 350.0, interface));
    }
    ZoneNetwork::new(
        rooms,
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
fn fan(from: Endpoint, to: Endpoint, height: f64) -> FanVent {
    FanVent {
        from: FanPort {
            endpoint: from,
            height: Length::new(height),
        },
        to: FanPort {
            endpoint: to,
            height: Length::new(height),
        },
        model: model(),
        controls: FanControl::default(),
        closed: false,
    }
}
fn energy(rooms: &[Room]) -> f64 {
    rooms
        .iter()
        .map(|r| r.gas.physical_internal_energy().value())
        .sum()
}
fn masses(rooms: &[Room]) -> Vec<f64> {
    let n = rooms[0].gas.gas().species_count();
    (0..n)
        .map(|s| {
            rooms
                .iter()
                .map(|r| {
                    r.gas.layer(Layer::Lower).species_masses()[s]
                        + r.gas.layer(Layer::Upper).species_masses()[s]
                })
                .sum()
        })
        .collect()
}

#[test]
fn network_solves_actual_fan_pressure_and_conserves_species_and_shaft_work() {
    let mut network = network(true, 1.0);
    let fan = fan(Endpoint::Room(0), Endpoint::Room(1), 2.5);
    network.set_fan_vents(vec![fan.clone()], || false).unwrap();
    let before = network.rooms().to_vec();
    let species = masses(&before);
    let report = network
        .step(
            Time::new(0.02),
            &[[LayerRates::default(); 2]; 2],
            NetworkControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.halvings, 0);
    let transfer = network.last_fan_transfers()[0].unwrap();
    assert_eq!(transfer.from_layer, Layer::Upper);
    assert_eq!(transfer.to_layer, Layer::Upper);
    let point = transfer.operating_point.unwrap();
    let dt = report.elapsed.value();
    close(
        energy(network.rooms()) - energy(&before),
        dt * point.gas_power.value(),
        2e-5,
    );
    close(
        report.external_physical_energy_j,
        dt * point.gas_power.value(),
        1e-10,
    );
    for (a, b) in masses(network.rooms()).iter().zip(species) {
        close(*a, b, 2e-13);
    }
    assert!(report.external_species_change_kg.iter().all(|x| *x == 0.0));
    assert!(report.mass_residual_kg.abs() < 1e-12);
    assert!(report.physical_energy_residual_j.abs() < 2e-5);
    let pressures: Vec<_> = network
        .rooms()
        .iter()
        .map(|r| r.gas.pressure().value())
        .collect();
    let inlet = super::super::column(
        &before,
        network.ambient(),
        network.ambient_density,
        &pressures,
        fan.from.endpoint,
    )
    .unwrap();
    let outlet = super::super::column(
        &before,
        network.ambient(),
        network.ambient_density,
        &pressures,
        fan.to.endpoint,
    )
    .unwrap();
    let reconstructed = fan
        .model
        .operating_point(
            inlet.pressure_at(fan.from.height, network.gravity).unwrap(),
            outlet.pressure_at(fan.to.height, network.gravity).unwrap(),
            inlet.density(crate::stratified::Layer::Upper),
            fan.controls,
            || false,
        )
        .unwrap();
    close(
        point.mass_rate.value(),
        reconstructed.mass_rate.value(),
        2e-8,
    );
    // The endpoint solve differs materially from freezing the fan at old pressure.
    let initial = fan
        .model
        .operating_point(
            super::super::column(
                &before,
                network.ambient(),
                network.ambient_density,
                &[101325.0; 2],
                fan.from.endpoint,
            )
            .unwrap()
            .pressure_at(fan.from.height, network.gravity)
            .unwrap(),
            super::super::column(
                &before,
                network.ambient(),
                network.ambient_density,
                &[101325.0; 2],
                fan.to.endpoint,
            )
            .unwrap()
            .pressure_at(fan.to.height, network.gravity)
            .unwrap(),
            Density::new(before[0].gas.layer(Layer::Upper).density_kg_m3().unwrap()),
            fan.controls,
            || false,
        )
        .unwrap();
    assert!((initial.mass_rate.value() - point.mass_rate.value()).abs() > 1e-3);
}

#[test]
fn ambient_supply_exhaust_and_absent_layer_ports_use_original_enthalpy() {
    for supply in [true, false] {
        let mut network = network(false, 0.0); // only upper layer is present
        let (from, to) = if supply {
            (Endpoint::Ambient, Endpoint::Room(0))
        } else {
            (Endpoint::Room(0), Endpoint::Ambient)
        };
        network
            .set_fan_vents(vec![fan(from, to, 0.0)], || false)
            .unwrap();
        let before = network.rooms().to_vec();
        let before_mass = masses(&before);
        let report = network
            .step(
                Time::new(0.005),
                &[[LayerRates::default(); 2]],
                NetworkControl::default(),
                || false,
            )
            .unwrap();
        let transfer = network.last_fan_transfers()[0].unwrap();
        assert_eq!(
            if supply {
                transfer.to_layer
            } else {
                transfer.from_layer
            },
            Layer::Upper
        );
        let point = transfer.operating_point.unwrap();
        let carried = if supply {
            before[0]
                .gas
                .inflow(
                    point.mass_rate,
                    &network.ambient().mass_fractions,
                    network.ambient().temperature,
                )
                .unwrap()
        } else {
            before[0]
                .gas
                .carried_rate(Layer::Upper, point.mass_rate)
                .unwrap()
        };
        let physical_enthalpy = carried.referenced_energy_w
            + carried
                .species_kg_s
                .iter()
                .zip(before[0].gas.gas().energy_references())
                .map(|(rate, reference)| rate * reference)
                .sum::<f64>();
        let sign = if supply { 1.0 } else { -1.0 };
        let expected_energy = report.elapsed.value()
            * (sign * physical_enthalpy + if supply { point.gas_power.value() } else { 0.0 });
        close(report.external_physical_energy_j, expected_energy, 1e-9);
        close(
            energy(network.rooms()) - energy(&before),
            expected_energy,
            1e-5,
        );
        for (s, (after, old)) in masses(network.rooms()).iter().zip(before_mass).enumerate() {
            close(
                after - old,
                sign * report.elapsed.value() * carried.species_kg_s[s],
                2e-13,
            );
        }
        assert!(point.shaft_power.value() > 0.0);
        close(
            network.rooms()[0].gas.layer(Layer::Lower).mass_kg(),
            0.0,
            0.0,
        );
    }
}

#[test]
fn stopped_isolation_geometry_updates_and_cancellation_are_transactional() {
    let mut network = network(true, 1.0);
    let mut closed = fan(Endpoint::Room(0), Endpoint::Room(1), 2.0);
    closed.closed = true;
    network
        .set_fan_vents(vec![closed.clone()], || false)
        .unwrap();
    let e = energy(network.rooms());
    let m = masses(network.rooms());
    let report = network
        .step(
            Time::new(0.02),
            &[[LayerRates::default(); 2]; 2],
            NetworkControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(
        network.last_fan_transfers()[0].unwrap().operating_point,
        None
    );
    close(report.external_physical_energy_j, 0.0, 0.0);
    close(energy(network.rooms()), e, 1e-5);
    assert_eq!(masses(network.rooms()), m);
    let vents = network.fan_vents().to_vec();
    let diagnostics = network.last_fan_transfers().to_vec();
    let mut bad = closed.clone();
    bad.to.height = Length::new(2.1);
    assert!(network.set_fan_vents(vec![bad], || false).is_err());
    assert_eq!(network.fan_vents(), vents);
    assert_eq!(network.last_fan_transfers(), diagnostics);
    assert_eq!(
        network.set_fan_vents(vec![], || true),
        Err(NetworkError::Cancelled)
    );
    assert_eq!(network.fan_vents(), vents);
    closed.closed = false;
    network.set_fan_vents(vec![closed], || false).unwrap();
    let before = network.rooms().to_vec();
    let diagnostics = network.last_fan_transfers().to_vec();
    let cancel = Cell::new(false);
    let result = network.step_with_stages(
        Time::new(0.01),
        NetworkControl::default(),
        |_, _, _, _| Ok(()),
        |_, _, _, _, _| {
            cancel.set(true);
            Ok(())
        },
        || cancel.get(),
    );
    assert_eq!(result, Err(NetworkError::Cancelled));
    assert_eq!(network.rooms(), before);
    assert_eq!(network.last_fan_transfers(), diagnostics);
}

#[test]
fn failed_full_interval_recomputes_and_publishes_only_accepted_fan_transfer() {
    let mut network = network(true, 1.0);
    network
        .set_fan_vents(vec![fan(Endpoint::Room(0), Endpoint::Room(1), 2.5)], || {
            false
        })
        .unwrap();
    let before = network.rooms().to_vec();
    let attempts = Cell::new(0);
    let (report, ()) = network
        .step_with_stages(
            Time::new(0.02),
            NetworkControl::default(),
            |_, _, _, _| {
                attempts.set(attempts.get() + 1);
                Ok(())
            },
            |_, dt, _, _, _| {
                if dt.value() > 0.01 {
                    Err(NetworkError::Invalid("force first whole-step rejection"))
                } else {
                    Ok(())
                }
            },
            || false,
        )
        .unwrap();
    assert_eq!(attempts.get(), 2);
    assert_eq!(report.halvings, 1);
    assert_eq!(report.elapsed, Time::new(0.01));
    let point = network.last_fan_transfers()[0]
        .unwrap()
        .operating_point
        .unwrap();
    close(
        energy(network.rooms()) - energy(&before),
        0.01 * point.gas_power.value(),
        2e-5,
    );
    close(
        before[0].gas.layer(Layer::Upper).mass_kg()
            - network.rooms()[0].gas.layer(Layer::Upper).mass_kg(),
        0.01 * point.mass_rate.value(),
        5e-14,
    );
}
