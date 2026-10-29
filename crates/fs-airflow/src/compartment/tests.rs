use super::*;
use crate::stratified::{DischargeCoefficient, Layer as VentLayer};
use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::{CaloricDomain, MixtureGas};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Area, MolarMass};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};

fn model(name: &str, mass: f64, cp: f64, linear: f64, formation: f64) -> Nasa9StandardStateModelV1 {
    let coefficients = [0.0, 0.0, cp, linear, 0.0, 0.0, 0.0, formation, 0.0];
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100_000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 2_000.0),
        sources: vec![],
        provenance: Provenance {
            source: "analytic compartment-network caloric fixture".into(),
            license: "fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100_000.0),
            ElementalReferenceIdV1::new("compartment-network-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn gas(equal: bool) -> MixtureGas {
    let models = vec![
        model("A", 0.028, 3.5, if equal { 0.0 } else { 0.0002 }, 1_000.0),
        model(
            "B",
            if equal { 0.028 } else { 0.032 },
            if equal { 3.5 } else { 4.0 },
            0.0,
            -2_000.0,
        ),
    ];
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1, 1]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 0]).unwrap(),
        models,
    )
    .unwrap();
    MixtureGas::new(
        bank,
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(2_000.0)),
        || false,
    )
    .unwrap()
}
fn room(gas: MixtureGas, pressure: f64, temperatures: [f64; 2], fractions: [[f64; 2]; 2]) -> Room {
    Room {
        gas: TwoZone::new(
            gas,
            Area::new(10.0),
            Length::new(3.0),
            Length::new(1.2),
            Pressure::new(pressure),
            temperatures.map(Temperature::new),
            fractions.map(Vec::from),
            || false,
        )
        .unwrap(),
        floor_elevation: Length::new(0.0),
    }
}
fn ambient() -> Ambient {
    Ambient {
        pressure: Pressure::new(100_000.0),
        reference_height: Length::new(0.0),
        temperature: Temperature::new(300.0),
        mass_fractions: vec![0.75, 0.25],
    }
}
fn opening(width: f64) -> WallOpening {
    WallOpening::new(
        Length::new(0.0),
        Length::new(2.5),
        Length::new(width),
        DischargeCoefficient::sharp_edged().unwrap(),
    )
    .unwrap()
}
fn controls() -> NetworkControl {
    NetworkControl {
        pressure_tolerance_pa: 1e-6,
        zone: ZoneSolve {
            absolute_energy_tolerance_j: 1e-8,
            relative_energy_tolerance: 1e-12,
            ..ZoneSolve::default()
        },
        ..NetworkControl::default()
    }
}
fn inventories(rooms: &[Room]) -> ([f64; 2], f64) {
    let mut species = [0.0; 2];
    let mut energy = 0.0;
    for room in rooms {
        for layer in [Layer::Lower, Layer::Upper] {
            for (i, mass) in room.gas.layer(layer).species_masses().iter().enumerate() {
                species[i] += mass;
            }
        }
        energy += room.gas.physical_internal_energy().value();
    }
    (species, energy)
}
fn assert_admissible(rooms: &[Room]) {
    for room in rooms {
        assert!(room.gas.pressure().value() > 0.0);
        assert!((0.0..=3.0).contains(&room.gas.interface_height().value()));
        let mut volume = 0.0;
        for layer in [Layer::Lower, Layer::Upper] {
            let state = room.gas.layer(layer);
            assert!(
                state
                    .species_masses()
                    .iter()
                    .all(|m| m.is_finite() && *m >= 0.0)
            );
            volume += state.volume().value();
            if state.mass_kg() > 0.0 {
                let sum: f64 = state
                    .species_masses()
                    .iter()
                    .map(|m| m / state.mass_kg())
                    .sum();
                assert!((sum - 1.0).abs() < 2e-15);
                assert!((200.0..=2_000.0).contains(&state.temperature().value()));
            }
        }
        assert!((volume - 30.0).abs() < 1e-12);
    }
}
fn assert_balance(report: &NetworkStepReport) {
    assert!(
        report.species_residual_kg[..2]
            .iter()
            .all(|r| r.abs() < 2e-11)
    );
    assert!(report.mass_residual_kg.abs() < 4e-11);
    assert!(
        report.physical_energy_residual_j.abs() < 1e-4,
        "energy residual {report:?}"
    );
    assert!(report.maximum_pressure_residual_pa <= controls().pressure_tolerance_pa);
}
fn communicating_network() -> ZoneNetwork {
    let source = gas(false);
    ZoneNetwork::new(
        vec![
            room(
                source.clone(),
                100_020.0,
                [300.0, 500.0],
                [[0.8, 0.2], [0.1, 0.9]],
            ),
            room(source, 100_000.0, [310.0, 380.0], [[0.3, 0.7], [0.6, 0.4]]),
        ],
        ambient(),
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Room(1),
            opening: opening(0.2),
        }],
        vec![],
        Acceleration::new(0.0),
        || false,
    )
    .unwrap()
}

#[test]
fn communicating_rooms_close_each_species_and_energy_while_pressure_relaxes() {
    let mut network = communicating_network();
    let before = network.rooms().to_vec();
    let (species, energy) = inventories(&before);
    let report = network
        .step(
            Time::new(0.02),
            &[[LayerRates::default(); 2]; 2],
            controls(),
            || false,
        )
        .unwrap();
    assert_balance(&report);
    assert_admissible(network.rooms());
    let (after_species, after_energy) = inventories(network.rooms());
    for i in 0..2 {
        assert!((after_species[i] - species[i]).abs() < 2e-11);
    }
    assert!((after_energy - energy).abs() < 1e-4);
    let old_gap = before[0].gas.pressure().value() - before[1].gas.pressure().value();
    let new_gap =
        network.rooms()[0].gas.pressure().value() - network.rooms()[1].gas.pressure().value();
    assert!(new_gap > 0.0 && new_gap < old_gap);
    assert!(
        network.rooms()[0].gas.layer(Layer::Upper).mass_kg()
            < before[0].gas.layer(Layer::Upper).mass_kg()
    );
    assert!(
        network.rooms()[1].gas.layer(Layer::Lower).mass_kg()
            > before[1].gas.layer(Layer::Lower).mass_kg()
    );
    assert_eq!(network.rooms()[0].gas.time(), report.elapsed);
    assert_eq!(network.rooms()[1].gas.time(), report.elapsed);
    assert!(report.iterations > 1);
}

#[test]
fn buoyant_ambient_exchange_retains_both_directions_and_their_species_enthalpies() {
    let initial = room(
        gas(false),
        100_000.0,
        [300.0, 650.0],
        [[0.9, 0.1], [0.2, 0.8]],
    );
    let mut network = ZoneNetwork::new(
        vec![initial.clone()],
        ambient(),
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Ambient,
            opening: opening(0.8),
        }],
        vec![],
        Acceleration::new(9.81),
        || false,
    )
    .unwrap();
    let report = network
        .step(
            Time::new(0.05),
            &[[LayerRates::default(); 2]],
            controls(),
            || false,
        )
        .unwrap();
    assert_balance(&report);
    assert_admissible(network.rooms());
    // Rebuild the accepted frozen-T/composition hydrostatic profiles independently
    // from the original room and the accepted endpoint pressure.
    let endpoint_pressure = network.rooms()[0].gas.pressure().value();
    let scale = endpoint_pressure / initial.gas.pressure().value();
    let inside = TwoLayerColumn::new(
        Length::new(0.0),
        Pressure::new(endpoint_pressure),
        initial.gas.interface_height(),
        Density::new(initial.gas.layer(Layer::Lower).density_kg_m3().unwrap() * scale),
        Density::new(initial.gas.layer(Layer::Upper).density_kg_m3().unwrap() * scale),
    )
    .unwrap();
    let molar_inverse: f64 = network
        .ambient()
        .mass_fractions
        .iter()
        .zip(initial.gas.gas().source().models())
        .map(|(y, source)| y / source.molar_mass().value())
        .sum();
    let rho =
        100_000.0 / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 300.0 * molar_inverse);
    let outside = TwoLayerColumn::new(
        Length::new(0.0),
        Pressure::new(100_000.0),
        Length::new(0.0),
        Density::new(rho),
        Density::new(rho),
    )
    .unwrap();
    let slabs = stratified::wall_vent_flow(
        &opening(0.8),
        inside,
        outside,
        Acceleration::new(9.81),
        || false,
    )
    .unwrap();
    assert!(slabs.slabs().iter().any(|s| s.mass_rate.value() < 0.0));
    assert!(slabs.slabs().iter().any(|s| s.mass_rate.value() > 0.0));
    let mut expected = [LayerRates::default(); 2];
    for slab in slabs.slabs() {
        let layer = match slab.left_layer {
            VentLayer::Lower => Layer::Lower,
            VentLayer::Upper => Layer::Upper,
        };
        let change = if slab.mass_rate.value() < 0.0 {
            initial
                .gas
                .inflow(
                    MassFlowRate::new(-slab.mass_rate.value()),
                    &network.ambient().mass_fractions,
                    network.ambient().temperature,
                )
                .unwrap()
        } else {
            initial
                .gas
                .carried_rate(layer, slab.mass_rate)
                .unwrap()
                .opposite()
        };
        expected[layer.index()].add(&change).unwrap();
    }
    let dt = report.elapsed.value();
    for layer in [Layer::Lower, Layer::Upper] {
        for i in 0..2 {
            let actual = network.rooms()[0].gas.layer(layer).species_masses()[i]
                - initial.gas.layer(layer).species_masses()[i];
            // The final pressure residual allows at most a tiny difference from
            // the last implicit pressure iterate used for the accepted flux.
            assert!((actual - dt * expected[layer.index()].species_kg_s[i]).abs() < 2e-8);
        }
    }
    let expected_heat = expected[0].referenced_energy_w + expected[1].referenced_energy_w;
    let actual_heat = network.rooms()[0].gas.sensible_internal_energy().value()
        - initial.gas.sensible_internal_energy().value();
    assert!((actual_heat - dt * expected_heat).abs() < 0.02);
    assert!(
        network.rooms()[0].gas.layer(Layer::Upper).species_masses()[1]
            < initial.gas.layer(Layer::Upper).species_masses()[1]
    );
}

#[test]
fn mechanical_transfer_preserves_donor_species_and_accounts_for_external_fan_work() {
    let source = gas(false);
    let rooms = vec![
        room(
            source.clone(),
            100_000.0,
            [300.0, 500.0],
            [[0.8, 0.2], [0.1, 0.9]],
        ),
        room(source, 100_000.0, [310.0, 380.0], [[0.3, 0.7], [0.6, 0.4]]),
    ];
    let before = rooms.clone();
    let (_, energy) = inventories(&before);
    let mut network = ZoneNetwork::new(
        rooms,
        ambient(),
        vec![],
        vec![MechanicalVent {
            from: Endpoint::Room(0),
            from_layer: Layer::Upper,
            to: Endpoint::Room(1),
            to_layer: Layer::Lower,
            mass_rate: MassFlowRate::new(0.02),
            deposited_power: Power::new(100.0),
        }],
        Acceleration::new(9.81),
        || false,
    )
    .unwrap();
    let report = network
        .step(
            Time::new(0.2),
            &[[LayerRates::default(); 2]; 2],
            controls(),
            || false,
        )
        .unwrap();
    assert_balance(&report);
    let dt = report.elapsed.value();
    for i in 0..2 {
        let moved = dt * 0.02 * before[0].gas.layer(Layer::Upper).species_masses()[i]
            / before[0].gas.layer(Layer::Upper).mass_kg();
        let lost = before[0].gas.layer(Layer::Upper).species_masses()[i]
            - network.rooms()[0].gas.layer(Layer::Upper).species_masses()[i];
        let gained = network.rooms()[1].gas.layer(Layer::Lower).species_masses()[i]
            - before[1].gas.layer(Layer::Lower).species_masses()[i];
        assert!((lost - moved).abs() < 3e-14);
        assert!((gained - moved).abs() < 3e-14);
    }
    assert!((inventories(network.rooms()).1 - energy - 100.0 * dt).abs() < 1e-4);
    assert_admissible(network.rooms());
}

#[test]
fn sealed_heated_room_recovers_constant_cp_pressure_and_keeps_all_species() {
    let initial = room(
        gas(true),
        100_000.0,
        [300.0, 450.0],
        [[0.8, 0.2], [0.3, 0.7]],
    );
    let energy = initial.gas.physical_internal_energy().value();
    let mut network = ZoneNetwork::new(
        vec![initial.clone()],
        ambient(),
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Ambient,
            opening: opening(0.0),
        }],
        vec![],
        Acceleration::new(9.81),
        || false,
    )
    .unwrap();
    let mut rates = [LayerRates::default(); 2];
    rates[1].add_heat(Power::new(5_000.0)).unwrap();
    let report = network
        .step(Time::new(0.2), &[rates], controls(), || false)
        .unwrap();
    assert_balance(&report);
    let after = &network.rooms()[0].gas;
    let heat = 5_000.0 * report.elapsed.value();
    assert!((after.pressure().value() - (100_000.0 + 0.4 * heat / 30.0)).abs() < 5e-6);
    assert!((after.physical_internal_energy().value() - energy - heat).abs() < 1e-4);
    for layer in [Layer::Lower, Layer::Upper] {
        assert_eq!(
            after.layer(layer).species_masses(),
            initial.gas.layer(layer).species_masses()
        );
    }
    assert!(after.interface_height().value() < initial.gas.interface_height().value());
    assert_admissible(network.rooms());
}

#[test]
fn cancelled_pressure_trials_publish_nothing_and_workspace_reuse_is_deterministic() {
    let mut baseline = communicating_network();
    let sources = [[LayerRates::default(); 2]; 2];
    let mut calls = 0;
    let expected = baseline
        .step(Time::new(0.02), &sources, controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 10);
    for stop in [1, calls / 2, calls] {
        let mut network = communicating_network();
        let original = network.rooms().to_vec();
        let mut seen = 0;
        let refused = network.step(Time::new(0.02), &sources, controls(), || {
            seen += 1;
            seen == stop
        });
        assert!(matches!(refused, Err(NetworkError::Cancelled)));
        assert_eq!(network.rooms(), original);
        let actual = network
            .step(Time::new(0.02), &sources, controls(), || false)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(network.rooms(), baseline.rooms());
    }
}

#[test]
fn invalid_topology_opening_and_source_depletion_leave_public_state_unchanged() {
    let mut network = communicating_network();
    let original = network.rooms().to_vec();
    let invalid_opening = WallOpening::new(
        Length::new(-0.1),
        Length::new(2.5),
        Length::new(0.2),
        DischargeCoefficient::sharp_edged().unwrap(),
    )
    .unwrap();
    assert!(network.set_wall_opening(0, invalid_opening).is_err());
    assert_eq!(network.rooms(), original);
    let mut rates = [[LayerRates::default(); 2]; 2];
    rates[0][0].species_kg_s[0] = -1e12;
    let control = NetworkControl {
        maximum_halvings: 0,
        ..controls()
    };
    assert!(
        network
            .step(Time::new(1.0), &rates, control, || false)
            .is_err()
    );
    assert_eq!(network.rooms(), original);
    network
        .step(
            Time::new(0.02),
            &[[LayerRates::default(); 2]; 2],
            controls(),
            || false,
        )
        .unwrap();
    assert_admissible(network.rooms());
    assert!(
        ZoneNetwork::new(
            vec![original[0].clone()],
            ambient(),
            vec![WallVent {
                left: Endpoint::Room(0),
                right: Endpoint::Room(0),
                opening: opening(0.1),
            }],
            vec![],
            Acceleration::new(9.81),
            || false
        )
        .is_err()
    );
}

#[test]
fn wall_candidates_follow_retried_duration_and_close_joint_gas_wall_energy() {
    use super::wall::{
        GasContact, HeatTransferCoefficient, SpecificEnthalpy, ThermalConductivity, Wall1d,
        WallBoundary, WallGeometry, WallSolve,
    };
    use fs_qty::Density;
    use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

    struct WallLaw;
    impl ThermalLaw for WallLaw {
        fn temperature_k(&self, enthalpy: f64) -> Result<f64, DepositionError> {
            let temperature = 300.0 + enthalpy / 1_000.0;
            if !(250.0..=2_000.0).contains(&temperature) {
                return Err(DepositionError::MaterialDomain);
            }
            Ok(temperature)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1e-3
        }
    }
    let law = WallLaw;
    let wall = Wall1d::new(
        &law,
        WallGeometry {
            area: Area::new(2.0),
            thickness: Length::new(0.04),
            density: Density::new(500.0),
            conductivity: ThermalConductivity::new(0.2),
            cells: 4,
        },
        SpecificEnthalpy::new(0.0),
        || false,
    )
    .unwrap();
    let original_cells = wall.simulation().cells().clone();
    let original_wall_energy = wall.simulation().energy_j();
    let mut workspace = wall.workspace().unwrap();
    let initial = room(
        gas(false),
        100_000.0,
        [400.0, 600.0],
        [[0.8, 0.2], [0.2, 0.8]],
    );
    // Internal transfer retains the identical donor species and enthalpy on
    // both sides. At the requested one-second duration it depletes the donor,
    // forcing the network to discard a successfully solved wall candidate.
    let transfer_rate = MassFlowRate::new(1.5 * initial.gas.layer(Layer::Lower).mass_kg());
    let mut network = ZoneNetwork::new(
        vec![initial],
        ambient(),
        vec![],
        vec![],
        Acceleration::new(9.81),
        || false,
    )
    .unwrap();
    let before = network.rooms().to_vec();
    let (initial_species, initial_gas_energy) = inventories(&before);
    let mut attempted_durations = Vec::new();
    let (report, candidate) = network
        .step_with_sources(
            Time::new(1.0),
            controls(),
            |accepted, duration, sources, cancelled| {
                attempted_durations.push(duration.value());
                assert_eq!(accepted, before);
                assert_eq!(wall.time(), Time::new(0.0));
                assert_eq!(wall.simulation().cells(), &original_cells);
                let gas = &accepted[0].gas;
                let lower_area = wall.geometry().area.value() * gas.interface_height().value()
                    / gas.height().value();
                let areas = [lower_area, wall.geometry().area.value() - lower_area];
                let contacts = [Layer::Lower, Layer::Upper].map(|layer| {
                    let i = layer.index();
                    GasContact {
                        area: Area::new(areas[i]),
                        temperature: gas.layer(layer).temperature(),
                        htc: HeatTransferCoefficient::new([10.0, 20.0][i]),
                    }
                });
                let candidate = wall.candidate(
                    duration,
                    WallBoundary {
                        gas: contacts,
                        ambient_temperature: Temperature::new(280.0),
                        ambient_htc: HeatTransferCoefficient::new(8.0),
                        external_radiant_power: Power::new(200.0),
                    },
                    WallSolve::default(),
                    &mut workspace,
                    cancelled,
                )?;
                let carried = gas.carried_rate(Layer::Lower, transfer_rate)?;
                sources[0] = [carried.opposite(), carried];
                for (rate, power) in sources[0].iter_mut().zip(candidate.report.gas_heat_power) {
                    rate.add_heat(power)?;
                }
                Ok(candidate)
            },
            || false,
        )
        .unwrap();

    assert!(report.halvings > 0);
    assert_eq!(attempted_durations.len(), report.halvings + 1);
    assert_eq!(attempted_durations[0], 1.0);
    for pair in attempted_durations.windows(2) {
        assert_eq!(pair[1], 0.5 * pair[0]);
    }
    assert_eq!(*attempted_durations.last().unwrap(), report.elapsed.value());
    assert!(report.iterations > 1);
    assert_eq!(report.maximum_split_pressure_change_pa, 0.0);
    assert_eq!(candidate.report.elapsed, report.elapsed);
    assert_eq!(candidate.wall.time(), report.elapsed);
    assert_eq!(network.rooms()[0].gas.time(), report.elapsed);
    assert_eq!(wall.time(), Time::new(0.0));
    assert_eq!(wall.simulation().cells(), &original_cells);
    assert_eq!(
        candidate.wall.simulation().mass_kg(),
        wall.simulation().mass_kg()
    );
    assert_balance(&report);
    assert_admissible(network.rooms());

    let (final_species, final_gas_energy) = inventories(network.rooms());
    for i in 0..2 {
        assert!((final_species[i] - initial_species[i]).abs() < 2e-11);
        assert_eq!(report.external_species_change_kg[i], 0.0);
    }
    let gas_exchange = candidate
        .report
        .gas_heat_energy
        .iter()
        .map(|energy| energy.value())
        .sum::<f64>();
    assert!(gas_exchange < 0.0);
    assert!((report.external_physical_energy_j - gas_exchange).abs() < 1e-7);
    assert!((final_gas_energy - initial_gas_energy - gas_exchange).abs() < 1e-4);
    assert_eq!(
        candidate.report.external_radiant_energy.value(),
        200.0 * report.elapsed.value()
    );
    assert!(candidate.report.ambient_heat_energy.value() < 0.0);
    let combined_residual = final_gas_energy - initial_gas_energy
        + candidate.wall.simulation().energy_j()
        - original_wall_energy
        - candidate.report.ambient_heat_energy.value()
        - candidate.report.external_radiant_energy.value();
    assert!(
        combined_residual.abs() < 1e-4,
        "joint energy residual: {combined_residual:e}"
    );
}
