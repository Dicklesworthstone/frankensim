use super::*;
use crate::compartment::radiation::walls::closed_room::*;
use crate::compartment::{Endpoint, WallVent};
use crate::stratified::{self, DischargeCoefficient, TwoLayerColumn, WallOpening};
use fs_gas::multicomponent::compartment::{LayerRates, ZoneSolve};
use fs_qty::MassFlowRate;

fn room(temperatures: [f64; 2], interface: f64) -> TwoZone {
    let (gas, fractions) = air::source().unwrap();
    TwoZone::new(
        gas,
        Area::new(1.0),
        Length::new(1.0),
        Length::new(interface),
        Pressure::new(101325.0),
        temperatures.map(Temperature::new),
        [fractions.clone(), fractions],
        || false,
    )
    .unwrap()
}

fn model(room: &TwoZone) -> ParticipatingEnclosure {
    ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(1.0),
            segments: vec![
                OpticalSegment {
                    layer: Layer::Lower,
                    length: room.interface_height(),
                },
                OpticalSegment {
                    layer: Layer::Upper,
                    length: Length::new(room.height().value() - room.interface_height().value()),
                },
            ],
        }],
        "normal-chord fixture rebuilt from actual endpoint layer depths",
        1e-12,
        || false,
    )
    .unwrap()
}

fn run(
    original: &TwoZone,
    walls: &[Wall1d<'static>; 2],
    dt: f64,
    kappa: [f64; 2],
    h: f64,
    controls: ClosedRoomWallSolve,
    cancel: impl FnMut() -> bool,
) -> Result<ClosedRoomWallCandidate<'static, 'static>, ClosedRoomWallError> {
    run_mode(original, walls, dt, kappa, h, controls, false, cancel)
}

#[allow(clippy::too_many_arguments)]
fn run_mode(
    original: &TwoZone,
    walls: &[Wall1d<'static>; 2],
    dt: f64,
    kappa: [f64; 2],
    h: f64,
    controls: ClosedRoomWallSolve,
    split: bool,
    cancel: impl FnMut() -> bool,
) -> Result<ClosedRoomWallCandidate<'static, 'static>, ClosedRoomWallError> {
    let mut optical_work = model(original).wall_workspace(walls).unwrap();
    let absorption = absorption(kappa);
    let mut gas_work = ClosedRoomWallWorkspace::new(original).unwrap();
    let evaluate = |endpoint: &TwoZone, cancelled: &mut dyn FnMut() -> bool| {
        let boundaries = [0, 1].map(|_| {
            let mut b = insulated();
            b.gas[0] = GasContact {
                area: Area::new(if endpoint.layer(Layer::Lower).mass_kg() > 0.0 {
                    1.0
                } else {
                    0.0
                }),
                temperature: endpoint.layer(Layer::Lower).temperature(),
                htc: HeatTransferCoefficient::new(h),
            };
            b
        });
        model(endpoint).candidate_walls(
            walls,
            Time::new(dt),
            &boundaries,
            [
                endpoint.layer(Layer::Lower).temperature(),
                endpoint.layer(Layer::Upper).temperature(),
            ],
            WallOptics::Gray {
                absorption: &absorption,
                solve: RadiationSolve::default(),
            },
            |_, t| Ok(emissivity(1.0, t)),
            None,
            EnclosureWallSolve::default(),
            &mut optical_work,
            cancelled,
        )
    };
    if split {
        candidate_split_room_walls(
            original,
            Time::new(dt),
            None,
            controls,
            &mut gas_work,
            evaluate,
            cancel,
        )
    } else {
        candidate_closed_room_walls(
            original,
            Time::new(dt),
            None,
            controls,
            &mut gas_work,
            evaluate,
            cancel,
        )
    }
}

// Independent scalar backward-Euler energy equation, using the original
// constitutive caloric curve but neither TwoZone::advance nor a radiation solve.
fn single_layer_oracle(original: &TwoZone, dt: f64, kappa: f64, h: f64) -> (f64, f64, f64) {
    let initial = original.layer(Layer::Lower).temperature();
    let energy_change = |temperature: f64| {
        original
            .gas()
            .source()
            .models()
            .iter()
            .zip(original.layer(Layer::Lower).species_masses())
            .map(|(model, mass)| {
                mass / model.molar_mass().value()
                    * model
                        .molar_internal_energy_change(initial, Temperature::new(temperature))
                        .unwrap()
                        .value()
            })
            .sum::<f64>()
    };
    let response = dt / 10_000.0 + 1.0 / 100.0;
    let (mut lo, mut hi) = (300.0, initial.value());
    for _ in 0..90 {
        let temperature = 0.5 * (lo + hi);
        let p = -energy_change(temperature) / (2.0 * dt);
        let face = 300.0 + p * response;
        let flux =
            -(-kappa).exp_m1() * STEFAN_BOLTZMANN_W_M2_K4 * (temperature.powi(4) - face.powi(4))
                + h * (temperature - face);
        if p - flux > 0.0 {
            lo = temperature;
        } else {
            hi = temperature;
        }
    }
    let temperature = 0.5 * (lo + hi);
    let p = -energy_change(temperature) / (2.0 * dt);
    (temperature, 300.0 + p * response, p)
}

#[test]
fn g1_stiff_finite_gas_radiation_and_convection_match_independent_energy_equations() {
    let original = room([900.0, 300.0], 1.0);
    let before = original.clone();
    let walls = [wall(300.0, 1.0, 1), wall(300.0, 1.0, 1)];
    for (kappa, h) in [(0.7, 0.0), (0.0, 50.0), (0.7, 50.0)] {
        let dt = 300.0;
        let candidate = run(
            &original,
            &walls,
            dt,
            [kappa, 0.0],
            h,
            ClosedRoomWallSolve::default(),
            || false,
        )
        .unwrap();
        let (temperature, face, p) = single_layer_oracle(&original, dt, kappa, h);
        close(
            candidate.room.layer(Layer::Lower).temperature().value(),
            temperature,
            3e-6,
        );
        for w in &candidate.walls.walls {
            close(w.report.surface_temperatures[0].value(), face, 3e-6);
            close(w.report.audit.stored_energy_change_j, dt * p, 1e-4);
            assert_eq!(w.wall.time(), Time::new(dt));
        }
        close(
            candidate.walls.report.gas_heat_power[0].value(),
            -2.0 * p,
            1e-5,
        );
        assert_eq!(candidate.walls.report.gas_heat_power[1].value(), 0.0);
        assert!(candidate.report.combined_energy_residual_j.abs() < 1e-4);
        assert!(
            candidate
                .report
                .endpoint_power_residual_w
                .iter()
                .all(|v| v.abs() < 1e-4)
        );
        assert_eq!(
            candidate.room.layer(Layer::Lower).species_masses(),
            original.layer(Layer::Lower).species_masses()
        );
        assert_eq!(original, before);
        assert_eq!(walls[0].time(), Time::new(0.0));
    }
    // Holding the initial gas spectrum for this same interval requests more
    // energy than the finite inventory contains. The new driver takes the
    // full requested interval using endpoint temperatures instead.
    let frozen = model(&original)
        .candidate_walls(
            &walls,
            Time::new(300.0),
            &[insulated(); 2],
            [Temperature::new(900.0), Temperature::new(300.0)],
            WallOptics::Gray {
                absorption: &absorption([0.7, 0.0]),
                solve: RadiationSolve::default(),
            },
            |_, t| Ok(emissivity(1.0, t)),
            None,
            EnclosureWallSolve::default(),
            &mut model(&original).wall_workspace(&walls).unwrap(),
            || false,
        )
        .unwrap();
    let mut rates = [LayerRates::default(); 2];
    frozen.report.add_gas_heat(&mut rates).unwrap();
    assert!(
        original
            .clone()
            .advance(Time::new(300.0), rates, ZoneSolve::default(), || false)
            .is_err()
    );
}

#[test]
fn g3_two_layer_pressure_work_and_geometry_feedback_refine_at_fixed_inventory() {
    let start = room([550.0, 900.0], 0.5);
    let evolve = |steps: usize| {
        let mut room = start.clone();
        let mut walls = [wall(300.0, 1.0, 1), wall(350.0, 1.0, 1)];
        let initial = room.sensible_internal_energy().value()
            + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
        for _ in 0..steps {
            let candidate = run(
                &room,
                &walls,
                8.0 / steps as f64,
                [0.2, 0.7],
                8.0,
                ClosedRoomWallSolve::default(),
                || false,
            )
            .unwrap();
            room = candidate.room;
            for (wall, next) in walls.iter_mut().zip(candidate.walls.walls) {
                *wall = next.wall;
            }
            assert!(candidate.report.combined_energy_residual_j.abs() < 1e-4);
        }
        let final_energy = room.sensible_internal_energy().value()
            + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
        close(final_energy, initial, 0.003);
        for layer in [Layer::Lower, Layer::Upper] {
            assert_eq!(
                room.layer(layer).species_masses(),
                start.layer(layer).species_masses()
            );
        }
        assert!((room.interface_height().value() - start.interface_height().value()).abs() > 1e-4);
        [
            room.layer(Layer::Lower).temperature().value(),
            room.layer(Layer::Upper).temperature().value(),
        ]
    };
    let reference = evolve(32);
    let errors: Vec<f64> = [4, 8, 16]
        .into_iter()
        .map(|n| {
            evolve(n)
                .into_iter()
                .zip(reference)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max)
        })
        .collect();
    assert!(errors[1] < 0.75 * errors[0], "refinement errors {errors:?}");
    assert!(errors[2] < 0.75 * errors[1], "refinement errors {errors:?}");
}

fn network_with_opening(original: &TwoZone, width: f64) -> (ZoneNetwork, WallOpening) {
    let opening = WallOpening::new(
        Length::new(0.0),
        Length::new(1.0),
        Length::new(width),
        DischargeCoefficient::sharp_edged().unwrap(),
    )
    .unwrap();
    let state = original.layer(Layer::Lower);
    let fractions = state
        .species_masses()
        .iter()
        .map(|m| m / state.mass_kg())
        .collect();
    let network = ZoneNetwork::new(
        vec![Room {
            gas: original.clone(),
            floor_elevation: Length::new(0.0),
        }],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(300.0),
            mass_fractions: fractions,
        },
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Ambient,
            opening: opening.clone(),
        }],
        vec![],
        Acceleration::new(0.0),
        || false,
    )
    .unwrap();
    (network, opening)
}

fn apply_thermal(
    rooms: &mut [Room],
    duration: Time,
    walls: &[Wall1d<'static>; 2],
    additional: &mut LayerRates,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<
    (
        EnclosureWallCandidate<'static, 'static>,
        ClosedRoomWallReport,
    ),
    NetworkError,
> {
    let result = run_mode(
        &rooms[0].gas,
        walls,
        duration.value(),
        [0.7, 0.0],
        0.0,
        ClosedRoomWallSolve::default(),
        true,
        cancelled,
    )?;
    for power in result.walls.report.gas_heat_power {
        additional.add_heat(power)?;
    }
    rooms[0].gas = result.room;
    Ok((result.walls, result.report))
}

#[test]
fn g1_source_inclusive_network_recovers_sealed_implicit_radiation() {
    let original = room([900.0, 300.0], 1.0);
    let (mut network, _) = network_with_opening(&original, 0.0);
    let walls = [wall(300.0, 1.0, 1), wall(300.0, 1.0, 1)];
    let dt = 300.0;
    let (balance, (), (candidate, endpoint)) = network
        .step_with_pressure_sources(
            Time::new(dt),
            NetworkControl {
                maximum_halvings: 0,
                ..NetworkControl::default()
            },
            |_, _, _, _| Ok(()),
            |rooms, duration, _, additional, cancel| {
                apply_thermal(rooms, duration, &walls, additional, cancel)
            },
            || false,
        )
        .unwrap();
    let (temperature, face, _) = single_layer_oracle(&original, dt, 0.7, 0.0);
    close(
        network.rooms()[0]
            .gas
            .layer(Layer::Lower)
            .temperature()
            .value(),
        temperature,
        3e-6,
    );
    for w in &candidate.walls {
        close(w.report.surface_temperatures[0].value(), face, 3e-6);
        assert_eq!(w.wall.time(), network.rooms()[0].gas.time());
    }
    assert_eq!(network.rooms()[0].gas.time(), Time::new(dt));
    assert_eq!(balance.halvings, 0);
    assert_eq!(balance.maximum_split_pressure_change_pa, 0.0);
    assert!(balance.maximum_pressure_residual_pa <= 1e-5);
    assert!(endpoint.combined_energy_residual_j.abs() < 1e-4);
}

#[test]
fn g1_vent_transfer_uses_radiatively_cooled_pressure_and_closes_physical_energy() {
    let original = room([900.0, 300.0], 1.0);
    let (mut network, opening) = network_with_opening(&original, 0.03);
    let walls = [wall(300.0, 1.0, 1), wall(300.0, 1.0, 1)];
    let (balance, (), (candidate, _)) = network
        .step_with_pressure_sources(
            Time::new(0.05),
            NetworkControl {
                maximum_halvings: 0,
                ..NetworkControl::default()
            },
            |_, _, _, _| Ok(()),
            |rooms, duration, _, additional, cancel| {
                apply_thermal(rooms, duration, &walls, additional, cancel)
            },
            || false,
        )
        .unwrap();
    let endpoint = &network.rooms()[0].gas;
    assert!(endpoint.pressure().value() < network.ambient().pressure.value());
    assert!(endpoint.layer(Layer::Lower).temperature().value() < 900.0);
    assert!(endpoint.layer(Layer::Lower).mass_kg() > original.layer(Layer::Lower).mass_kg());
    assert_eq!(endpoint.time(), balance.elapsed);
    assert!(
        candidate
            .walls
            .iter()
            .all(|w| w.wall.time() == endpoint.time())
    );
    assert_eq!(balance.maximum_split_pressure_change_pa, 0.0);
    assert!(balance.maximum_pressure_residual_pa <= 1e-5);

    // Independently reconstruct opening transfer at the accepted cooled pressure,
    // using the vent owner's declared frozen-start donor temperatures. At start
    // both pressures are equal and gravity is zero: omitting radiation from the
    // pressure equation would predict exactly zero vent mass transfer.
    let initial_density = original.layer(Layer::Lower).density_kg_m3().unwrap();
    let inside_density =
        initial_density * endpoint.pressure().value() / original.pressure().value();
    let outside_density = initial_density * 900.0 / 300.0;
    let inside = TwoLayerColumn::new(
        Length::new(0.0),
        endpoint.pressure(),
        Length::new(1.0),
        Density::new(inside_density),
        Density::new(inside_density),
    )
    .unwrap();
    let outside = TwoLayerColumn::new(
        Length::new(0.0),
        network.ambient().pressure,
        Length::new(0.0),
        Density::new(outside_density),
        Density::new(outside_density),
    )
    .unwrap();
    let flow =
        stratified::wall_vent_flow(&opening, inside, outside, Acceleration::new(0.0), || false)
            .unwrap();
    let incoming: f64 = -flow
        .slabs()
        .iter()
        .map(|s| s.mass_rate.value())
        .sum::<f64>();
    let dt = balance.elapsed.value();
    close(
        endpoint.layer(Layer::Lower).mass_kg() - original.layer(Layer::Lower).mass_kg(),
        dt * incoming,
        2e-9,
    );
    let carried = original
        .inflow(
            MassFlowRate::new(incoming),
            &network.ambient().mass_fractions,
            network.ambient().temperature,
        )
        .unwrap();
    let physical_ambient = dt
        * (carried.referenced_energy_w
            + carried
                .species_kg_s
                .iter()
                .zip(original.gas().energy_references())
                .map(|(m, e)| m * e)
                .sum::<f64>());
    let heat = dt
        * candidate
            .report
            .gas_heat_power
            .iter()
            .map(|p| p.value())
            .sum::<f64>();
    close(
        balance.external_physical_energy_j - heat,
        physical_ambient,
        1e-3,
    );
    let combined = endpoint.physical_internal_energy().value()
        - original.physical_internal_energy().value()
        + candidate.report.wall_storage_change.value()
        - physical_ambient;
    assert!(combined.abs() < 1e-3, "combined energy {combined:e}");
    assert!(balance.physical_energy_residual_j.abs() < 1e-4);
    assert!(balance.species_residual_kg.iter().all(|v| v.abs() < 1e-12));
}

#[test]
fn g4_pressure_trials_retry_and_cancel_without_publishing_gas_or_walls() {
    use std::cell::Cell;
    let original = room([900.0, 300.0], 1.0);
    let walls = [wall(300.0, 1.0, 1), wall(300.0, 1.0, 1)];
    let (mut network, _) = network_with_opening(&original, 0.03);
    let before = network.rooms().to_vec();
    let cancel_after_thermal = Cell::new(false);
    let cancelled = network.step_with_pressure_sources(
        Time::new(0.08),
        NetworkControl::default(),
        |_, _, _, _| Ok(()),
        |rooms, dt, _, additional, cancel| {
            let result = apply_thermal(rooms, dt, &walls, additional, cancel)?;
            cancel_after_thermal.set(true);
            Ok(result)
        },
        || cancel_after_thermal.get(),
    );
    assert!(matches!(cancelled, Err(NetworkError::Cancelled)));
    assert_eq!(network.rooms(), before);
    assert!(walls.iter().all(|w| w.time() == Time::new(0.0)));
    let (balance, (), (candidate, _)) = network
        .step_with_pressure_sources(
            Time::new(0.08),
            NetworkControl::default(),
            |_, _, _, _| Ok(()),
            |rooms, dt, _, additional, cancel| {
                let result = apply_thermal(rooms, dt, &walls, additional, cancel)?;
                if dt.value() > 0.04 {
                    return Err(NetworkError::Invalid(
                        "forced rejected thermal source trial",
                    ));
                }
                Ok(result)
            },
            || false,
        )
        .unwrap();
    assert_eq!(balance.elapsed, Time::new(0.04));
    assert_eq!(balance.halvings, 1);
    assert!(
        candidate
            .walls
            .iter()
            .all(|w| w.wall.time() == balance.elapsed)
    );
    assert_eq!(network.rooms()[0].gas.time(), balance.elapsed);
    assert!(walls.iter().all(|w| w.time() == Time::new(0.0)));
    assert_eq!(
        NetworkError::from(ClosedRoomWallError::Cancelled),
        NetworkError::Cancelled
    );
}

#[test]
fn g4_cancellation_budget_and_wrong_wall_clocks_leave_original_owners_unchanged() {
    let original = room([550.0, 900.0], 0.5);
    let before = original.clone();
    let walls = [wall(300.0, 1.0, 1), wall(300.0, 1.0, 1)];
    let mut polls = 0;
    let result = run(
        &original,
        &walls,
        10.0,
        [0.2, 0.7],
        0.0,
        ClosedRoomWallSolve::default(),
        || {
            polls += 1;
            polls >= 40
        },
    );
    assert!(matches!(result, Err(ClosedRoomWallError::Cancelled)));
    assert!(matches!(
        run(
            &original,
            &walls,
            10.0,
            [0.2, 0.7],
            0.0,
            ClosedRoomWallSolve {
                maximum_evaluations: 1,
                ..ClosedRoomWallSolve::default()
            },
            || false
        ),
        Err(ClosedRoomWallError::Budget { .. })
    ));
    let model = model(&original);
    let result = candidate_closed_room_walls(
        &original,
        Time::new(2.0),
        None,
        ClosedRoomWallSolve::default(),
        &mut ClosedRoomWallWorkspace::new(&original).unwrap(),
        |gas, cancel| {
            model.candidate_walls(
                &walls,
                Time::new(1.0),
                &[insulated(); 2],
                [
                    gas.layer(Layer::Lower).temperature(),
                    gas.layer(Layer::Upper).temperature(),
                ],
                WallOptics::Gray {
                    absorption: &absorption([0.2, 0.7]),
                    solve: RadiationSolve::default(),
                },
                |_, t| Ok(emissivity(1.0, t)),
                None,
                EnclosureWallSolve::default(),
                &mut model.wall_workspace(&walls).unwrap(),
                cancel,
            )
        },
        || false,
    );
    assert!(matches!(result, Err(ClosedRoomWallError::Invalid(_))));
    assert_eq!(original, before);
    assert!(walls.iter().all(|w| w.time() == Time::new(0.0)));
    run(
        &original,
        &walls,
        10.0,
        [0.2, 0.7],
        0.0,
        ClosedRoomWallSolve::default(),
        || false,
    )
    .unwrap();
}
