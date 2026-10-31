use super::*;
use crate::compartment::radiation::spectral::{BandSurfaceEmissivity, SpectralBand};
use crate::compartment::radiation::{AbsorptionCoefficient, OpticalPath, OpticalSegment};
use crate::compartment::wall::{
    GasContact, HeatTransferCoefficient, SpecificEnthalpy, ThermalConductivity, WallGeometry,
};
use crate::compartment::{Ambient, NetworkControl, NetworkError, Room, ZoneNetwork};
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY, ViewFactorMatrix,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Acceleration, Area, Density, Length, Pressure};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

#[path = "closed_room/tests.rs"]
mod closed_room;

// Original admitted low-temperature NASA9 air source, without a second caloric bank.
#[path = "../../../../examples/compartment_fire/air.rs"]
mod air;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / 1000.0;
        if (250.0..=2200.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
}
static LAW: Sensible = Sensible;

fn wall(t: f64, area: f64, cells: usize) -> Wall1d<'static> {
    Wall1d::new(
        &LAW,
        WallGeometry {
            area: Area::new(area),
            thickness: Length::new(0.02),
            density: Density::new(500.0),
            conductivity: ThermalConductivity::new(1.0),
            cells,
        },
        SpecificEnthalpy::new((t - 300.0) * 1000.0),
        || false,
    )
    .unwrap()
}
fn emissivity(value: f64, temperature: Temperature) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 2200.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored analytic gray boundary".into(),
                license: "internal-test".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "analytic wall".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    SurfaceEmissivity::from_card(
        "test",
        &card,
        temperature.value(),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap()
}
fn absorption(values: [f64; 2]) -> [LayerAbsorption; 2] {
    values.map(|v| {
        LayerAbsorption::new(
            AbsorptionCoefficient::new(v),
            "authored optical fixture, not air spectroscopy",
        )
        .unwrap()
    })
}
fn optics(area: f64) -> ParticipatingEnclosure {
    ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(area).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(area),
            segments: vec![
                OpticalSegment {
                    layer: Layer::Lower,
                    length: Length::new(0.5),
                },
                OpticalSegment {
                    layer: Layer::Upper,
                    length: Length::new(1.0),
                },
            ],
        }],
        "authored reciprocal normal-chord fixture",
        1e-12,
        || false,
    )
    .unwrap()
}
fn insulated() -> WallBoundary {
    WallBoundary {
        gas: [GasContact {
            area: Area::new(0.0),
            temperature: Temperature::new(300.0),
            htc: HeatTransferCoefficient::new(0.0),
        }; 2],
        ambient_temperature: Temperature::new(300.0),
        ambient_htc: HeatTransferCoefficient::new(0.0),
        external_radiant_power: Power::new(0.0),
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.15e} != {b:.15e}, tolerance={tolerance:e}"
    );
}
fn account(candidate: &EnclosureWallCandidate<'_, '_>) {
    let report = &candidate.report;
    let dt = report.elapsed.value();
    let mut direct = report
        .gas_heat_energy
        .iter()
        .map(|e| e.value())
        .sum::<f64>();
    for (i, wall) in candidate.walls.iter().enumerate() {
        direct += wall.report.audit.stored_energy_change_j
            - wall.report.ambient_heat_energy.value()
            - wall.report.external_radiant_energy.value();
        assert_eq!(
            wall.report.coupled_radiant_energy.value(),
            dt * report.radiation.surface_heat_power()[i].value()
        );
        assert_eq!(wall.wall.time(), report.elapsed);
    }
    close(direct, report.combined_energy_residual_j, 2e-9);
    assert!(report.combined_energy_residual_j.abs() < 1e-4);
    assert!(report.maximum_temperature_residual_k <= 1e-7);
    assert!(report.maximum_power_residual_w <= 1e-4);
}

#[test]
fn implicit_faces_match_independent_one_cell_black_plate_solution_when_frozen_flux_fails() {
    let model = optics(1.0);
    let walls = [wall(900.0, 1.0, 1), wall(300.0, 1.0, 1)];
    let kappa = absorption([0.0; 2]);
    let dt = 300.0;
    let frozen = STEFAN_BOLTZMANN_W_M2_K4 * (900.0_f64.powi(4) - 300.0_f64.powi(4));
    assert!(
        walls[0]
            .candidate_with_exchange(
                Time::new(dt),
                insulated(),
                Power::new(-frozen),
                WallSolve::default(),
                &mut walls[0].workspace().unwrap(),
                || false
            )
            .is_err()
    );
    let result = model
        .candidate_walls(
            &walls,
            Time::new(dt),
            &[insulated(); 2],
            [Temperature::new(300.0); 2],
            WallOptics::Gray {
                absorption: &kappa,
                solve: RadiationSolve::default(),
            },
            |_, t| Ok(emissivity(1.0, t)),
            None,
            EnclosureWallSolve::default(),
            &mut model.wall_workspace(&walls).unwrap(),
            || false,
        )
        .unwrap();
    // Independent one-cell energy and half-cell face equations:
    // T_hot_face = 900 - P(dt/C + 1/G); T_cold_face = 300 + P(dt/C + 1/G).
    let capacity = 500.0 * 0.02 * 1000.0;
    let half_conductance = 1.0 / 0.01;
    let response = dt / capacity + 1.0 / half_conductance;
    let (mut lo, mut hi) = (0.0, 300.0 / response);
    for _ in 0..100 {
        let p = 0.5 * (lo + hi);
        let balance = STEFAN_BOLTZMANN_W_M2_K4
            * ((900.0 - response * p).powi(4) - (300.0 + response * p).powi(4))
            - p;
        if balance > 0.0 {
            lo = p;
        } else {
            hi = p;
        }
    }
    let power = 0.5 * (lo + hi);
    close(
        result.report.radiation.surface_heat_power()[0].value(),
        -power,
        2e-5,
    );
    close(
        result.report.radiation.surface_heat_power()[1].value(),
        power,
        2e-5,
    );
    close(
        result.walls[0].report.surface_temperatures[0].value(),
        900.0 - response * power,
        1e-7,
    );
    close(
        result.walls[1].report.surface_temperatures[0].value(),
        300.0 + response * power,
        1e-7,
    );
    close(
        result.walls[0].wall.cell_temperature(0).unwrap().value(),
        900.0 - dt * power / capacity,
        1e-7,
    );
    assert_eq!(result.report.gas_heat_power, [Power::new(0.0); 2]);
    assert!(result.report.iterations > 1);
    assert_eq!(walls[0].time(), Time::new(0.0));
    account(&result);
}

#[test]
fn exact_equilibrium_and_isolated_wall_recovery_use_the_original_owners() {
    let model = optics(1.0);
    let walls = [wall(500.0, 1.0, 3), wall(500.0, 1.0, 2)];
    let mut boundary = insulated();
    boundary.gas[1] = GasContact {
        area: Area::new(1.0),
        temperature: Temperature::new(500.0),
        htc: HeatTransferCoefficient::new(10.0),
    };
    boundary.ambient_temperature = Temperature::new(500.0);
    boundary.ambient_htc = HeatTransferCoefficient::new(3.0);
    let kappa = absorption([0.3, 1.0]);
    let result = model
        .candidate_walls(
            &walls,
            Time::new(100.0),
            &[boundary; 2],
            [Temperature::new(500.0); 2],
            WallOptics::Gray {
                absorption: &kappa,
                solve: RadiationSolve::default(),
            },
            |_, t| Ok(emissivity(0.8, t)),
            None,
            EnclosureWallSolve::default(),
            &mut model.wall_workspace(&walls).unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(result.report.iterations, 1);
    assert_eq!(result.report.gas_heat_power, [Power::new(0.0); 2]);
    assert_eq!(result.report.wall_storage_change.value(), 0.0);
    for (i, w) in result.walls.iter().enumerate() {
        assert_eq!(w.wall.simulation().cells(), walls[i].simulation().cells());
    }
    account(&result);

    // Equal walls see no inter-wall radiation in a transparent enclosure;
    // each endpoint must recover its original isolated thermal solve.
    let walls = [wall(300.0, 1.0, 2), wall(300.0, 1.0, 2)];
    let kappa = absorption([0.0; 2]);
    let expected = walls[0]
        .candidate(
            Time::new(2.0),
            boundary,
            WallSolve::default(),
            &mut walls[0].workspace().unwrap(),
            || false,
        )
        .unwrap();
    let result = model
        .candidate_walls(
            &walls,
            Time::new(2.0),
            &[boundary; 2],
            [Temperature::new(500.0); 2],
            WallOptics::Gray {
                absorption: &kappa,
                solve: RadiationSolve::default(),
            },
            |_, t| Ok(emissivity(0.8, t)),
            None,
            EnclosureWallSolve::default(),
            &mut model.wall_workspace(&walls).unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(
        result.walls[0].wall.simulation().cells(),
        expected.wall.simulation().cells()
    );
    assert_eq!(
        result.walls[1].wall.simulation().cells(),
        expected.wall.simulation().cells()
    );
    assert_eq!(
        result.report.radiation.surface_heat_power(),
        &[Power::new(0.0); 2]
    );
    account(&result);
}

fn network() -> ZoneNetwork {
    let (gas, fractions) = air::source().unwrap();
    let room = TwoZone::new(
        gas,
        Area::new(2.0),
        Length::new(2.0),
        Length::new(1.0),
        Pressure::new(101325.0),
        [Temperature::new(400.0), Temperature::new(900.0)],
        [fractions.clone(), fractions.clone()],
        || false,
    )
    .unwrap();
    ZoneNetwork::new(
        vec![Room {
            gas: room,
            floor_elevation: Length::new(0.0),
        }],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(300.0),
            mass_fractions: fractions,
        },
        vec![],
        vec![],
        Acceleration::new(9.81),
        || false,
    )
    .unwrap()
}
fn contacts(room: &TwoZone, external: f64) -> WallBoundary {
    WallBoundary {
        gas: [Layer::Lower, Layer::Upper].map(|l| GasContact {
            area: Area::new(1.0),
            temperature: room.layer(l).temperature(),
            htc: HeatTransferCoefficient::new(8.0),
        }),
        ambient_temperature: Temperature::new(280.0),
        ambient_htc: HeatTransferCoefficient::new(3.0),
        external_radiant_power: Power::new(external),
    }
}

#[test]
fn actual_network_candidate_closes_gas_wall_energy_and_retries_all_owners_together() {
    let model = optics(2.0);
    let walls = [wall(700.0, 2.0, 3), wall(300.0, 2.0, 3)];
    let kappa = absorption([0.2, 0.6]);
    let mut network = network();
    let initial = network.rooms()[0].gas.physical_internal_energy().value()
        + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
    let initial_mass: Vec<_> = [Layer::Lower, Layer::Upper]
        .iter()
        .map(|&l| network.rooms()[0].gas.layer(l).species_masses().to_vec())
        .collect();
    let mut work = model.wall_workspace(&walls).unwrap();
    let mut calls = 0;
    let (gas_report, candidate) = network
        .step_with_sources(
            Time::new(0.1),
            NetworkControl::default(),
            |rooms, dt, rates, cancel| {
                let candidate = model.candidate_room_walls(
                    &rooms[0].gas,
                    &walls,
                    dt,
                    &[contacts(&rooms[0].gas, 25.0), contacts(&rooms[0].gas, 40.0)],
                    WallOptics::Gray {
                        absorption: &kappa,
                        solve: RadiationSolve::default(),
                    },
                    |i, t| Ok(emissivity(if i == 0 { 0.7 } else { 0.9 }, t)),
                    None,
                    EnclosureWallSolve::default(),
                    &mut work,
                    &mut *cancel,
                )?;
                calls += 1;
                candidate.report.add_gas_heat(&mut rates[0])?;
                if calls == 1 {
                    return Err(NetworkError::Invalid(
                        "force one complete coupled-wall retry",
                    ));
                }
                Ok(candidate)
            },
            || false,
        )
        .unwrap();
    assert_eq!(gas_report.halvings, 1);
    assert_eq!(gas_report.elapsed, Time::new(0.05));
    assert!(calls >= 2);
    let final_energy = network.rooms()[0].gas.physical_internal_energy().value()
        + candidate
            .walls
            .iter()
            .map(|w| w.wall.simulation().energy_j())
            .sum::<f64>();
    let external = candidate.report.ambient_heat_energy.value()
        + candidate.report.external_radiant_energy.value();
    close(final_energy - initial, external, 2e-5);
    for (i, l) in [Layer::Lower, Layer::Upper].into_iter().enumerate() {
        for (a, b) in network.rooms()[0]
            .gas
            .layer(l)
            .species_masses()
            .iter()
            .zip(&initial_mass[i])
        {
            close(*a, *b, 1e-13);
        }
    }
    for w in &candidate.walls {
        assert_eq!(w.wall.time(), network.rooms()[0].gas.time());
    }
    assert_eq!(walls[0].time(), Time::new(0.0));
    account(&candidate);
    // The identical accepted duration/source state must reproduce its wall result.
    let direct = model
        .candidate_room_walls(
            &self::network().rooms()[0].gas,
            &walls,
            gas_report.elapsed,
            &[
                contacts(&self::network().rooms()[0].gas, 25.0),
                contacts(&self::network().rooms()[0].gas, 40.0),
            ],
            WallOptics::Gray {
                absorption: &kappa,
                solve: RadiationSolve::default(),
            },
            |i, t| Ok(emissivity(if i == 0 { 0.7 } else { 0.9 }, t)),
            None,
            EnclosureWallSolve::default(),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(candidate.report, direct.report);
    for (a, b) in candidate.walls.iter().zip(&direct.walls) {
        assert_eq!(a.wall.simulation().cells(), b.wall.simulation().cells());
    }
}

#[test]
fn spectral_feedback_retains_original_reports_and_endpoint_checks() {
    let model = optics(1.0);
    let walls = [wall(650.0, 1.0, 1), wall(300.0, 1.0, 1)];
    let band = SpectralOptics::new(
        vec![
            SpectralBand {
                wavelengths: [Length::new(1e-7), Length::new(5e-6)],
                absorption: absorption([0.0, 0.7]),
            },
            SpectralBand {
                wavelengths: [Length::new(5e-6), Length::new(1e-3)],
                absorption: absorption([0.0, 0.1]),
            },
        ],
        [None, None],
        "authored wavelength-independent surface emissivities",
    )
    .unwrap();
    let solve = SpectralSolve {
        heat_tolerance_w: 1e-6,
        radiosity_tolerance_w_m2: 1e-6,
        ..SpectralSolve::default()
    };
    let result = model
        .candidate_walls(
            &walls,
            Time::new(5.0),
            &[insulated(); 2],
            [Temperature::new(300.0), Temperature::new(900.0)],
            WallOptics::Spectral {
                optics: &band,
                solve,
            },
            |_, t| Ok(emissivity(0.8, t)),
            None,
            EnclosureWallSolve::default(),
            &mut model.wall_workspace(&walls).unwrap(),
            || false,
        )
        .unwrap();
    let WallRadiation::Spectral(report) = &result.report.radiation else {
        panic!("spectral result lost");
    };
    assert_eq!(report.optics, band);
    assert!(report.wavelength_evaluations > 0);
    assert!(result.report.work.wavelength_evaluations > report.wavelength_evaluations);
    assert!(result.report.work.spectral_factorizations > report.factorizations);
    assert!(result.report.gas_heat_power[1].value() < 0.0);
    for (source, t) in report
        .surface_emissivities
        .iter()
        .zip(&result.report.radiation_surface_temperatures)
    {
        assert_eq!(source.temperature_k(), t.value());
    }
    account(&result);
}

#[test]
fn selective_surface_bands_heat_finite_walls_without_dummy_gray_cards() {
    let model = optics(1.0);
    let walls = [wall(900.0, 1.0, 1), wall(300.0, 1.0, 1)];
    let bands = SpectralOptics::new_with_surface_bands(
        vec![
            SpectralBand {
                wavelengths: [Length::new(1e-7), Length::new(5e-6)],
                absorption: absorption([0.0; 2]),
            },
            SpectralBand {
                wavelengths: [Length::new(5e-6), Length::new(1e-3)],
                absorption: absorption([0.0; 2]),
            },
        ],
        [None; 2],
        [[0.9, 0.2], [0.1, 0.8]]
            .map(|row| {
                row.map(|value| {
                    BandSurfaceEmissivity::new(
                        value,
                        [Temperature::new(250.0), Temperature::new(2200.0)],
                        "authored selective opaque surface fixture",
                    )
                    .unwrap()
                })
                .to_vec()
            })
            .to_vec(),
    )
    .unwrap();
    let result = model
        .candidate_walls(
            &walls,
            Time::new(5.0),
            &[insulated(); 2],
            [Temperature::new(300.0); 2],
            WallOptics::Spectral {
                optics: &bands,
                solve: SpectralSolve::default(),
            },
            |_, _| panic!("band sources must not query an unrelated total emissivity"),
            None,
            EnclosureWallSolve::default(),
            &mut model.wall_workspace(&walls).unwrap(),
            || false,
        )
        .unwrap();
    let WallRadiation::Spectral(report) = &result.report.radiation else {
        panic!("spectral result lost");
    };
    assert!(report.surface_emissivities.is_empty());
    assert_eq!(report.optics, bands);
    assert!(report.surface_heat_power[1].value() > 100.0);
    assert!(result.report.radiation_surface_temperatures[0].value() < 900.0);
    assert!(result.report.radiation_surface_temperatures[1].value() > 300.0);
    assert_eq!(report.gas_heat_power, [Power::new(0.0); 2]);
    account(&result);
}

#[test]
fn cancellation_budgets_and_bad_ownership_leave_inputs_unchanged_and_work_reusable() {
    let model = optics(1.0);
    let walls = [wall(800.0, 1.0, 2), wall(300.0, 1.0, 2)];
    let kappa = absorption([0.2, 0.5]);
    let mode = WallOptics::Gray {
        absorption: &kappa,
        solve: RadiationSolve::default(),
    };
    let mut work = model.wall_workspace(&walls).unwrap();
    let mut count = 0;
    let cancel = std::cell::Cell::new(false);
    let failed = model.candidate_walls(
        &walls,
        Time::new(1.0),
        &[insulated(); 2],
        [Temperature::new(500.0), Temperature::new(1000.0)],
        mode,
        |_, t| {
            count += 1;
            if count == 3 {
                cancel.set(true);
            }
            Ok(emissivity(0.8, t))
        },
        None,
        EnclosureWallSolve::default(),
        &mut work,
        || cancel.get(),
    );
    assert!(matches!(failed, Err(EnclosureWallError::Cancelled)));
    assert_eq!(walls[0].time(), Time::new(0.0));
    let controls = EnclosureWallSolve {
        maximum_wall_evaluations: 1,
        ..EnclosureWallSolve::default()
    };
    let failed = model.candidate_walls(
        &walls,
        Time::new(1.0),
        &[insulated(); 2],
        [Temperature::new(500.0), Temperature::new(1000.0)],
        mode,
        |_, t| Ok(emissivity(0.8, t)),
        None,
        controls,
        &mut work,
        || false,
    );
    assert!(matches!(
        failed,
        Err(EnclosureWallError::Budget { owner: "wall", .. })
    ));
    let mut bad = insulated();
    bad.gas[0] = GasContact {
        area: Area::new(1.0),
        temperature: Temperature::new(300.0),
        htc: HeatTransferCoefficient::new(1.0),
    };
    assert!(matches!(
        model.candidate_walls(
            &walls,
            Time::new(1.0),
            &[bad; 2],
            [Temperature::new(500.0), Temperature::new(1000.0)],
            mode,
            |_, t| Ok(emissivity(0.8, t)),
            None,
            EnclosureWallSolve::default(),
            &mut work,
            || false
        ),
        Err(EnclosureWallError::Invalid(
            "convective and radiative gas temperatures differ"
        ))
    ));
    assert!(
        model
            .wall_workspace(&[wall(300.0, 2.0, 1), wall(300.0, 1.0, 1)])
            .is_err()
    );
    let run = |work: &mut EnclosureWallWorkspace| {
        model
            .candidate_walls(
                &walls,
                Time::new(1.0),
                &[insulated(); 2],
                [Temperature::new(500.0), Temperature::new(1000.0)],
                mode,
                |_, t| Ok(emissivity(0.8, t)),
                None,
                EnclosureWallSolve::default(),
                work,
                || false,
            )
            .unwrap()
    };
    let mut result = run(&mut work);
    let fresh = run(&mut model.wall_workspace(&walls).unwrap());
    assert_eq!(result.report, fresh.report);
    for (a, b) in result.walls.iter().zip(&fresh.walls) {
        assert_eq!(a.wall.simulation().cells(), b.wall.simulation().cells());
    }
    let mut rates = [LayerRates::default(); 2];
    result.report.add_gas_heat(&mut rates).unwrap();
    for (r, q) in rates.iter().zip(result.report.gas_heat_power) {
        assert_eq!(r.referenced_energy_w, q.value());
    }
    let saved = rates;
    result.report.gas_heat_power[1] = Power::new(f64::INFINITY);
    assert!(result.report.add_gas_heat(&mut rates).is_err());
    assert_eq!(saved, rates);
}
