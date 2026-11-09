//! Explicit original caloric receivers for the same nonlinear interface loss.

use super::super::adaptive::{ThermalFilmAdaptiveControls, ThermalFilmQuery};
use super::*;

const WALL: CellId = CellId(90);

fn flow(b: f64, wall: f64) -> TransientCarreauFilm {
    let mut value = TransientCarreauFilm::from_profile(
        source(),
        300.0,
        SPEC,
        vec![wall + 0.5 * (b + 0.25), wall + 0.5 * (b + 0.75)],
        0.0,
    )
    .unwrap();
    value.set_wall_slip_length_m(b).unwrap();
    value
}

fn thermal<'a>(
    fluid: &'a Caloric,
    wall: Option<&'a Caloric>,
    area: f64,
) -> DepositionSimulation<'a> {
    let original = inventory(fluid, [300.0; 2], area, 0.0, 0.0);
    let mut events: Vec<_> = original
        .cells()
        .values()
        .cloned()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    let mut laws = vec![(MaterialId(7), fluid as &dyn ThermalLaw)];
    if let Some(wall) = wall {
        laws.push((MaterialId(8), wall as &dyn ThermalLaw));
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: WALL,
                segment_id: 90,
                body_id: 2,
                material_id: MaterialId(8),
                mass_kg: 0.7 * area,
                specific_h_j_kg: wall.datum,
                geometry: CellGeometry {
                    start_m: [0.0, -0.4, 0.0],
                    end_m: [0.0, 0.0, 0.0],
                    volume_m3: 0.4 * area,
                    surface_area_m2: 2.0 * area,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 4096,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn make<'a>(
    caloric: &'a Caloric,
    wall: Option<&'a Caloric>,
    b: f64,
    velocity: f64,
    area: f64,
    routing: ThermalFilmSlipHeat,
) -> ThermalCarreauFilm<'a> {
    ThermalCarreauFilm::new_with_wall_slip_heat(
        flow(b, velocity),
        thermal(caloric, wall, area),
        &IDS,
        area,
        [300.0; 2],
        routing,
        || false,
    )
    .unwrap()
}

fn rate_at_stress(temperature: f64, stress: f64) -> f64 {
    // Independent monotone bisection using the original law, not production's
    // inverse or cached film stresses.
    let (mut lo, mut hi) = (0.0, 20.0);
    for _ in 0..80 {
        let rate = 0.5 * (lo + hi);
        let mu = source()
            .evaluate(temperature, rate)
            .unwrap()
            .dynamic_viscosity_pa_s();
        if mu * rate < stress.abs() {
            lo = rate;
        } else {
            hi = rate;
        }
    }
    (0.5 * (lo + hi)).copysign(stress)
}

#[test]
fn native_partitions_match_independent_local_nonlinear_heat_and_do_not_double_deposit() {
    let caloric = Caloric {
        cp: 0.05,
        datum: 1000.0,
    };
    let (area, b, wall, dt) = (2.5, 0.4, -0.1, 0.04);
    for routing in [
        ThermalFilmSlipHeat {
            fluid_fraction: 0.0,
            remainder: ThermalFilmSlipHeatReceiver::Surroundings,
        },
        ThermalFilmSlipHeat {
            fluid_fraction: 0.35,
            remainder: ThermalFilmSlipHeatReceiver::Surroundings,
        },
        ThermalFilmSlipHeat {
            fluid_fraction: 1.0,
            remainder: ThermalFilmSlipHeatReceiver::Surroundings,
        },
        ThermalFilmSlipHeat {
            fluid_fraction: 0.35,
            remainder: ThermalFilmSlipHeatReceiver::Cell(IDS[0]),
        },
    ] {
        let mut model = make(&caloric, None, b, wall, area, routing);
        let initial = model.flow().velocities_m_s().to_vec();
        let drive = FilmDrive {
            wall_velocity_m_s: wall,
            acceleration_m_s2: 0.05,
            ..drive()
        };
        let report = model
            .advance(
                dt,
                drive,
                LIMITS,
                TransientCarreauControls::default(),
                controls(),
                &[0.0; 2],
                insulated,
                || false,
            )
            .unwrap();
        let t = IDS.map(|id| model.thermal().temperature_k(id).unwrap());
        let mid: Vec<_> = initial
            .iter()
            .zip(model.flow().velocities_m_s())
            .map(|(u, v)| 0.5 * (u + v))
            .collect();
        let rates = [
            (mid[0] - wall) / (0.25 + b),
            (mid[1] - mid[0]) / 0.5,
            rate_at_stress(t[1], drive.surface_shear_pa),
        ];
        let temperatures = [t[0], 0.5 * (t[0] + t[1]), t[1]];
        let stress: [f64; 3] = std::array::from_fn(|i| {
            source()
                .evaluate(temperatures[i], rates[i].abs())
                .unwrap()
                .dynamic_viscosity_pa_s()
                * rates[i]
        });
        let bulk: [f64; 2] = std::array::from_fn(|i| {
            area * dt * 0.25 * (stress[i] * rates[i] + stress[i + 1] * rates[i + 1])
        });
        let slip = area * dt * stress[0] * b * rates[0];
        let into_fluid = routing.fluid_fraction * slip;
        let into_receiver = if matches!(routing.remainder, ThermalFilmSlipHeatReceiver::Cell(_)) {
            (1.0 - routing.fluid_fraction) * slip
        } else {
            0.0
        };
        for i in 0..2 {
            let expected = bulk[i]
                + if i == 0 {
                    into_fluid + into_receiver
                } else {
                    0.0
                };
            near(report.deposited_cell_heat_j[i], expected, 3e-12);
            near(0.5 * area * caloric.cp * (t[i] - 300.0), expected, 5e-10);
        }
        near(report.deposited_viscous_heat_j, bulk.iter().sum(), 3e-12);
        near(report.wall_slip_heat_j, slip, 3e-12);
        near(report.deposited_fluid_slip_heat_j, into_fluid, 3e-12);
        near(report.deposited_receiver_slip_heat_j, into_receiver, 3e-12);
        near(
            report.outward_wall_slip_heat_j,
            slip - into_fluid - into_receiver,
            3e-12,
        );
        near(report.combined_energy_residual_j, 0.0, 5e-10);
        near(report.wall_slip_partition_residual_j, 0.0, 3e-14);
        near(
            report.flow.flow.midpoint_wall_fluid_velocity_m_s,
            wall + b * rates[0],
            3e-12,
        );
        assert_eq!(model.flow().wall_slip_length_m(), b);
        assert_eq!(model.wall_slip_heat_routing(), routing);
        assert_eq!(model.flow().material().law(), &source());
    }
}

#[test]
fn existing_finite_wall_receives_the_complement_with_its_own_caloric_storage() {
    let fluid = Caloric {
        cp: 0.05,
        datum: 1000.0,
    };
    let wall = Caloric {
        cp: 2.0,
        datum: 500.0,
    };
    let (area, dt, power) = (2.5, 0.04, 0.03);
    let routing = ThermalFilmSlipHeat {
        fluid_fraction: 0.25,
        remainder: ThermalFilmSlipHeatReceiver::Cell(WALL),
    };
    let mut model = make(&fluid, Some(&wall), 0.4, -0.1, area, routing);
    let original_wall = model.thermal().cells()[&WALL].clone();
    let report = model
        .advance(
            dt,
            FilmDrive {
                wall_velocity_m_s: -0.1,
                ..drive()
            },
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 3],
            |candidate, cells| {
                let mut exchange = insulated(candidate, cells)?;
                exchange.heat.bulk.sources_w.insert(WALL, power);
                Ok::<_, DepositionError>(exchange)
            },
            || false,
        )
        .unwrap();
    let accepted = &model.thermal().cells()[&WALL];
    assert_eq!(accepted.material_id, original_wall.material_id);
    assert_eq!(accepted.mass_kg, original_wall.mass_kg);
    assert_eq!(accepted.geometry, original_wall.geometry);
    let stored = accepted.mass_kg * (accepted.specific_h_j_kg - original_wall.specific_h_j_kg);
    near(
        stored,
        dt * power + report.deposited_receiver_slip_heat_j,
        5e-10,
    );
    near(
        model.thermal().temperature_k(WALL).unwrap() - 300.0,
        stored / (accepted.mass_kg * wall.cp),
        5e-12,
    );
    near(
        report.deposited_receiver_slip_heat_j,
        0.75 * report.wall_slip_heat_j,
        3e-14,
    );
    near(
        report.deposited_cell_heat_j.iter().sum(),
        report.deposited_viscous_heat_j + report.deposited_fluid_slip_heat_j,
        3e-14,
    );
    near(report.other_cell_heat_j[&WALL], power * dt, 3e-14);
    near(report.other_heat_j, power * dt, 3e-14);
    assert_eq!(report.outward_wall_slip_heat_j, 0.0);
    near(
        report.thermal.stored_energy_change_j + report.kinetic_energy_change_j,
        report.other_heat_j + report.mechanical_work_j,
        5e-10,
    );
}

#[test]
fn explicit_routing_and_zero_slip_preserve_admission_and_original_step() {
    let caloric = Caloric {
        cp: 0.05,
        datum: 1000.0,
    };
    assert!(matches!(
        ThermalCarreauFilm::new(
            flow(0.3, 0.0),
            thermal(&caloric, None, 1.0),
            &IDS,
            1.0,
            [300.0; 2],
            || false
        ),
        Err(ThermalFilmError::Invalid(_))
    ));
    for routing in [
        ThermalFilmSlipHeat {
            fluid_fraction: -0.1,
            remainder: ThermalFilmSlipHeatReceiver::Surroundings,
        },
        ThermalFilmSlipHeat {
            fluid_fraction: f64::NAN,
            remainder: ThermalFilmSlipHeatReceiver::Surroundings,
        },
        ThermalFilmSlipHeat {
            fluid_fraction: 1.0,
            remainder: ThermalFilmSlipHeatReceiver::Cell(CellId(404)),
        },
    ] {
        assert!(matches!(
            ThermalCarreauFilm::new_with_wall_slip_heat(
                flow(0.3, 0.0),
                thermal(&caloric, None, 1.0),
                &IDS,
                1.0,
                [300.0; 2],
                routing,
                || false
            ),
            Err(ThermalFilmError::Invalid(_))
        ));
    }
    let mut original = ThermalCarreauFilm::new(
        flow(0.0, 0.0),
        thermal(&caloric, None, 1.0),
        &IDS,
        1.0,
        [300.0; 2],
        || false,
    )
    .unwrap();
    let mut explicit = make(
        &caloric,
        None,
        0.0,
        0.0,
        1.0,
        ThermalFilmSlipHeat {
            fluid_fraction: 0.35,
            remainder: ThermalFilmSlipHeatReceiver::Cell(IDS[0]),
        },
    );
    let a = original
        .advance(
            0.04,
            drive(),
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 2],
            insulated,
            || false,
        )
        .unwrap();
    let b = explicit
        .advance(
            0.04,
            drive(),
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 2],
            insulated,
            || false,
        )
        .unwrap();
    assert_eq!(a.flow, b.flow);
    assert_eq!(a.deposited_cell_heat_j, b.deposited_cell_heat_j);
    assert_eq!(a.endpoint_evaluations, b.endpoint_evaluations);
    assert_eq!(
        a.momentum_material_evaluations,
        b.momentum_material_evaluations
    );
    for id in IDS {
        assert_eq!(
            original.thermal().cells()[&id].specific_h_j_kg.to_bits(),
            explicit.thermal().cells()[&id].specific_h_j_kg.to_bits()
        );
    }
    assert_eq!(b.wall_slip_heat_j, 0.0);
}

#[test]
fn adaptive_exact_state_still_resolves_interface_work_and_both_heat_destinations() {
    let caloric = Caloric {
        cp: 0.05,
        datum: 1000.0,
    };
    let (length, eta, wall) = (0.4, 0.35, 0.25);
    let routing = ThermalFilmSlipHeat {
        fluid_fraction: eta,
        remainder: ThermalFilmSlipHeatReceiver::Surroundings,
    };
    let viscosity = |time: f64| {
        source()
            .evaluate(300.0 + 10.0 * time, 0.5)
            .unwrap()
            .dynamic_viscosity_pa_s()
    };
    let forcing = |q: ThermalFilmQuery| {
        Ok::<_, DepositionError>(FilmDrive {
            wall_velocity_m_s: wall,
            surface_shear_pa: 0.5 * viscosity(q.endpoint_s),
            ..FilmDrive::default()
        })
    };
    let exchange = |q: ThermalFilmQuery, candidate: &DepositionSimulation<'_>, cells: &[CellId]| {
        let mut result = insulated(candidate, cells)?;
        for (i, &id) in cells.iter().enumerate() {
            result.heat.bulk.sources_w.insert(
                id,
                0.5 * caloric.cp * 10.0
                    - 0.125 * viscosity(q.endpoint_s)
                    - if i == 0 {
                        eta * length * 0.25 * viscosity(q.endpoint_s)
                    } else {
                        0.0
                    },
            );
        }
        Ok::<_, DepositionError>(result)
    };
    let native = ThermalFilmControls {
        maximum_endpoint_evaluations: 100_000,
        maximum_momentum_material_evaluations: 50_000_000,
        ..controls()
    };
    let accuracy = |tolerance| ThermalFilmAdaptiveControls {
        initial_dt_s: 0.08,
        minimum_dt_s: 1e-6,
        maximum_dt_s: 0.08,
        velocity_absolute_tolerance_m_s: 1e-6,
        temperature_absolute_tolerance_k: 1e-5,
        energy_absolute_tolerance_j: tolerance,
        transport_absolute_tolerance_m2: 1e-6,
        relative_tolerance: 0.0,
        window_energy_tolerance_j: 2e-8,
        ..ThermalFilmAdaptiveControls::default()
    };
    let mut loose = make(&caloric, None, length, wall, 1.0, routing);
    let mut tight = make(&caloric, None, length, wall, 1.0, routing);
    let coarse = loose
        .advance_adaptive(
            0.08,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            native,
            accuracy(5e-6),
            forcing,
            exchange,
            || false,
        )
        .unwrap();
    let fine = tight
        .advance_adaptive(
            0.08,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            native,
            accuracy(5e-8),
            forcing,
            exchange,
            || false,
        )
        .unwrap();
    let n = 2048;
    let dt = 0.08 / n as f64;
    let integral = dt / 3.0
        * (0..=n)
            .map(|i| {
                let weight = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                weight * viscosity(i as f64 * dt)
            })
            .sum::<f64>();
    assert!(fine.accepted_steps > coarse.accepted_steps && fine.accuracy_rejections > 0);
    assert!(
        (fine.wall_slip_heat_j - 0.25 * length * integral).abs()
            < (coarse.wall_slip_heat_j - 0.25 * length * integral).abs() / 3.0
    );
    near(fine.viscous_heat_j, 0.25 * integral, 5e-6);
    near(fine.wall_slip_heat_j, 0.25 * length * integral, 2e-6);
    near(
        fine.deposited_fluid_slip_heat_j,
        eta * 0.25 * length * integral,
        7e-7,
    );
    near(
        fine.outward_wall_slip_heat_j,
        (1.0 - eta) * 0.25 * length * integral,
        1.3e-6,
    );
    near(fine.wall_work_j, -0.125 * integral, 2.5e-6);
    near(fine.surface_work_j, 0.475 * integral, 9.5e-6);
    near(fine.stored_energy_change_j, 0.04, 1e-9);
    near(fine.kinetic_energy_change_j, 0.0, 1e-10);
    near(fine.combined_energy_residual_j, 0.0, 2e-8);
    near(fine.wall_slip_partition_residual_j, 0.0, 2e-12);
    for (i, id) in IDS.iter().enumerate() {
        near(tight.flow().velocities_m_s()[i], [0.575, 0.825][i], 1e-10);
        near(tight.thermal().temperature_k(*id).unwrap(), 300.8, 1e-8);
        near(
            fine.other_cell_heat_j[id] + fine.deposited_cell_heat_j[i],
            0.02,
            1e-9,
        );
    }
    assert_eq!(tight.flow().wall_slip_length_m(), length);
    assert_eq!(tight.wall_slip_heat_routing(), routing);
    assert_eq!(fine.wall_slip_heat_routing, routing);
    assert_eq!(fine.deposited_receiver_slip_heat_j, 0.0);
}

#[test]
fn late_cancellation_and_adaptive_source_refusal_keep_slip_receivers_and_both_states() {
    let fluid = Caloric {
        cp: 0.05,
        datum: 1000.0,
    };
    let wall = Caloric {
        cp: 2.0,
        datum: 500.0,
    };
    let routing = ThermalFilmSlipHeat {
        fluid_fraction: 0.3,
        remainder: ThermalFilmSlipHeatReceiver::Cell(WALL),
    };
    let make = || make(&fluid, Some(&wall), 0.4, 0.0, 1.0, routing);
    let mut reference = make();
    let mut polls = 0;
    let expected = reference
        .advance(
            0.02,
            drive(),
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 3],
            insulated,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut victim = make();
    let old_velocity = victim.flow().velocities_m_s().to_vec();
    let old_enthalpy: Vec<_> = victim
        .thermal()
        .cells()
        .values()
        .map(|c| c.specific_h_j_kg.to_bits())
        .collect();
    let mut seen = 0;
    assert!(matches!(
        victim.advance(
            0.02,
            drive(),
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 3],
            insulated,
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(ThermalFilmError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(victim.flow().velocities_m_s(), old_velocity);
    assert_eq!(
        victim
            .thermal()
            .cells()
            .values()
            .map(|c| c.specific_h_j_kg.to_bits())
            .collect::<Vec<_>>(),
        old_enthalpy
    );
    assert_eq!(victim.flow().wall_slip_length_m(), 0.4);
    assert_eq!(victim.wall_slip_heat_routing(), routing);
    let retry = victim
        .advance(
            0.02,
            drive(),
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 3],
            insulated,
            || false,
        )
        .unwrap();
    assert_eq!(retry.flow, expected.flow);
    assert_eq!(retry.deposited_cell_heat_j, expected.deposited_cell_heat_j);
    assert_eq!(
        retry.deposited_receiver_slip_heat_j,
        expected.deposited_receiver_slip_heat_j
    );
    let mut victim = make();
    let calls = Cell::new(0usize);
    let result = victim.advance_adaptive(
        0.06,
        &[],
        LIMITS,
        TransientCarreauControls::default(),
        ThermalFilmControls {
            maximum_endpoint_evaluations: 100_000,
            maximum_momentum_material_evaluations: 50_000_000,
            ..controls()
        },
        ThermalFilmAdaptiveControls {
            initial_dt_s: 0.02,
            maximum_dt_s: 0.02,
            velocity_absolute_tolerance_m_s: 1.0,
            temperature_absolute_tolerance_k: 10.0,
            energy_absolute_tolerance_j: 1.0,
            transport_absolute_tolerance_m2: 1.0,
            window_energy_tolerance_j: 2e-8,
            ..ThermalFilmAdaptiveControls::default()
        },
        |q| {
            calls.set(calls.get() + 1);
            if q.start_s >= 0.02 {
                Err(DepositionError::MaterialDomain)
            } else {
                Ok(drive())
            }
        },
        |_, candidate, cells| insulated(candidate, cells),
        || false,
    );
    assert!(result.is_err() && calls.get() > 3);
    assert_eq!(victim.flow().accepted_time_s(), 0.0);
    assert_eq!(victim.flow().velocities_m_s(), old_velocity);
    assert_eq!(
        victim
            .thermal()
            .cells()
            .values()
            .map(|c| c.specific_h_j_kg.to_bits())
            .collect::<Vec<_>>(),
        old_enthalpy
    );
    assert_eq!(victim.flow().wall_slip_length_m(), 0.4);
    assert_eq!(victim.wall_slip_heat_routing(), routing);
}
