//! Source geometry, physical work and one-clock rendering under bow migration.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringMember,
    multiple::{
        stations::BowStationMotionLimits,
        performance::bow_stations::{BowStationKnot, BowStationProgram},
        ribbon::UniformBowRibbon,
    },
};
use fs_couple::bowed_string::runtime::reciprocal::InteriorStringConnection;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::ModalAttachment;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::surfaces::geometry::{
    FrictionGeometry, FrictionGeometryConfig,
};
use fs_math::c64::C64;

fn limits() -> BowStationMotionLimits {
    BowStationMotionLimits {
        maximum_speed_m_s: 20.0,
        maximum_step_m: 0.001,
        max_update_terms: 1_048_576,
    }
}
fn program() -> BowStationProgram {
    BowStationProgram {
        knots: vec![
            BowStationKnot {
                sample: 128,
                station_fractions: vec![0.215, 0.29],
            },
            BowStationKnot {
                sample: 320,
                station_fractions: vec![0.205, 0.31],
            },
        ],
        limits: limits(),
        max_knots: 8,
        max_coordinate_values: 16,
    }
}
fn path(sample: u64) -> [f64; 2] {
    if sample < 128 {
        let t = sample as f64 / 128.0;
        [0.2 + (0.215 - 0.2) * t, 0.3 + (0.29 - 0.3) * t]
    } else if sample == 128 {
        [0.215, 0.29]
    } else if sample < 320 {
        let t = (sample - 128) as f64 / 192.0;
        [0.215 + (0.205 - 0.215) * t, 0.29 + (0.31 - 0.29) * t]
    } else {
        [0.205, 0.31]
    }
}

/// Independent sine-field displacement of a uniform material string. This does
/// not call the source projector or the bow adapter under test.
fn sine(card: &fs_couple::bowed_string::BowedStringCard, station: f64) -> Vec<f64> {
    let modal_mass = card.linear_density_kg_m * card.length_m * 0.5;
    (1..=card.mode_count)
        .map(|n| (n as f64 * core::f64::consts::PI * station).sin() / modal_mass.sqrt())
        .collect()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-12 * a.abs().max(b.abs()) + 1e-15,
        "{a:e} != {b:e}"
    );
}

#[test]
fn unchanged_source_stations_preserve_the_original_coupled_steps_bitwise() {
    with_cx(|cx, gate| {
        let mut original = reacting_surface(joint(cx, gate, 15600.0), gate);
        let mut moving = reacting_surface(joint(cx, gate, 15600.0), gate)
            .with_bow_station_motion(limits(), gate)
            .unwrap();
        let zero = vec![0.0; original.system().mode_count()];
        let drives = [drive(0.01, 0.1), drive(-0.02, 0.05)];
        for _ in 0..64 {
            let a = original.step(&zero, &drives).unwrap().clone();
            let b = moving
                .step_with_bow_stations(&zero, &drives, &[0.2, 0.3], None)
                .unwrap()
                .clone();
            assert_eq!(a, b);
            assert_eq!(states(&original), states(&moving));
        }
    });
}

#[test]
fn midpoint_material_maps_match_an_independent_sine_reference_and_close_physical_work() {
    with_cx(|cx, gate| {
        let mut actual = reacting_surface(joint(cx, gate, 15600.0), gate)
            .with_bow_station_motion(limits(), gate)
            .unwrap();
        let (_, mut reference) = reacting_surface(joint(cx, gate, 15600.0), gate).into_parts();
        let cards = actual.cards().to_vec();
        let original_modes: Vec<_> = actual
            .system()
            .components()
            .iter()
            .map(|m| m.modes().to_vec())
            .collect();
        let zero = vec![0.0; actual.system().mode_count()];
        let drives = [drive(0.3, 0.0002), drive(-0.2, 0.0001)];
        let mut previous = [0.2, 0.3];
        let mut total_work = 0.0;
        let mut total_loss = 0.0;
        for step in 0..96 {
            let end = [
                0.2 + 0.00002 * (step + 1) as f64,
                0.3 - 0.00003 * (step + 1) as f64,
            ];
            let geometry: Vec<_> = cards
                .iter()
                .enumerate()
                .map(|(i, card)| FrictionGeometry {
                    left: ModalAttachment {
                        component: i,
                        shapes: sine(card, previous[i] + 0.5 * (end[i] - previous[i])),
                    },
                    right: reference.port(i).unwrap().right.clone(),
                })
                .collect();
            let before = states(&actual);
            let frame = actual
                .step_with_bow_stations(&zero, &drives, &end, None)
                .unwrap()
                .clone();
            let expected = reference
                .step_with_friction_geometry(
                    &zero,
                    &drives,
                    &geometry,
                    None,
                    None,
                    FrictionGeometryConfig {
                        max_update_terms: limits().max_update_terms,
                    },
                    gate,
                )
                .unwrap()
                .clone();
            assert_eq!(frame, expected);
            assert_eq!(actual.bow_station_fractions(), Some(end.as_slice()));
            let after = states(&actual);
            assert_eq!(
                after,
                reference
                    .components()
                    .iter()
                    .flat_map(|m| m.states().iter().copied())
                    .collect::<Vec<_>>()
            );
            for (port, geometry) in geometry.iter().enumerate() {
                assert_eq!(actual.system().port(port).unwrap().left, geometry.left);
                let increment = |attachment: &ModalAttachment| {
                    let offset: usize = actual.system().components()[..attachment.component]
                        .iter()
                        .map(|m| m.modes().len())
                        .sum();
                    attachment
                        .shapes
                        .iter()
                        .enumerate()
                        .map(|(i, b)| {
                            b * (after[offset + i].displacement_m_sqrt_kg
                                - before[offset + i].displacement_m_sqrt_kg)
                        })
                        .sum::<f64>()
                };
                let travel =
                    increment(&geometry.left) - geometry.right.as_ref().map_or(0.0, increment);
                let report = frame.ports[port];
                let drive_travel = drives[port].speed_m_s / 48000.0;
                near(report.slip_distance_m, drive_travel - travel);
                near(report.drive_work_j, report.traction_n * drive_travel);
                near(
                    report.dissipation_j,
                    report.traction_n * (drive_travel - travel),
                );
            }
            assert_eq!(frame.obstacle_work_j, 0.0);
            total_work += frame.external_work_j + frame.drive_work_j;
            total_loss += frame.network_dissipation_j + frame.friction_dissipation_j;
            near(frame.network_energy_j + total_loss, total_work);
            previous = end;
        }
        // A HOLD must use the physical endpoint, not the preceding midpoint.
        let geometry: Vec<_> = cards
            .iter()
            .enumerate()
            .map(|(i, card)| FrictionGeometry {
                left: ModalAttachment {
                    component: i,
                    shapes: sine(card, previous[i]),
                },
                right: reference.port(i).unwrap().right.clone(),
            })
            .collect();
        let frame = actual.step(&zero, &drives).unwrap().clone();
        let expected = reference
            .step_with_friction_geometry(
                &zero,
                &drives,
                &geometry,
                None,
                None,
                FrictionGeometryConfig {
                    max_update_terms: limits().max_update_terms,
                },
                gate,
            )
            .unwrap()
            .clone();
        assert_eq!(frame, expected);
        for (port, map) in geometry.iter().enumerate() {
            assert_eq!(actual.system().port(port).unwrap().left, map.left);
        }
        assert_eq!(
            actual
                .system()
                .components()
                .iter()
                .map(|m| m.modes().to_vec())
                .collect::<Vec<_>>(),
            original_modes
        );
        assert_eq!(actual.source().bows()[0].location.station_fraction, 0.2);
        assert_eq!(actual.source().bows()[1].location.station_fraction, 0.3);
        assert!(
            actual.system().components()[2..]
                .iter()
                .flat_map(|m| m.states())
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0)
        );
    });
}

#[test]
fn refused_motion_and_cancellation_restore_geometry_state_pressure_and_retry() {
    with_cx(|cx, gate| {
        let make = || {
            joint(cx, gate, 7800.0)
                .with_bow_station_motion(limits(), gate)
                .unwrap()
        };
        let mut actual = make();
        let mut reference = make();
        let zero = vec![0.0; actual.system().mode_count()];
        let drives = [drive(0.01, 0.1), drive(-0.02, 0.05)];
        for _ in 0..8 {
            actual
                .step_with_bow_stations(&zero, &drives, &[0.2001, 0.2999], None)
                .unwrap();
            reference
                .step_with_bow_stations(&zero, &drives, &[0.2001, 0.2999], None)
                .unwrap();
        }
        let snapshot = states(&actual);
        let frame = actual.system().last_frame().unwrap().clone();
        let ports: Vec<_> = (0..2)
            .map(|i| actual.system().port(i).unwrap().clone())
            .collect();
        for end in [
            vec![0.201],
            vec![0.4, 0.3],
            vec![0.0, 0.3],
            vec![f64::NAN, 0.3],
        ] {
            assert!(
                actual
                    .step_with_bow_stations(&zero, &drives, &end, None)
                    .is_err()
            );
        }
        // This error follows staged map rebuilding, so rollback must restore
        // the generic owner's coupled port caches as well as source positions.
        assert!(
            actual
                .step_with_bow_stations(&[], &drives, &[0.2002, 0.2998], None)
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            actual
                .step_with_bow_stations_under_gate(
                    &zero,
                    &drives,
                    &[0.2002, 0.2998],
                    None,
                    &cancelled
                )
                .is_err()
        );
        assert_eq!(states(&actual), snapshot);
        assert_eq!(actual.system().last_frame(), Some(&frame));
        assert_eq!(
            actual.bow_station_fractions(),
            Some([0.2001, 0.2999].as_slice())
        );
        for (i, port) in ports.iter().enumerate() {
            assert_eq!(actual.system().port(i), Some(port));
        }
        assert_eq!(
            actual.step(&zero, &drives).unwrap(),
            reference.step(&zero, &drives).unwrap()
        );
        assert_eq!(states(&actual), states(&reference));
    });
}

#[test]
fn source_trajectory_uses_the_original_gesture_clock_and_is_callback_invariant() {
    with_cx(|cx, gate| {
        let create = || {
            joint(cx, gate, 7800.0)
                .into_scheduled(schedule(), tracks(), performance_config())
                .unwrap()
        };
        let events = create().events().to_vec();
        let mut direct = joint(cx, gate, 7800.0)
            .with_bow_station_motion(limits(), gate)
            .unwrap();
        let zero = vec![0.0; direct.system().mode_count()];
        let mut drives = [drive(0.0, 0.0); 2];
        let mut cursor = 0;
        let mut expected = Vec::new();
        for sample in 0..513 {
            while events.get(cursor).is_some_and(|e| e.sample == sample) {
                drives[events[cursor].port] = events[cursor].drive;
                cursor += 1;
            }
            expected.push(
                direct
                    .step_with_bow_stations(&zero, &drives, &path(sample + 1), None)
                    .unwrap()
                    .observer_pressure_pa
                    .to_bits(),
            );
        }
        for block in [1, 37, 128] {
            let mut performance = create().with_bow_station_motion(program(), gate).unwrap();
            assert_eq!(performance.bow_station_program(), Some(&program()));
            let mut pressure = vec![0.0; 513];
            for chunk in pressure.chunks_mut(block) {
                performance.block(chunk).unwrap();
            }
            assert_eq!(
                pressure.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
                expected
            );
            assert_eq!(
                performance.bank().bow_station_fractions(),
                Some([0.205, 0.31].as_slice())
            );
            assert_eq!(states(performance.bank()), states(&direct));
        }
        let mut fixed = create();
        let mut fixed_pressure = vec![0.0; 513];
        for chunk in fixed_pressure.chunks_mut(128) {
            fixed.block(chunk).unwrap();
        }
        assert!(
            fixed_pressure
                .iter()
                .map(|p| p.to_bits())
                .zip(expected)
                .any(|(a, b)| a != b)
        );
    });
}

#[test]
fn bounded_program_admits_all_endpoints_and_refuses_unrepresentable_or_incomplete_motion() {
    let gate = CancelGate::new_clock_free();
    let initial = [0.2, 0.3];
    let lengths = [0.5, 0.5];
    program()
        .validate(&initial, &lengths, 513, 48000, &gate)
        .unwrap();
    for case in 0..10 {
        let mut bad = program();
        match case {
            0 => bad.knots[0].sample = 0,
            1 => bad.knots[1].sample = 128,
            2 => bad.knots[1].sample = 514,
            3 => bad.knots[0].station_fractions.clear(),
            4 => bad.knots[0].station_fractions[0] = f64::NAN,
            5 => bad.knots[0].station_fractions[0] = 1.0,
            6 => bad.max_knots = 1,
            7 => bad.max_coordinate_values = 3,
            8 => bad.limits.maximum_speed_m_s = 1e-5,
            _ => bad.limits.maximum_step_m = 1e-8,
        }
        assert!(
            bad.validate(&initial, &lengths, 513, 48000, &gate).is_err(),
            "case {case}"
        );
    }
    with_cx(|cx, gate| {
        assert!(
            joint(cx, gate, 7800.0)
                .with_bow_station_motion(
                    BowStationMotionLimits {
                        max_update_terms: 1,
                        ..limits()
                    },
                    gate
                )
                .is_err()
        );
        let bank = joint(cx, gate, 7800.0)
            .with_bow_station_motion(limits(), gate)
            .unwrap();
        assert!(bank.with_bow_station_motion(limits(), gate).is_err());
    });
}

#[test]
fn finite_ribbon_moves_as_one_interval_and_keeps_width_quadrature() {
    with_cx(|cx, gate| {
        let mut bank = base(cx, gate, 7800.0)
            .with_uniform_bow_ribbon(0.01, 3, ribbon_limits(), gate)
            .unwrap()
            .with_bow_station_motion(limits(), gate)
            .unwrap();
        let ribbon = bank.source().ribbon().unwrap().clone();
        let zero = vec![0.0; bank.system().mode_count()];
        let mut drives = vec![drive(0.0, 0.0); 3];
        ribbon.project_drive(drive(0.01, 0.1), &mut drives).unwrap();
        for sample in 0..16 {
            let end = UniformBowRibbon::try_new(
                StringBowLocation {
                    string: 0,
                    station_fraction: 0.2 + 0.00001 * (sample + 1) as f64,
                },
                bank.cards()[0].length_m,
                0.01,
                3,
            )
            .unwrap();
            bank.step_with_bow_stations(&zero, &drives, end.stations(), None)
                .unwrap();
            assert_eq!(bank.bow_station_fractions(), Some(end.stations()));
            assert_eq!(bank.source().ribbon(), Some(&ribbon));
        }
        let before = states(&bank);
        let stations = bank.bow_station_fractions().unwrap().to_vec();
        let mut deformed = stations.clone();
        deformed[1] += 0.0001;
        assert!(
            bank.step_with_bow_stations(&zero, &drives, &deformed, None)
                .is_err()
        );
        assert_eq!(states(&bank), before);
        assert_eq!(bank.bow_station_fractions(), Some(stations.as_slice()));
        let bad = BowStationProgram {
            knots: vec![BowStationKnot {
                sample: 1,
                station_fractions: deformed,
            }],
            limits: limits(),
            max_knots: 1,
            max_coordinate_values: 3,
        };
        let source = base(cx, gate, 7800.0)
            .with_uniform_bow_ribbon(0.01, 3, ribbon_limits(), gate)
            .unwrap();
        let gesture = GestureSchedule::try_new(
            1000,
            vec![GestureTrack {
                id: "ribbon".into(),
                target: GestureTarget::BowStroke { string: 0 },
                initial: value(0.01, 0.1, 0.2),
                events: vec![],
            }],
        )
        .unwrap();
        assert!(
            source
                .into_scheduled_ribbon(gesture, "ribbon".into(), performance_config())
                .unwrap()
                .with_bow_station_motion(bad, gate)
                .is_err()
        );
    });
}

fn narrow_bank(cx: &Cx<'_>, gate: &CancelGate, rate: u32, fluid: bool) -> MultiBowStringBank {
    let body = ModalAcousticTimeModel::try_new(
        rate,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 900.0,
            damping_ratio: 0.01,
            pressure_per_modal_velocity: C64::new(1.0, 0.0),
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let mut options = options();
    options.sample_rate_hz = rate;
    options.network.max_connections = 32;
    let member = MaterialStringMember {
        source: source(7800.0),
        connections: vec![InteriorStringConnection {
            station_fraction: 0.8,
            body: ModalAttachment {
                component: 0,
                shapes: vec![1.0],
            },
            stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1,
            rest_extension_m: 0.0,
        }],
    };
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let mut bank =
        ReciprocalStringBank::new(cx, vec![member], gas, vec![body], law(), options, gate).unwrap();
    if fluid {
        bank = bank
            .with_cylinder_gas_loading(
                fs_couple::air_path::cylinder::CylinderFlowResolution {
                    minimum_wavenumber_radius: 0.01,
                    maximum_wavenumber_radius: 10.0,
                    intervals: 2,
                },
                gate,
            )
            .unwrap();
    }
    bank.with_additional_bows(vec![], joint_limits(), gate)
        .unwrap()
}

#[test]
fn moving_string_port_preserves_causal_gas_coordinates_and_source_binding() {
    with_cx(|cx, gate| {
        let mut actual = narrow_bank(cx, gate, 48000, true)
            .with_bow_station_motion(limits(), gate)
            .unwrap();
        let (_, mut reference) = narrow_bank(cx, gate, 48000, true).into_parts();
        let card = actual.cards()[0].clone();
        let zero = vec![0.0; actual.system().mode_count()];
        let drives = [drive(0.3, 0.0001)];
        let mut start = 0.2;
        for sample in 0..64 {
            let end = 0.2 + 0.00001 * (sample + 1) as f64;
            let geometry = [FrictionGeometry {
                left: ModalAttachment {
                    component: 0,
                    shapes: sine(&card, start + 0.5 * (end - start)),
                },
                right: None,
            }];
            let a = actual
                .step_with_bow_stations(&zero, &drives, &[end], None)
                .unwrap()
                .clone();
            let b = reference
                .step_with_friction_geometry(
                    &zero,
                    &drives,
                    &geometry,
                    None,
                    None,
                    FrictionGeometryConfig {
                        max_update_terms: limits().max_update_terms,
                    },
                    gate,
                )
                .unwrap()
                .clone();
            assert_eq!(a, b);
            start = end;
        }
        assert_eq!(
            actual
                .source()
                .bank()
                .cylinder_gas_loading()
                .unwrap()
                .strings()
                .len(),
            1
        );
        assert!(
            actual.system().components()[2..]
                .iter()
                .flat_map(|m| m.states())
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0)
        );
        assert_eq!(
            states(&actual),
            reference
                .components()
                .iter()
                .flat_map(|m| m.states().iter().copied())
                .collect::<Vec<_>>()
        );
    });
}

#[test]
fn smooth_migration_has_a_second_order_time_refinement_path() {
    with_cx(|cx, gate| {
        let mut final_states = Vec::new();
        for rate in [24000, 48000, 96000, 192000] {
            let mut bank = narrow_bank(cx, gate, rate, false)
                .with_bow_station_motion(limits(), gate)
                .unwrap();
            let zero = vec![0.0; bank.system().mode_count()];
            let drives = [drive(0.5, 0.00002)];
            let steps = rate / 200;
            for sample in 0..steps {
                let position = 0.2 + 0.1 * f64::from(sample + 1) / f64::from(steps);
                let frame = bank
                    .step_with_bow_stations(&zero, &drives, &[position], None)
                    .unwrap();
                assert!(frame.ports[0].slip_distance_m > 0.0);
                near(frame.ports[0].traction_n, 0.4 * drives[0].normal_force_n);
            }
            final_states.push(states(&bank));
        }
        let difference = |left: &[ModalAcousticState], right: &[ModalAcousticState]| {
            left.iter()
                .zip(right)
                .map(|(a, b)| {
                    (a.displacement_m_sqrt_kg - b.displacement_m_sqrt_kg).powi(2)
                        + 1e-6 * (a.velocity_m_sqrt_kg_per_s - b.velocity_m_sqrt_kg_per_s).powi(2)
                })
                .sum::<f64>()
                .sqrt()
        };
        let errors: Vec<_> = final_states
            .windows(2)
            .map(|pair| difference(&pair[0], &pair[1]))
            .collect();
        assert!(errors.iter().all(|e| *e > 0.0));
        for pair in errors.windows(2) {
            let order = (pair[0] / pair[1]).log2();
            assert!(
                (order - 2.0).abs() < 0.15,
                "observed midpoint order {order}; errors {errors:?}"
            );
        }
    });
}

#[test]
fn migrating_bows_and_independent_stops_share_work_clock_and_callback_state() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
        performance::{
            motion::{ObstacleMotionLimits, StopMotionKnot, StopMotionProgram},
            stations::{StopStationKnot, StopStationProgram},
        },
        stops::{MaterialStringStop, StringStopSide, stations::StopStationMotionLimits},
    };
    use fs_couple::render::schedule::force::coupled::contact::{
        ModalContactConfig, friction::multiple::obstacles::FrictionObstacleConfig,
    };
    use fs_dcontact::Obstacle;
    with_cx(|cx, gate| {
        let make = || {
            joint(cx, gate, 7800.0)
                .with_string_stops(
                    vec![MaterialStringStop {
                        location: StringBowLocation {
                            string: 0,
                            station_fraction: 0.43,
                        },
                        side: StringStopSide::Positive,
                        body: None,
                        law: Obstacle::new(
                            vec![-1.0],
                            1,
                            1,
                            vec![0.0],
                            vec![1.0],
                            5000.0,
                            1.0,
                            "synthetic independent moving finger".into(),
                        )
                        .unwrap()
                        .with_internal_loss(0.2)
                        .unwrap(),
                        limits: ModalContactConfig {
                            max_iterations: 96,
                            maximum_force_n: 10.0,
                            maximum_penetration_m: 0.001,
                            force_absolute_tolerance_n: 1e-11,
                            force_relative_tolerance: 1e-9,
                        },
                    }],
                    FrictionObstacleConfig {
                        max_contacts: 1,
                        max_setup_terms: 1_048_576,
                    },
                    gate,
                )
                .unwrap()
        };
        let stops = StopStationProgram {
            knots: vec![StopStationKnot {
                sample: 256,
                station_fractions: vec![0.434],
            }],
            limits: StopStationMotionLimits {
                maximum_speed_m_s: 2.0,
                maximum_step_m: 0.0001,
                max_update_terms: 1_048_576,
            },
            max_knots: 1,
            max_coordinate_values: 1,
        };
        let press = StopMotionProgram {
            knots: vec![
                StopMotionKnot {
                    sample: 128,
                    offsets_m: vec![-1e-5],
                },
                StopMotionKnot {
                    sample: 320,
                    offsets_m: vec![0.0],
                },
            ],
            limits: ObstacleMotionLimits {
                maximum_abs_offset_m: 0.001,
                maximum_speed_m_s: 0.1,
            },
            max_knots: 2,
            max_coordinate_values: 2,
        };
        let events = make()
            .into_scheduled(schedule(), tracks(), performance_config())
            .unwrap()
            .events()
            .to_vec();
        let mut direct = make()
            .with_stop_motion(press.limits, gate)
            .unwrap()
            .with_stop_station_motion(stops.limits, gate)
            .unwrap()
            .with_bow_station_motion(limits(), gate)
            .unwrap();
        let zero = vec![0.0; direct.system().mode_count()];
        let mut drives = [drive(0.0, 0.0); 2];
        let mut cursor = 0;
        let mut expected = Vec::new();
        let (mut old_station, mut old_offset) = (0.43, 0.0);
        let (mut work, mut loss, mut migration_work) = (0.0, 0.0, 0.0);
        for sample in 0..513 {
            while events.get(cursor).is_some_and(|e| e.sample == sample) {
                drives[events[cursor].port] = events[cursor].drive;
                cursor += 1;
            }
            let endpoint = sample + 1;
            let stop = if endpoint < 256 {
                0.43 + (0.434 - 0.43) * (endpoint as f64 / 256.0)
            } else {
                0.434
            };
            let offset = if endpoint < 128 {
                -1e-5 * (endpoint as f64 / 128.0)
            } else if endpoint < 320 {
                -1e-5 + 1e-5 * ((endpoint - 128) as f64 / 192.0)
            } else {
                0.0
            };
            let before = direct.system().components()[0].states().to_vec();
            let frame = direct
                .step_with_bow_and_stop_stations_under_gate(
                    &zero,
                    &drives,
                    &path(endpoint),
                    Some(&[stop]),
                    Some(&[offset]),
                    gate,
                )
                .unwrap()
                .clone();
            expected.push(frame.observer_pressure_pa.to_bits());
            let contact = direct.system().last_obstacle_frames().unwrap()[0];
            let old_map = sine(&direct.cards()[0], old_station);
            let end_map = sine(&direct.cards()[0], stop);
            let migration = contact.normal_force_n
                * old_map
                    .iter()
                    .zip(end_map)
                    .zip(before.iter().zip(direct.system().components()[0].states()))
                    .map(|((&a, b), (q0, q1))| {
                        (b - a) * 0.5 * (q0.displacement_m_sqrt_kg + q1.displacement_m_sqrt_kg)
                    })
                    .sum::<f64>();
            let support = -contact.normal_force_n * (offset - old_offset);
            near(contact.attachment_work_j, migration);
            near(frame.obstacle_work_j, migration + support);
            work += frame.external_work_j + frame.drive_work_j + migration + support;
            loss += frame.network_dissipation_j
                + frame.friction_dissipation_j
                + frame.obstacle_dissipation_j;
            assert!((direct.system().total_energy_j().unwrap() + loss - work).abs() < 1e-10);
            migration_work += migration.abs();
            old_station = stop;
            old_offset = offset;
        }
        assert!(migration_work > 0.0);
        for block in [1, 37, 128] {
            let mut performance = make()
                .into_scheduled(schedule(), tracks(), performance_config())
                .unwrap()
                .with_stop_motion(press.clone(), gate)
                .unwrap()
                .with_stop_station_motion(stops.clone(), gate)
                .unwrap()
                .with_bow_station_motion(program(), gate)
                .unwrap();
            let mut pressure = vec![0.0; 513];
            for chunk in pressure.chunks_mut(block) {
                performance.block(chunk).unwrap();
            }
            assert_eq!(
                pressure.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
                expected
            );
            assert_eq!(states(performance.bank()), states(&direct));
            assert_eq!(
                performance.bank().stop_station_fractions(),
                Some([0.434].as_slice())
            );
            assert_eq!(
                performance.bank().bow_station_fractions(),
                Some([0.205, 0.31].as_slice())
            );
        }
        // A refusal must publish neither endpoint, even when stop projection
        // succeeds before the new bow endpoint reaches the mechanical owner.
        let before = states(&direct);
        let frame = direct.system().last_frame().unwrap().clone();
        assert!(
            direct
                .step_with_bow_and_stop_stations_under_gate(
                    &[],
                    &drives,
                    &[0.20501, 0.31001],
                    Some(&[0.43401]),
                    Some(&[0.0]),
                    gate
                )
                .is_err()
        );
        assert_eq!(states(&direct), before);
        assert_eq!(direct.system().last_frame(), Some(&frame));
        assert_eq!(direct.stop_station_fractions(), Some([0.434].as_slice()));
        assert_eq!(
            direct.bow_station_fractions(),
            Some([0.205, 0.31].as_slice())
        );
    });
}
