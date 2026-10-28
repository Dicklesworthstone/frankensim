//! The canonical schedule and PCM route must use the same finite heat owner.
use super::*;
use fs_couple::pcm_wav::{
    encode_pcm16_wav,
    stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16},
};
use std::io::Cursor;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    performance::{MultiBowPerformance, motion::{StopMotionKnot, StopMotionProgram, ObstacleMotionLimits},
        stations::{StopStationKnot, StopStationProgram, StopStationMotionLimits}},
    stops::{MaterialStringStop, StringStopSide},
};
use fs_couple::render::schedule::force::coupled::contact::{
    ModalContactConfig, friction::multiple::obstacles::FrictionObstacleConfig,
};
use fs_dcontact::Obstacle;

fn finger_bank(cx: &Cx<'_>, gate: &CancelGate, card: &MaterialCard) -> MultiBowStringBank {
    instrument(cx, gate, card).with_string_stops(vec![MaterialStringStop {
        location: StringBowLocation { string: 0, station_fraction: 0.43 },
        side: StringStopSide::Positive, body: None,
        law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 5000.0, 1.0,
            "manufactured compliant finger, not measured tissue".into())
            .unwrap().with_internal_loss(0.2).unwrap(),
        limits: ModalContactConfig { max_iterations: 96, maximum_force_n: 10.0,
            maximum_penetration_m: 0.001, force_absolute_tolerance_n: 1e-11,
            force_relative_tolerance: 1e-9 },
    }], FrictionObstacleConfig { max_contacts: 1, max_setup_terms: 1_048_576 }, gate).unwrap()
}
fn press() -> StopMotionProgram {
    StopMotionProgram {
        knots: [(12, -3e-6), (36, -3e-6), (72, -1e-6)].into_iter()
            .map(|(sample, offset)| StopMotionKnot { sample, offsets_m: vec![offset] }).collect(),
        limits: ObstacleMotionLimits { maximum_abs_offset_m: 0.001, maximum_speed_m_s: 0.1 },
        max_knots: 3, max_coordinate_values: 3,
    }
}
fn migrate() -> StopStationProgram {
    StopStationProgram {
        knots: vec![StopStationKnot { sample: 128, station_fractions: vec![0.4312] }],
        limits: StopStationMotionLimits { maximum_speed_m_s: 1.0, maximum_step_m: 0.00003,
            max_update_terms: 1_048_576 },
        max_knots: 1, max_coordinate_values: 1,
    }
}
fn finger_performance(cx: &Cx<'_>, gate: &CancelGate, card: &MaterialCard) -> MultiBowPerformance {
    finger_bank(cx, gate, card).into_scheduled(held_schedule(), tracks(), performance_config())
        .unwrap().with_stop_motion(press(), gate).unwrap()
        .with_stop_station_motion(migrate(), gate).unwrap()
}

fn held_schedule() -> GestureSchedule {
    GestureSchedule::try_new(
        700,
        vec![
            GestureTrack {
                id: "bow-0".into(),
                target: GestureTarget::BowStroke { string: 0 },
                initial: value(0.001, 0.1, 0.2),
                events: vec![
                    GestureEvent {
                        time_s: 0.002,
                        transition_s: 0.001,
                        value: value(-0.001, 0.1, 0.2),
                    },
                    GestureEvent {
                        time_s: 0.006,
                        transition_s: 0.001,
                        value: value(0.001, 0.1, 0.2),
                    },
                ],
            },
            GestureTrack {
                id: "bow-1".into(),
                target: GestureTarget::BowStroke { string: 1 },
                initial: value(0.0, 0.0, 0.3),
                events: vec![],
            },
        ],
    )
    .unwrap()
}

#[test]
fn scheduled_heat_preserves_canonical_controls_callback_history_and_original_pcm() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, false);
        let make = || finger_performance(cx, gate, &card);
        let mut one = make();
        let mut chunks = make();
        let mut reference = make();
        let mut pcm = make();
        let mut a = one
            .with_contact_thermal_network(
                thermal(caloric.curve(), true),
                binding(&card, caloric.curve()).to_vec(),
                vec![],
                controls(3),
                gate,
                other,
            )
            .unwrap();
        let mut b = chunks
            .with_contact_thermal_network(
                thermal(caloric.curve(), true),
                binding(&card, caloric.curve()).to_vec(),
                vec![],
                controls(3),
                gate,
                other,
            )
            .unwrap();
        let mut expected = vec![0.0; 513];
        let mut actual = expected.clone();
        let mut ordinary = expected.clone();
        for sample in expected.chunks_mut(1) {
            a.block(sample).unwrap();
        }
        for block in actual.chunks_mut(37) {
            b.block(block).unwrap();
        }
        for block in ordinary.chunks_mut(128) {
            reference.block(block).unwrap();
        }
        // Constant source coefficients: thermal/Jacobian trials may not age the
        // mechanical history or change the original scheduled pressure at all.
        assert_eq!(expected, actual);
        assert_eq!(expected, ordinary);
        assert_eq!(states(a.performance().bank()), states(reference.bank()));
        assert_eq!(a.thermal().cells(), b.thermal().cells());
        assert_eq!(a.thermal().time_s(), b.thermal().time_s());
        assert_eq!(a.performance().bank().stop_station_fractions(), Some([0.4312].as_slice()));
        assert_eq!(a.performance().bank().system().obstacle_offsets_m(), Some([-1e-6].as_slice()));
        assert_eq!(a.performance().bank().source().stops()[0].location.station_fraction, 0.43);
        assert_eq!(a.last_frame().unwrap().mechanics.obstacle_work_j, 0.0);
        assert_eq!(
            a.performance().bank().system().material_shear_state(0),
            b.performance().bank().system().material_shear_state(0)
        );
        assert!(total_heat(a.thermal()) > 0.0);
        assert!(a.thermal().temperature_k(CellId(11)).unwrap() > 300.0);
        assert_eq!(
            a.performance()
                .bank()
                .bow_shear_layer(0)
                .unwrap()
                .temperature_kelvin(),
            a.thermal().temperature_k(FILM).unwrap()
        );
        assert!(
            a.last_frame().unwrap().total_energy_residual_j.abs()
                <= controls(3).total_energy_tolerance_j
        );
        let recovered = a.into_thermal();
        assert_eq!(recovered.time_s(), 513.0 * (1.0 / 48_000.0));
        assert_eq!(one.samples_rendered(), 513);

        let (bytes, clipped) = encode_pcm16_wav(&expected, 48_000, 0.01).unwrap();
        let mut renderer = pcm
            .with_contact_thermal_network(
                thermal(caloric.curve(), true),
                binding(&card, caloric.curve()).to_vec(),
                vec![],
                controls(3),
                gate,
                other,
            )
            .unwrap();
        let mut writer = Pcm16WavStream::new(Cursor::new(Vec::new()), 48_000, 0.01, 128).unwrap();
        let progress =
            render_pressure_pcm16(&mut renderer, &mut writer, gate, &mut [0.0; 37], 513).unwrap();
        assert!(matches!(progress, ScheduledWavProgress::Completed { .. }));
        let (output, report) = writer.finish().unwrap();
        assert_eq!(output.into_inner(), bytes);
        assert_eq!(report.samples, 513);
        assert_eq!(report.clipped_samples, clipped as u64);
        assert_eq!(renderer.thermal().cells(), b.thermal().cells());
    });
}

#[test]
fn temperature_dependent_source_changes_scheduled_pressure_and_keeps_claim_identity() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(true, false);
        let (constant, constant_caloric) = coating(false, false);
        let mut a = finger_performance(cx, gate, &card);
        let mut b = finger_performance(cx, gate, &constant);
        let initial = a.bank().bow_shear_layer(0).unwrap().clone();
        let mut warm = a
            .with_contact_thermal_network(
                thermal(caloric.curve(), false),
                binding(&card, caloric.curve()).to_vec(),
                vec![],
                controls(1),
                gate,
                other,
            )
            .unwrap();
        let mut fixed = b
            .with_contact_thermal_network(
                thermal(constant_caloric.curve(), false),
                binding(&constant, constant_caloric.curve()).to_vec(),
                vec![],
                controls(1),
                gate,
                other,
            )
            .unwrap();
        let mut actual = vec![0.0; 513];
        let mut baseline = actual.clone();
        for block in actual.chunks_mut(37) {
            warm.block(block).unwrap();
        }
        for block in baseline.chunks_mut(128) {
            fixed.block(block).unwrap();
        }
        assert!(
            actual
                .iter()
                .zip(&baseline)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>()
                > 1e-14
        );
        let current = warm.performance().bank().bow_shear_layer(0).unwrap();
        assert!(current.temperature_kelvin() > initial.temperature_kelvin());
        assert!(current.plastic_resistance_n_s_m() < initial.plastic_resistance_n_s_m());
        for (now, old) in current
            .resolved()
            .properties()
            .iter()
            .zip(initial.resolved().properties())
        {
            assert_eq!(now.answer().receipt.selected, old.answer().receipt.selected);
        }
        assert_eq!(
            current.resolved().card_identity(),
            initial.resolved().card_identity()
        );
    });
}

#[test]
fn failed_thermal_source_keeps_accepted_clocks_history_and_output_then_poisons_callback() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(true, false);
        let mut performance = finger_performance(cx, gate, &card);
        let before = states(performance.bank());
        let mut renderer = performance
            .with_contact_thermal_network(
                thermal(caloric.curve(), false),
                binding(&card, caloric.curve()).to_vec(),
                vec![],
                controls(1),
                gate,
                |_| Err::<StepCoefficients<()>, _>("original source refused"),
            )
            .unwrap();
        let cells = renderer.thermal().cells().clone();
        let mut output = [123.0];
        assert!(renderer.block(&mut output).is_err());
        assert_eq!(output, [123.0]);
        assert_eq!(renderer.samples_rendered(), 0);
        assert_eq!(renderer.thermal().time_s(), 0.0);
        assert_eq!(renderer.thermal().cells(), &cells);
        assert_eq!(states(renderer.performance().bank()), before);
        assert_eq!(renderer.performance().bank().stop_station_fractions(), Some([0.43].as_slice()));
        assert_eq!(renderer.performance().bank().system().obstacle_offsets_m(), Some([0.0].as_slice()));
        assert!(renderer.last_frame().is_none());
        assert!(matches!(
            renderer.block(&mut output),
            Err(RenderError::Poisoned)
        ));
    });
}

#[test]
fn moving_finger_heat_keeps_product_rule_work_and_rejected_geometry_can_retry() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, false);
        let make = || finger_bank(cx, gate, &card)
            .with_stop_motion(press().limits, gate).unwrap()
            .with_stop_station_motion(migrate().limits, gate).unwrap();
        let mut bank = make();
        let mut ordinary = make();
        let mut heat = thermal(caloric.curve(), false);
        let bound = binding(&card, caloric.curve());
        let zero = vec![0.0; bank.system().mode_count()];
        let drives = [drive(0.001, 0.1), drive(0.0, 0.0)];
        let initial_energy = bank.system().total_energy_j().unwrap();
        let (mut work, mut loss, mut migrated, mut deposited) = (0.0, 0.0, 0.0, 0.0);
        let (mut old_station, mut old_offset) = (0.43, 0.0);
        for endpoint in 1..=160 {
            let station = if endpoint < 128 { 0.43 + (0.4312 - 0.43) * (endpoint as f64 / 128.0) }
                else { 0.4312 };
            let offset = if endpoint < 12 { -3e-6 * (endpoint as f64 / 12.0) }
                else if endpoint <= 36 { -3e-6 }
                else if endpoint < 72 { -3e-6 + 2e-6 * ((endpoint - 36) as f64 / 36.0) }
                else { -1e-6 };
            let before = states(&bank);
            let guess = [total_heat(&heat)];
            if endpoint == 48 {
                let cells = heat.cells().clone();
                let accepted = bank.system().last_frame().cloned();
                let film = bank.system().material_shear_state(0).cloned();
                let column = bank.system().obstacle_normal_column(0).unwrap().to_vec();
                assert!(bank.step_with_contact_thermal_network_and_stops(
                    &mut heat, &bound, &[], controls(1), &guess, &zero, &drives,
                    Some(&[station]), Some(&[offset]), gate,
                    |_| Err::<StepCoefficients<()>, _>("source refuses after geometry staging"),
                ).is_err());
                assert_eq!(states(&bank), before);
                assert_eq!(heat.cells(), &cells);
                assert_eq!(heat.time_s(), 47.0 * bank.system().sample_period_s());
                assert_eq!(bank.system().samples_rendered(), 47);
                assert_eq!(bank.system().last_frame(), accepted.as_ref());
                assert_eq!(bank.system().material_shear_state(0), film.as_ref());
                assert_eq!(bank.system().obstacle_normal_column(0).unwrap(), column);
                assert_eq!(bank.stop_station_fractions(), Some([old_station].as_slice()));
                assert_eq!(bank.system().obstacle_offsets_m(), Some([old_offset].as_slice()));
            }
            let frame = bank.step_with_contact_thermal_network_and_stops(
                &mut heat, &bound, &[], controls(1), &guess, &zero, &drives,
                Some(&[station]), Some(&[offset]), gate, other,
            ).unwrap();
            let reference = ordinary.step_with_stop_stations(
                &zero, &drives, &[station], Some(&[offset])).unwrap();
            assert_eq!(frame.mechanics.observer_pressure_pa, reference.observer_pressure_pa);
            assert_eq!(states(&bank), states(&ordinary));
            // Independent homogeneous-string mass-normalized sine basis, not
            // the contact owner's staged geometry/work report.
            let source = &bank.cards()[0];
            let normalizer = (0.5 * source.linear_density_kg_m * source.length_m).sqrt();
            let contact = bank.system().last_obstacle_frames().unwrap()[0];
            let migration = contact.normal_force_n * before.iter()
                .zip(bank.system().components()[0].states()).take(source.mode_count)
                .enumerate().map(|(n, (a, b))| {
                    let k = (n + 1) as f64 * core::f64::consts::PI;
                    ((k * station).sin() - (k * old_station).sin()) / normalizer
                        * 0.5 * (a.displacement_m_sqrt_kg + b.displacement_m_sqrt_kg)
                }).sum::<f64>();
            let actuator = migration - contact.normal_force_n * (offset - old_offset);
            assert!((contact.attachment_work_j - migration).abs() < 1e-15);
            assert!((frame.mechanics.obstacle_work_j - actuator).abs() < 1e-15);
            assert!(frame.total_energy_residual_j.abs() <= controls(1).total_energy_tolerance_j);
            migrated += migration.abs();
            deposited += frame.routed_plastic_heat_j;
            work += frame.mechanics.external_work_j + frame.mechanics.drive_work_j + actuator;
            loss += frame.mechanics.network_dissipation_j + frame.mechanics.friction_dissipation_j
                + frame.mechanics.obstacle_dissipation_j - frame.routed_plastic_heat_j;
            assert!((bank.system().total_energy_j().unwrap() + total_heat(&heat)
                - initial_energy + loss - work).abs() < 1e-10);
            old_station = station;
            old_offset = offset;
        }
        assert!(migrated > 0.0 && deposited > 0.0);
        assert!((total_heat(&heat) - deposited).abs() < 1e-12);
        assert_eq!(bank.system().last_frame().unwrap().obstacle_work_j, 0.0);
    });
}
