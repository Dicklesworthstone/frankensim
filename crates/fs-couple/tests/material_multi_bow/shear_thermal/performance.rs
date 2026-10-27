//! The canonical schedule and PCM route must use the same finite heat owner.
use super::*;
use fs_couple::pcm_wav::{
    encode_pcm16_wav,
    stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16},
};
use std::io::Cursor;

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
        let make = || {
            instrument(cx, gate, &card)
                .into_scheduled(held_schedule(), tracks(), performance_config())
                .unwrap()
        };
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
        let mut a = instrument(cx, gate, &card)
            .into_scheduled(held_schedule(), tracks(), performance_config())
            .unwrap();
        let mut b = instrument(cx, gate, &constant)
            .into_scheduled(held_schedule(), tracks(), performance_config())
            .unwrap();
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
        let mut performance = instrument(cx, gate, &card)
            .into_scheduled(held_schedule(), tracks(), performance_config())
            .unwrap();
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
        assert!(renderer.last_frame().is_none());
        assert!(matches!(
            renderer.block(&mut output),
            Err(RenderError::Poisoned)
        ));
    });
}
