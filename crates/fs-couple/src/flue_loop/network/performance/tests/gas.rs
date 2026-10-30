//! Synthetic constitutive integration controls; no measured flute claim.
use super::*;
use crate::air_path::GasAbsorptionModel;
use crate::bernoulli_aperture::fractional::FractionalTubePolicy;
use crate::flue_loop::network::performance::file::FlueFilePerformance;
use crate::pcm_wav::baffled::gas::GasPropagation;
use crate::pcm_wav::baffled::openings::CircularOpening;
use fs_exec::CancelGate;
use fs_material::gas_species::input::GasRecord;

const MOLECULAR: &str = "molecular 4096 1000000 50000000 declared 2 N2 1 vibration 1 1000 1 power 10 300 0 250 350 O2 0";

fn source_text(law: Option<&str>) -> String {
    let mut text = String::from(concat!(
        "frankensim-flue-performance-v1\nclock 96000 48000\noutput 768 5 1000000\n",
        "gas 293.15 101325 0\nexterior mixture 300 101325 mole 2 N2 0.8 O2 0.2\n",
        "network 4 3\nduct_node inlet\nduct_node junction\n",
        "duct_node reflection -0.9\nduct_node reflection -0.8\n",
        "duct_section 0 1 0.15 0.006 1e-10\n",
        "duct_section 1 2 0.12 0.006 1e-10\n",
        "duct_section 1 3 0.18 0.004 1e-10\n",
        "flue 0.010 0.001 0.004 0 0 0.02\nwindway none\njet_card tonal-interim\n",
        "jet 90 1e-6 3 0\njet_limits 1 0.3 8192\nwave_memory_bytes 1048576\n",
        "fractional 1000 0.01\nviscothermal none\nreceiver 0.02 0.01 0.2 4 8 1000\n",
        "outlets 2 1\noutlet 3 0.05 0\nmouth_aperture rectangle 0.012 0.009\n",
        "mouth_radiation -0.15 0 4 2\n",
    ));
    if let Some(law) = law {
        let receiver = "receiver 0.02 0.01 0.2 4 8 1000\n";
        text = text.replace(receiver, &format!("{receiver}observer_gas {law}\n"));
    }
    text.push_str("controls 4096 100000\ngesture\n");
    text.push_str(std::str::from_utf8(&schedule(300.0).to_canonical_bytes()).unwrap());
    text
}

fn read(text: &str, block: usize) -> FlueFilePerformance {
    FlueFilePerformance::from_bytes(text.as_bytes(), block, &CancelGate::new_clock_free()).unwrap()
}

fn render(
    mut performance: FlueOutletPerformance,
    block: usize,
) -> (Vec<f64>, FlueOutletPerformance) {
    let mut pressure = vec![0.0; performance.config().samples as usize];
    for output in pressure.chunks_mut(block) {
        performance.block(output).unwrap();
    }
    (pressure, performance)
}

#[test]
fn molecular_gas_reaches_every_original_aperture_and_preserves_source_motion() {
    let gate = CancelGate::new_clock_free();
    let input = source_text(Some(MOLECULAR));
    let file = read(&input, 37);
    let exterior = file.gas_records()[1].clone();
    assert!(matches!(file.gas_records()[0], GasRecord::Legacy(_)));
    let renderer = file.renderer();
    let terminal_report = renderer.observer_gas().unwrap();
    let mouth_report = renderer.mouth_observer_gas().unwrap();
    assert_eq!(terminal_report.law, GasAbsorptionModel::MolecularRelaxation);
    assert_eq!(&terminal_report.gas, exterior.state());
    assert_ne!(&terminal_report.gas, renderer.voice().gas());
    assert_eq!(terminal_report.gas, mouth_report.gas);
    assert_eq!(terminal_report.surface_filter_indices.len(), 2 * 4 * 8);
    assert_eq!(mouth_report.surface_filter_indices.len(), 4 * 2);
    assert!(terminal_report.filters.len() > 1);
    assert_eq!(
        terminal_report.additional_latency_samples,
        mouth_report.additional_latency_samples
    );
    assert_eq!(
        renderer.receiver().filter_delay_samples(),
        renderer.mouth_receiver().unwrap().filter_delay_samples()
    );

    // Feed the actual source flows directly to the common physical receivers.
    // This checks all terminal signs, the independent mouth aperture, and every
    // 96 kHz mechanical sample entering the 48 kHz microphone.
    let mut direct = FlueVoice::from_fractional_network(
        Some(file.card()),
        90.0,
        1e-6,
        file.geometry(),
        renderer.voice().network().unwrap().source().clone(),
        None,
        FractionalTubePolicy {
            maximum_frequency_hz: 1000.0,
            maximum_summed_phase_error_rad: 0.01,
        },
        96_000,
        3,
        0,
        &gate,
    )
    .unwrap();
    direct
        .apply(FlueControl::SetBlowingPressure(300.0))
        .unwrap();
    let selection =
        GasPropagation::read_record_for_gas(&mut MOLECULAR.split_ascii_whitespace(), &exterior)
            .unwrap();
    let config = renderer.config().receiver;
    let (mut terminals, mut terminals_report) =
        BaffledPressure::circular_openings_multirate_gas_differentiating(
            &[
                CircularOpening {
                    center_m: [0.0, 0.0],
                    radius_m: 0.006,
                },
                CircularOpening {
                    center_m: [0.05, 0.0],
                    radius_m: 0.004,
                },
            ],
            96_000,
            48_000,
            config,
            exterior.state(),
            &selection,
            &mut || false,
        )
        .unwrap();
    let (mut mouth, mut mouth_report) =
        BaffledPressure::rectangular_outlet_multirate_gas_differentiating(
            [0.012, 0.009],
            96_000,
            48_000,
            RectangularOutletReceiver {
                center_m: [-0.15, 0.0],
                position_m: config.position_m,
                cells: [4, 2],
                maximum_frequency_hz: config.maximum_frequency_hz,
            },
            exterior.state(),
            &selection,
            &mut || false,
        )
        .unwrap();
    let latency = terminals_report
        .additional_latency_samples
        .max(mouth_report.additional_latency_samples);
    terminals
        .align_gas_latency(&mut terminals_report, latency)
        .unwrap();
    mouth.align_gas_latency(&mut mouth_report, latency).unwrap();
    let mut reference = vec![0.0; 768];
    for expected in &mut reference {
        let (mut outlet_trace, mut mouth_trace) = ([0.0; 4], [0.0; 2]);
        for frame in 0..2 {
            direct.step_block(&mut [0.0]).unwrap();
            let network = direct.network().unwrap();
            outlet_trace[2 * frame] = network.node_frame(2).unwrap().net_flow_into_node_m3_s;
            outlet_trace[2 * frame + 1] = network.node_frame(3).unwrap().net_flow_into_node_m3_s;
            mouth_trace[frame] = network.node_frame(0).unwrap().net_flow_into_node_m3_s;
        }
        *expected =
            terminals.step_trace(&outlet_trace).unwrap() + mouth.step_trace(&mouth_trace).unwrap();
    }
    let (actual, completed) = render(file.into_renderer(), 37);
    assert_eq!(actual, reference);
    assert!(actual.iter().any(|p| p.abs() > 1e-12));
    assert_eq!(
        completed.voice().network().unwrap().last_frame(),
        direct.network().unwrap().last_frame()
    );
    let (whole, _) = render(read(&input, 768).into_renderer(), 768);
    assert_eq!(
        actual, whole,
        "callback partition cannot reset path or multirate history"
    );

    let (unabsorbed, unchanged) = render(read(&source_text(None), 37).into_renderer(), 37);
    assert!(
        actual
            .iter()
            .zip(&unabsorbed)
            .any(|(a, b)| (a - b).abs() > 1e-12)
    );
    assert_eq!(
        completed.voice().network().unwrap().last_frame(),
        unchanged.voice().network().unwrap().last_frame()
    );
    assert_eq!(
        completed.voice().network().unwrap().stored_energy_j(),
        unchanged.voice().network().unwrap().stored_energy_j()
    );
    assert_eq!(
        completed.voice().mouth_flow_m3_s(),
        unchanged.voice().mouth_flow_m3_s()
    );

    let lossless = source_text(Some("lossless 4096 1000000 50000000"));
    let (explicit, _) = render(read(&lossless, 37).into_renderer(), 37);
    for (a, b) in explicit.iter().zip(&unabsorbed) {
        assert!(
            (a - b).abs() < 1e-11 * b.abs().max(1.0),
            "lossless control retains original receiver"
        );
    }
}

#[test]
fn published_flute_cards_use_their_declared_exterior_without_changing_calibration() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/instruments/descriptions/flute");
    for name in [
        "flute-d4.performance",
        "flute-g4.performance",
        "flute-a5.performance",
    ] {
        let bytes = std::fs::read(root.join(name)).expect("published flute source");
        let source = FlueFilePerformance::from_bytes(&bytes, 97, &CancelGate::new_clock_free())
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(
            source.info().full_scale_pa,
            5.0,
            "{name}: physical PCM calibration"
        );
        assert_eq!(source.info().samples, 120_000);
        assert_eq!(source.info().sample_rate_hz, 48_000);
        assert_eq!(source.info().mechanics_rate_hz, 288_000);
        let renderer = source.renderer();
        let exterior = renderer.radiation_input().unwrap().exterior;
        let terminal = renderer
            .observer_gas()
            .expect("published terminal propagation");
        let mouth = renderer
            .mouth_observer_gas()
            .expect("published embouchure propagation");
        assert_eq!(terminal.law, GasAbsorptionModel::StokesKirchhoff);
        assert_eq!(mouth.law, GasAbsorptionModel::StokesKirchhoff);
        assert_eq!(terminal.gas, exterior);
        assert_eq!(mouth.gas, exterior);
        assert_eq!(exterior, *source.gas_records()[0].state());
        assert_eq!(
            terminal.additional_latency_samples,
            mouth.additional_latency_samples
        );
        assert_eq!(
            renderer.receiver().filter_delay_samples(),
            renderer.mouth_receiver().unwrap().filter_delay_samples()
        );
        assert_eq!(renderer.voice().island().displacement_gain, 90.0);
        assert_eq!(renderer.voice().island().seed_amplitude_rel, 1e-6);
        assert_eq!(source.card().authority.label(), "X-Est");

        // Retain the complete original gesture and calibration; render only an
        // onset prefix so the admission test does not rerender each full phrase.
        let mut renderer = source.into_renderer();
        let mut pressure = vec![0.0; 1536];
        for block in pressure.chunks_mut(97) {
            renderer.block(block).unwrap();
        }
        assert!(pressure.iter().all(|p| p.is_finite()));
        assert!(
            pressure.iter().any(|p| p.abs() > 1e-12),
            "{name}: physical onset reaches the microphone"
        );
        assert_eq!(renderer.voice().samples_rendered(), 1536 * 6);
    }
}
