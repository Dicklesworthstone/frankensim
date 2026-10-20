//! Native source selection reaches the original shared exterior/jet state.
use super::*;
use fs_couple::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationArtifact;
use fs_couple::flue_loop::{FlueVoice, network::performance::FlueOutletPerformance};
use fs_couple::flue_loop::network::performance::file::FlueInputError;

const PERFORMANCE: &str = include_str!("../../../../examples/flue-mimo.performance");
const ARTIFACT: &str = include_str!("../../../../examples/flue-mutual.mimo-radiation");
const SELECTION: &str = "radiation mimo flue-mutual.mimo-radiation 32 1048576 257";
const COMPUTED: &str = "radiation computed-baffled 20 65 4 8 12 20 0.01 67108864 2 100 1e-8 64 1048576 257";

fn computed_source() -> String {
    native_source(false).replace("clock 96000 48000", "clock 8000 8000")
        .replace("output 1024 1", "output 256 1")
        .replace("duct_section 0 1 0.18 0.006", "duct_section 0 1 0.18 0.002")
        .replace("duct_section 0 3 0.04 0.004", "duct_section 0 3 0.08 0.0015")
        .replace("outlet 3 0.08 0", "outlet 3 0.015 0")
        .replace("windway unsteady 0.015 32 64", "windway none")
        .replace("viscothermal sections 2\nloss_section 1 4 300 500\nloss_section 2 4 300 500", "viscothermal none")
        .replace("fractional 1000 0.01", "fractional 500 1")
        .replace("receiver 0 0 0.5 4 16 1000", "receiver 0 0 0.5 4 16 500")
        .replace(SELECTION, COMPUTED)
}

fn load_computed(text: &str, block: usize) -> FlueFilePerformance {
    FlueFilePerformance::from_bytes_with_all_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None,
        |_| panic!("this physical fixture has no solid walls"),
        |_| panic!("computed radiation must not require a matrix artifact")).unwrap()
}

fn computed_pressure(file: FlueFilePerformance) -> Vec<f64> {
    let samples = file.info().samples as usize;
    let mut renderer = file.into_renderer();
    let mut result = vec![0.0; samples];
    for chunk in result.chunks_mut(renderer.config().max_block) { renderer.block(chunk).unwrap(); }
    assert_eq!(renderer.voice().samples_rendered(), samples as u64);
    assert!(result.iter().all(|p| p.is_finite()));
    assert!(result.iter().any(|p| p.abs() > 1e-15));
    result
}

#[test]
fn computed_geometry_and_atmosphere_reach_the_original_jet_and_microphone_without_an_artifact() {
    let text = computed_source();
    let file = load_computed(&text, 37);
    assert!(file.mimo_radiation_record().is_none());
    assert!(file.radiation_bindings().is_empty());
    let record = file.computed_radiation_record().unwrap();
    assert_eq!(record.band_hz(), [20.0, 4000.0]);
    assert_eq!(record.apertures()[0].radius_m, 0.002);
    assert_eq!(record.apertures()[1].radius_m, 0.0015);
    assert_eq!(record.apertures()[1].center_m, [0.015, 0.0]);
    let environment = file.renderer().voice().network().unwrap().mimo_environment().unwrap();
    assert_eq!(record.gas_record().state(), &environment.source().exterior);
    assert_eq!(record.apertures(), environment.source().apertures);
    let identity = record.realization_identity();
    let matrix = environment.source().impedance.eval_iw(core::f64::consts::TAU * 713.0);
    let reference = computed_pressure(file);
    assert_eq!(computed_pressure(load_computed(&text, 1)), reference);

    // No new matrix file or coefficient edits: source gas data alone causes a
    // new physical source sweep and changes both the load and rendered sound.
    let changed = load_computed(&text.replace("101325", "81060"), 37);
    let record = changed.computed_radiation_record().unwrap();
    assert_eq!(record.gas_record().state().pressure, 81060.0);
    assert_ne!(record.realization_identity(), identity);
    let load = &changed.renderer().voice().network().unwrap().mimo_environment().unwrap().source().impedance;
    let updated = load.eval_iw(core::f64::consts::TAU * 713.0);
    for (a, b) in matrix.iter().zip(updated) {
        assert!((b - a.scale(0.8)).abs() < 0.001 * a.abs());
    }
    assert_ne!(computed_pressure(changed), reference);
}

#[test]
fn computed_physical_source_reaches_cli_pcm_and_reports_actual_fit_scope() {
    let text = include_str!("../../../../examples/flue-computed-radiation.performance");
    let pressure = computed_pressure(load_computed(text, 37));
    let expected = encode_pcm16_wav(&pressure, 48_000, 5.0).unwrap().0;
    let dir = directory();
    let input = dir.join("computed.performance");
    let output = dir.join("computed.wav");
    std::fs::write(&input, text).unwrap();
    success(run(&input, &output, 137, false));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    for token in ["\"source\":\"computed-baffled\"", "\"band_hz\":[1.2e2,2.4e4]",
        "\"held_out_samples\":64", "\"sampled_passivity_class\":\"GridOnly\"",
        "\"radiation_feedback_added\":true", "\"rungs\":[", "no passivity repair"] {
        assert!(metadata.contains(token), "{token}: {metadata}");
    }
}

fn native_source(walls: bool) -> String {
    let mut head = PERFORMANCE.split_once("\ngesture\n").unwrap().0
        .replace("output 24000 5", "output 1024 1");
    if walls {
        head = head.replace(SELECTION, &format!(
            "walls 1 20 1000\nwall 4 isotropic wall.fsmcdpk 293.15 1000 0.0002 stress-free kelvin-voigt\n{SELECTION}"));
    }
    format!("{head}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t3e2\nevents\t0\n")
}
fn load_mimo(text: &str, artifact: &str, material: &NormalizedMaterialCardPack, block: usize)
    -> Result<FlueFilePerformance, FlueInputError>
{
    FlueFilePerformance::from_bytes_with_all_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None, |name| {
            assert_eq!(name, "wall.fsmcdpk");
            Ok(material.clone())
        }, |name| {
            assert_eq!(name, "flue-mutual.mimo-radiation");
            Ok(artifact.as_bytes().to_vec())
        })
}
fn pressure(mut renderer: FlueOutletPerformance, block: usize) -> Vec<f64> {
    let mut result = vec![0.0; renderer.config().samples as usize];
    for chunk in result.chunks_mut(block) { renderer.block(chunk).unwrap(); }
    assert_eq!(renderer.voice().samples_rendered(), 2 * result.len() as u64);
    result
}
struct MutualTrace {
    pressure: Vec<f64>,
    incoming: Vec<f64>,
    work: f64,
    wall_energy: f64,
}
fn mutual_trace(file: FlueFilePerformance) -> MutualTrace {
    let mut renderer = file.into_renderer();
    let mut trace = MutualTrace { pressure: Vec::new(), incoming: Vec::new(), work: 0.0, wall_energy: 0.0 };
    while renderer.remaining_samples() != 0 {
        let mut one = [0.0];
        renderer.block(&mut one).unwrap();
        trace.pressure.push(one[0]);
        let network = renderer.voice().network().unwrap();
        let frame = network.last_frame().unwrap();
        trace.incoming.push(frame.incoming_pressure_pa);
        trace.work += frame.coupled_terminal_work_j.abs();
        let scale = (frame.stored_energy_j + frame.storage_change_j.abs()
            + frame.inlet_work_j.abs() + frame.interior_loss_j + frame.terminal_loss_j
            + frame.coupled_terminal_work_j.abs()).max(1e-30);
        assert!(frame.balance_residual_j().abs() < 1e-9 * scale);
        if let Ok(wall) = network.wall_state(4) {
            trace.wall_energy = trace.wall_energy.max(wall.stored_energy_j);
        }
    }
    assert!(trace.pressure.iter().any(|p| p.abs() > 1e-12));
    assert!(trace.work > 0.0);
    trace
}

#[test]
fn native_mimo_matches_direct_owners_with_walls_gases_mouth_and_multirate_receiver() {
    let text = native_source(true);
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let file = load_mimo(&text, ARTIFACT, &material, 37).unwrap();
    let gate = CancelGate::new_clock_free();
    let decoded = MimoRadiationArtifact::from_bytes_with_sources(ARTIFACT.as_bytes(), None, &gate).unwrap();
    let record = file.mimo_radiation_record().unwrap();
    assert_eq!(record.artifact_hash(), decoded.content_hash());
    assert_eq!(record.declaration_identity(), [0x11; 32]);
    assert_eq!(record.realization_identity(), [0x22; 32]);
    assert_eq!(record.path(), "flue-mutual.mimo-radiation");
    assert!(file.radiation_bindings().is_empty());
    assert_eq!(file.wall_records()[0].pack(), &material);
    let renderer = file.renderer();
    let voice = renderer.voice();
    let network = voice.network().unwrap();
    let environment = network.mimo_environment().unwrap();
    assert_eq!(record.gas_record().state(), &environment.source().exterior);
    assert_eq!(record.apertures(), environment.source().apertures.as_slice());
    assert_eq!(network.source().graph.sections.len(), 4);
    assert_eq!(network.source().graph.inlet_section_index().unwrap(), 1);
    assert_eq!(network.source().gases[1].temperature, 308.15);
    assert_eq!(network.source().gases[3].temperature, 298.15);
    assert_eq!(network.material_walls()[0].source().node, 4);
    assert_eq!(network.source().viscothermal.len(), 2);
    assert!(network.graph().sections.len() > network.source().graph.sections.len());
    assert_eq!(renderer.outlet_nodes(), &[1, 3]);
    assert_eq!(voice.mouth_contraction().unwrap().gas(), voice.gas());
    let windway = voice.unsteady_windway().unwrap();
    assert_eq!(windway.channel().gas(), voice.gas());
    let island = voice.island();
    let direct = FlueVoice::from_fractional_network(Some(file.card()),
        island.displacement_gain, island.seed_amplitude_rel, file.geometry(),
        network.source().clone(), network.wall_input().cloned(),
        network.fractional_sections().unwrap().0, 96_000, 3, 0, &gate).unwrap()
        .with_mouth_contraction(voice.mouth_contraction().unwrap().spec().clone()).unwrap()
        .with_unsteady_windway(windway.channel().spec().length_m,
            windway.channel().spec().transverse_cells, windway.convection_cells()).unwrap()
        .with_mimo_environment(decoded.source().clone(), environment.radiation().clone(), record.limits(), &gate).unwrap();
    let direct = FlueOutletPerformance::with_multirate_radiation(direct,
        renderer.schedule().clone(), *renderer.config(), renderer.radiation_input().unwrap().clone(), 96_000).unwrap();
    let expected = pressure(direct, 37);
    assert_eq!(pressure(file.into_renderer(), 37), expected);
    for block in [1, 137, 512] {
        assert_eq!(pressure(load_mimo(&text, ARTIFACT, &material, block).unwrap().into_renderer(), block), expected);
    }
    let trace = mutual_trace(load_mimo(&text, ARTIFACT, &material, 1).unwrap());
    assert!(trace.wall_energy > 0.0, "MIMO selection must retain and advance the material wall");
}

#[test]
fn matrix_material_and_coherent_atmosphere_data_reach_feedback_and_microphone() {
    let text = native_source(true);
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let baseline = mutual_trace(load_mimo(&text, ARTIFACT, &material, 1).unwrap());
    let diagonal = ARTIFACT.replace("residue 60000000 0 12000000 0 12000000 0 80000000 0",
        "residue 60000000 0 0 0 0 0 80000000 0")
        .replace("direct 600000 120000 120000 800000", "direct 600000 0 0 800000");
    let changed = mutual_trace(load_mimo(&text, &diagonal, &material, 1).unwrap());
    assert_ne!(changed.incoming, baseline.incoming, "off-diagonal reaction must return to the live jet");
    assert_ne!(changed.pressure, baseline.pressure);
    let changed = mutual_trace(load_mimo(&text, ARTIFACT,
        &pack(1500.0, 4e6, None, 100_000.0), 1).unwrap());
    assert_ne!(changed.incoming, baseline.incoming);
    assert_ne!(changed.pressure, baseline.pressure);

    // An independently supplied lower-pressure operator scales this synthetic
    // impedance at fixed T/composition. Import does not perform that scaling.
    let atmosphere = text.replace("101325", "81060");
    assert!(load_mimo(&atmosphere, ARTIFACT, &material, 37).is_err(), "stale gas must not be relabeled");
    let updated = ARTIFACT.replace("101325", "81060")
        .replace("residue 60000000 0 12000000 0 12000000 0 80000000 0",
            "residue 48000000 0 9600000 0 9600000 0 64000000 0")
        .replace("direct 600000 120000 120000 800000", "direct 480000 96000 96000 640000")
        .replace(&"22".repeat(32), &"44".repeat(32));
    let changed = load_mimo(&atmosphere, &updated, &material, 1).unwrap();
    let record = changed.mimo_radiation_record().unwrap();
    assert_eq!(record.gas_record().state().pressure, 81060.0);
    assert_eq!(record.realization_identity(), [0x44; 32]);
    assert_eq!(changed.renderer().radiation_input().unwrap().exterior, *record.gas_record().state());
    let changed = mutual_trace(changed);
    assert_ne!(changed.incoming, baseline.incoming);
    assert_ne!(changed.pressure, baseline.pressure);
}

#[test]
fn explicit_artifact_and_runtime_limits_refuse_and_legacy_loader_remains_exact() {
    let text = native_source(false);
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let gate = CancelGate::new_clock_free();
    assert!(FlueFilePerformance::from_bytes_with_sources(text.as_bytes(), 37, &gate,
        None, |_| Ok(material.clone())).is_err(), "old loader cannot silently omit requested MIMO feedback");
    for selection in [
        "radiation mimo ../outside.mimo-radiation 32 1048576 257",
        "radiation mimo flue-mutual.mimo-radiation 0 1048576 257",
        "radiation mimo flue-mutual.mimo-radiation 32 1 257",
        "radiation mimo flue-mutual.mimo-radiation 32 1048576 16",
        "radiation mimo flue-mutual.mimo-radiation 32 1048576 257 extra",
    ] {
        assert!(load_mimo(&text.replace(SELECTION, selection), ARTIFACT, &material, 37).is_err());
    }
    for artifact in [
        ARTIFACT.replace("aperture 0.08 0 0.004", "aperture 0.08 0 0.005"),
        ARTIFACT.replace("band 1 48000", "band 1 24000"),
        ARTIFACT.replace("pole real -3000", "pole real NaN"),
        ARTIFACT.replace("improper 0 0 0 0", "improper 1 0 0 0"),
        format!("{ARTIFACT}ignored\n"),
    ] {
        assert!(load_mimo(&text, &artifact, &material, 37).is_err());
    }
    let legacy = text.replace(SELECTION, "radiation prescribed");
    let expected = FlueFilePerformance::from_bytes(legacy.as_bytes(), 37, &gate).unwrap();
    let actual = FlueFilePerformance::from_bytes_with_all_sources(legacy.as_bytes(), 37, &gate,
        None, |_| Err("no material source expected".into()),
        |_| panic!("omitted MIMO must not load an artifact")).unwrap();
    assert!(actual.mimo_radiation_record().is_none());
    assert_eq!(pressure(expected.into_renderer(), 37), pressure(actual.into_renderer(), 37));
}

#[test]
fn native_cli_loads_relative_matrix_artifact_and_preserves_pcm_partition() {
    let text = native_source(false);
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let expected = pressure(load_mimo(&text, ARTIFACT, &material, 37).unwrap().into_renderer(), 37);
    let expected_wav = encode_pcm16_wav(&expected, 48_000, 1.0).unwrap().0;
    let dir = directory();
    let input = dir.join("mutual.performance");
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("flue-mutual.mimo-radiation"), ARTIFACT).unwrap();
    for block in [1, 137] {
        let output = dir.join(format!("mutual-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), expected_wav);
    }
    std::fs::write(&input, text.replace("101325", "81060")).unwrap();
    let output = dir.join("stale-mutual.wav");
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}
