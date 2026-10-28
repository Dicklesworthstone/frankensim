//! One physical jet drives every branch, regional gas and radiating outlet.
use super::*;
use fs_couple::bernoulli_aperture::network::NetworkNode;
use fs_couple::flue_loop::{FlueVoice, network::performance::FlueOutletPerformance};

#[path = "network/mouth_reaction.rs"]
mod mouth_reaction;

const NETWORK: &str = include_str!("../../../examples/flue-network.performance");
fn source() -> String {
    let head = NETWORK.split_once("\ngesture\n").unwrap().0
        .replace("output 24000 5", "output 1024 1");
    format!("{head}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t3e2\nevents\t0\n")
}

#[test]
fn independent_aperture_card_reaches_pcm_and_reports_actual_geometry() {
    let baseline = articulated_source();
    let text = baseline.replace("mouth_radiation -0.3 0 4 4",
        "mouth_aperture rectangle 0.010 0.012\nmouth_radiation -0.3 0 4 4");
    let gate = CancelGate::new_clock_free();
    let file = FlueFilePerformance::from_bytes(text.as_bytes(), 37, &gate).unwrap();
    let original = FlueFilePerformance::from_bytes(baseline.as_bytes(), 37, &gate).unwrap();
    assert_eq!(file.geometry(), original.geometry());
    assert_eq!(file.renderer().mouth_radiation().unwrap().aperture_size_m, Some([0.010, 0.012]));
    let pressure = physical(&text).0;
    assert_ne!(pressure, physical(&baseline).0);
    let dir = directory();
    let input = dir.join("independent-aperture.performance");
    let output = dir.join("independent-aperture.wav");
    std::fs::write(&input, &text).unwrap();
    success(run(&input, &output, 37, false));
    assert_eq!(std::fs::read(&output).unwrap(), encode_pcm16_wav(&pressure, 48_000, 1.0).unwrap().0);
    let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(sidecar.contains("\"size_m\":[1e-2,1.2e-2]"));
    for invalid in [
        text.replace("rectangle 0.010 0.012", "rectangle 0 0.012"),
        text.replace("rectangle 0.010 0.012", "rectangle NaN 0.012"),
        text.replace("rectangle 0.010 0.012", "ellipse 0.010 0.012"),
        text.replace("mouth_radiation -0.3 0 4 4\n", ""),
    ] {
        assert!(FlueFilePerformance::from_bytes(invalid.as_bytes(), 37, &gate).is_err());
    }
}

#[test]
fn physical_network_and_actual_inlet_gas_survive_file_lowering_and_cli() {
    let text = source();
    let gate = CancelGate::new_clock_free();
    let file = FlueFilePerformance::from_bytes(text.as_bytes(), 37, &gate).unwrap();
    let p = file.renderer();
    let n = p.voice().network().unwrap();
    assert_eq!(n.source().graph.inlet_section_index().unwrap(), 1);
    assert_eq!(n.source().graph.nodes[2], NetworkNode::Inlet);
    assert_eq!(n.source().gases[1].temperature, 308.15);
    assert_eq!(n.source().gases[2].temperature, 298.15);
    assert_eq!(p.voice().gas(), &n.source().gases[1]);
    assert_eq!(p.voice().viscous_windway().unwrap().gas(), Some(&n.source().gases[1]));
    assert_eq!(p.radiation_input().unwrap().exterior.temperature, 293.15);
    assert_eq!(p.outlet_nodes(), &[1,3]);
    assert_eq!(n.source().viscothermal[0].section, 1);
    assert!(n.graph().sections.len() > n.source().graph.sections.len());
    // The file constructor composes the original source, propagation and
    // multirate receiver. Rebuilding those owners explicitly is bit-identical.
    let island = p.voice().island();
    let direct = FlueVoice::from_fractional_network(Some(file.card()),
        island.displacement_gain,island.seed_amplitude_rel,file.geometry(),
        n.source().clone(),None,n.fractional_sections().unwrap().0,
        96000,3,0,&gate).unwrap().with_viscous_windway(0.015).unwrap();
    let mut direct = FlueOutletPerformance::with_multirate_radiation(direct,
        p.schedule().clone(),*p.config(),p.radiation_input().unwrap().clone(),96000).unwrap();
    let mut expected = vec![0.0; 1024];
    for block in expected.chunks_mut(37) { direct.block(block).unwrap(); }
    assert_eq!(physical(&text).0, expected);
    let dir = directory();
    let input = dir.join("network.performance");
    std::fs::write(&input,&text).unwrap();
    let wav = encode_pcm16_wav(&expected,48000,1.0).unwrap().0;
    for block in [1,137] {
        let output = dir.join(format!("network-{block}.wav"));
        success(run(&input,&output,block,false));
        assert_eq!(std::fs::read(&output).unwrap(),wav);
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"duct_network\"","\"inlet_section\":1","\"node\":3",
            "\"kind\":\"impedance\"","\"temperature_k\":3.0815e2"] {
            assert!(sidecar.contains(token),"{token}: {sidecar}");
        }
    }
}

#[test]
fn branch_medium_geometry_and_outlet_selection_change_actual_pressure() {
    let text = source();
    let baseline = physical(&text).0;
    let energy: f64 = baseline.iter().map(|p| p*p).sum();
    for changed in [
        text.replace("0 3 0.04 0.004","0 3 0.07 0.004"),
        text.replace("gas 298.15 101325 0","gas mixture 298.15 101325 mole 2 N2 0.5 CO2 0.5"),
        text.replace("outlets 1 1\noutlet 3 0.08 0","outlets 1 0"),
    ] {
        let pressure = physical(&changed).0;
        let difference: f64 = pressure.iter().zip(&baseline).map(|(a,b)| (a-b).powi(2)).sum();
        assert!(difference > 1e-10*energy,"physical branch/outlet data must change receiver pressure");
    }
}

#[test]
fn incompatible_graph_gas_and_receiver_inputs_refuse_before_any_cli_output() {
    let text = source();
    let dir = directory();
    let input = dir.join("refused-network.performance");
    for (i, changed) in [
        text.replace("gas 298.15 101325 0","gas 298.15 80000 0"),
        text.replace("duct_node inlet","duct_node junction"),
        text.replace("duct_section 0 3","duct_section 0 9"),
        text.replace("outlet 3 0.08 0","outlet 0 0.08 0"),
        text.replace("outlet 3 0.08 0","outlet 1 0.08 0"),
        text.replace("outlet 3 0.08 0","outlet 3 0 0"),
        text.replace("loss_section 1 4","loss_section 9 4"),
        text.replace("impedance 1000000 100 0","impedance -1 100 0"),
        text.replace("outlets 1 1\noutlet 3 0.08 0\n",""),
    ].iter().enumerate() {
        std::fs::write(&input,changed).unwrap();
        let output = dir.join(format!("refused-{i}.wav"));
        assert!(!run(&input,&output,37,false).status.success());
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}

#[test]
fn unsteady_windway_data_reaches_regional_source_and_physical_cli_output() {
    let text = source().replace("windway slit 0.015", "windway unsteady 0.015 32 64");
    let file = FlueFilePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new_clock_free()).unwrap();
    let voice = file.renderer().voice();
    let field = voice.unsteady_windway().unwrap();
    assert_eq!(field.channel().gas(), &voice.network().unwrap().source().gases[1]);
    assert_eq!(field.channel().spec().length_m, 0.015);
    assert_eq!(field.channel().spec().transverse_cells, 32);
    assert_eq!(field.convection_cells(),64);
    assert_eq!(voice.jet_speed_m_s(),0.0);
    let expected = physical(&text).0;
    assert_ne!(expected,physical(&source()).0,"resolved source inertia must reach the receiver");
    let wav = encode_pcm16_wav(&expected,48000,1.0).unwrap().0;
    let dir = directory();
    let input = dir.join("unsteady.performance");
    std::fs::write(&input,&text).unwrap();
    for block in [1,137] {
        let output = dir.join(format!("unsteady-{block}.wav"));
        success(run(&input,&output,block,false));
        assert_eq!(std::fs::read(&output).unwrap(),wav);
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"windway_model\":\"unsteady-no-slip-slit-v1\"",
            "\"windway_transverse_cells\":32","\"jet_convection_cells\":64",
            "\"inlet_section\":1"] {
            assert!(sidecar.contains(token),"{token}: {sidecar}");
        }
    }
    for (i, row) in ["windway unsteady 0.015 1 64","windway unsteady 0.015 32 4097",
        "windway unsteady 0.015 32","windway unsteady -0.015 32 64"].iter().enumerate() {
        std::fs::write(&input,text.replace("windway unsteady 0.015 32 64",row)).unwrap();
        let output = dir.join(format!("refused-unsteady-{i}.wav"));
        assert!(!run(&input,&output,37,false).status.success());
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}

#[test]
fn developing_windway_data_holds_the_regional_jet_and_reports_its_validity_measures() {
    let text = source().replace("windway slit 0.015", "windway slit-developing 0.015");
    let file = FlueFilePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new_clock_free()).unwrap();
    let voice = file.renderer().voice();
    let (law, held) = voice.developing_windway().unwrap();
    assert_eq!(law.gas(), &voice.network().unwrap().source().gases[1]);
    assert_eq!(law.length_m(), 0.015);
    assert_eq!(held.pressure_pa, 3e2, "validity binds the schedule's held pressure");
    assert!(voice.viscous_windway().is_none() && voice.unsteady_windway().is_none());
    held.validity().unwrap();
    let expected = physical(&text).0;
    assert_ne!(expected,physical(&source()).0,"the entrance-shear speed must reach the receiver");
    let wav = encode_pcm16_wav(&expected,48000,1.0).unwrap().0;
    let dir = directory();
    let input = dir.join("developing.performance");
    std::fs::write(&input,&text).unwrap();
    for block in [1,137] {
        let output = dir.join(format!("developing-{block}.wav"));
        success(run(&input,&output,block,false));
        assert_eq!(std::fs::read(&output).unwrap(),wav);
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        let measures = [format!("\"x_plus\":{:e}",held.x_plus),
            format!("\"reynolds_hydraulic\":{:e}",held.reynolds_hydraulic),
            format!("\"reynolds_length\":{:e}",held.reynolds_length)];
        for token in ["\"windway_model\":\"developing-slit-shah1978-v1\"",
            "\"windway_length_m\":1.5e-2","\"pressure_pa\":3e2",
            "\"attack_covered\":false","\"inlet_section\":1"].map(String::from).into_iter().chain(measures) {
            assert!(sidecar.contains(&token),"{token}: {sidecar}");
        }
    }
    // A 0.5 mm slit holds x+ near 1e-4, below Shah's stated range; a 1 Pa note holds
    // Re_Dh near 34, where Chen's inlet term moves the drop by about 5 %.
    for (i, changed) in [
        text.replace("slit-developing 0.015","slit-developing 0.0005"),
        text.replace("initial\tpressure-pa\t3e2","initial\tpressure-pa\t1e0"),
        text.replace("initial\tpressure-pa\t3e2","initial\tpressure-pa\t0e0"),
        text.replace("slit-developing 0.015","slit-developing 0"),
        text.replace("slit-developing 0.015","slit-developing -0.015"),
        text.replace("slit-developing 0.015","slit-developing"),
    ].iter().enumerate() {
        std::fs::write(&input,changed).unwrap();
        let output = dir.join(format!("refused-developing-{i}.wav"));
        assert!(!run(&input,&output,37,false).status.success());
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}

#[test]
fn current_sample_windway_feedback_retains_regional_gas_pressure_and_pcm() {
    // A manufactured short cut-up makes transit less than the mechanical
    // sample period. Its current-mouth term must use the same physical
    // source, multirate receiver and PCM owner as the ordinary geometry.
    let text = source().replace("windway slit 0.015", "windway unsteady 0.015 32 64")
        .replace("flue 0.010 0.001 0.004", "flue 0.010 0.001 0.00002")
        .replace("jet 90 1e-6", "jet 0.01 1e-6");
    let changed = text.replace("N2 0.77 O2 0.2 CO2 0.03", "N2 0.5 O2 0.2 CO2 0.3");
    let gate = CancelGate::new_clock_free();
    let render = |input: &str, block: usize| {
        let mut renderer = FlueFilePerformance::from_bytes(input.as_bytes(), block, &gate)
            .unwrap().into_renderer();
        let mut pressure = vec![0.0; 1024];
        for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
        let v = renderer.voice();
        let field = v.unsteady_windway().unwrap();
        assert_eq!(field.channel().gas(), &v.network().unwrap().source().gases[1]);
        assert!(v.jet_speed_m_s()*v.island().convection_ratio/96_000.0 > 0.00002,
            "the fixture must exercise a current-mouth outlet coefficient");
        assert!(pressure.iter().all(|p| p.is_finite()));
        assert!(pressure.iter().any(|p| p.abs() > 1e-10));
        pressure
    };
    let expected = render(&text, 37);
    assert_eq!(expected, render(&text, 1));
    assert_eq!(expected, render(&text, 137));
    let altered = render(&changed, 37);
    assert_ne!(expected, altered, "source-gas data must reach actual receiver pressure");
    let wav = encode_pcm16_wav(&expected, 48_000, 1.0).unwrap().0;
    assert_ne!(wav, encode_pcm16_wav(&altered, 48_000, 1.0).unwrap().0);
    let dir = directory();
    let input = dir.join("current-windway.performance");
    let output = dir.join("current-windway.wav");
    std::fs::write(&input, &text).unwrap();
    success(run(&input, &output, 137, false));
    assert_eq!(std::fs::read(output).unwrap(), wav);
}


fn articulated_source() -> String {
    use fs_scenario::gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue};
    let head = include_str!("../../../examples/flue-embouchure.performance")
        .split_once("\ngesture\n").unwrap().0.replace("output 24000 5", "output 1024 1");
    let schedule = GestureSchedule::try_new(700, vec![
        GestureTrack { id: "blow".into(), target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(300.0), events: vec![] },
        GestureTrack { id: "angle".into(), target: GestureTarget::JetAngle,
            initial: GestureValue::AngleRad(0.0), events: vec![
                GestureEvent { time_s: 1.0/700.0, transition_s: 4.0/700.0,
                    value: GestureValue::AngleRad(0.035) },
                GestureEvent { time_s: 6.0/700.0, transition_s: 2.0/700.0,
                    value: GestureValue::AngleRad(-0.015) },
            ] },
    ]).unwrap();
    format!("{head}\ngesture\n{}", String::from_utf8(schedule.to_canonical_bytes()).unwrap())
}

#[test]
fn mouth_separation_and_angle_data_share_the_original_fluid_runtime_and_cli() {
    let text = articulated_source();
    let gate = CancelGate::new_clock_free();
    let file = FlueFilePerformance::from_bytes(text.as_bytes(), 37, &gate).unwrap();
    let p = file.renderer();
    let v = p.voice();
    let n = v.network().unwrap();
    let loss = v.mouth_contraction().unwrap();
    assert_eq!(loss.gas(), &n.source().gases[1]);
    assert_eq!(loss.width_m(), file.geometry().flue_width_m);
    assert_eq!(loss.opening_m(), file.geometry().cut_up_m);
    assert_eq!(loss.spec().provenance, "authored-fixed-window-demonstrator");
    assert_eq!(loss.spec().area_ratio, 0.85);
    assert_eq!(v.unsteady_windway().unwrap().channel().gas(), loss.gas());
    let island = v.island();
    let direct = FlueVoice::from_fractional_network(Some(file.card()),
        island.displacement_gain, island.seed_amplitude_rel, file.geometry(),
        n.source().clone(), None, n.fractional_sections().unwrap().0,
        96_000, 3, 0, &gate).unwrap()
        .with_mouth_contraction(loss.spec().clone()).unwrap()
        .with_unsteady_windway(0.015, 32, 64).unwrap();
    let mut direct = FlueOutletPerformance::with_multirate_radiation(direct,
        p.schedule().clone(), *p.config(), p.radiation_input().unwrap().clone(), 96_000).unwrap()
        .with_mouth_radiation(p.mouth_radiation().unwrap()).unwrap();
    let mut expected = vec![0.0; 1024];
    for block in expected.chunks_mut(37) { direct.block(block).unwrap(); }
    assert_eq!(physical(&text).0, expected);
    assert_eq!(direct.voice().jet_angle_rad(), -0.015);
    assert_eq!(direct.voice().samples_rendered(), 2048);
    let separation: f64 = direct.voice().diagnostics().iter().map(|d| d.separation_loss_j).sum();
    assert!(separation > 0.0, "the declared mouth law must be advanced in the source solve");
    assert!(direct.voice().unsteady_windway().unwrap().channel().velocities_m_s()
        .iter().any(|v| *v > 0.0));
    assert!(direct.voice().network().unwrap().stored_energy_j() > 0.0);
    let wav = encode_pcm16_wav(&expected, 48_000, 1.0).unwrap().0;
    let dir = directory();
    let input = dir.join("embouchure.performance");
    std::fs::write(&input, &text).unwrap();
    for block in [1, 137] {
        let output = dir.join(format!("embouchure-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), wav);
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert!(sidecar.contains("\"mouth_contraction\":{"));
        assert!(sidecar.contains("\"mouth_radiation\":{\"shape\":\"rectangle\""));
        assert!(sidecar.contains("\"source\":\"signed outward inlet volume flow\""));
        assert!(sidecar.contains("\"source\":\"authored-fixed-window-demonstrator\""));
        assert!(sidecar.contains("\"windway_model\":\"unsteady-no-slip-slit-v1\""));
    }
}

#[test]
fn mouth_source_can_be_selected_from_data_without_changing_the_jet_or_terminal_loads() {
    let text = articulated_source();
    let without = text.replace("mouth_radiation -0.3 0 4 4\n", "");
    let with_pressure = physical(&text).0;
    let without_pressure = physical(&without).0;
    let difference: f64 = with_pressure.iter().zip(&without_pressure)
        .map(|(a,b)| (a-b).powi(2)).sum();
    assert!(difference > 1e-12, "the source mouth must contribute actual pressure");
    let gate = CancelGate::new_clock_free();
    for invalid in ["0 0 4 4", "-0.3 0 0 4", "-0.3 0 513 4", "NaN 0 4 4"] {
        let invalid = text.replace("mouth_radiation -0.3 0 4 4",
            &format!("mouth_radiation {invalid}"));
        assert!(FlueFilePerformance::from_bytes(invalid.as_bytes(), 37, &gate).is_err());
    }
}

#[test]
fn data_only_angle_separation_and_atmosphere_changes_reach_microphone_pressure() {
    let text = articulated_source();
    let baseline = physical(&text).0;
    for changed in [
        text.replace("mouth_loss contraction 0.85 0.3 authored-fixed-window-demonstrator", "mouth_loss none"),
        text.replace("angle-rad\t3.5e-2", "angle-rad\t1e-2"),
        text.replace("O2 0.2 CO2 0.03", "O2 0.15 CO2 0.08"),
    ] {
        let pressure = physical(&changed).0;
        assert!(pressure.iter().zip(&baseline).any(|(a,b)| (a-b).abs() > 1.0e-15),
            "authored physical controls and gas must reach the microphone");
    }
    let no_loss = text.replace("mouth_loss contraction 0.85 0.3 authored-fixed-window-demonstrator\n", "");
    assert_eq!(physical(&no_loss).0, physical(&no_loss.replace("jet_card ",
        "mouth_loss none\njet_card ")).0);
    for row in [
        "mouth_loss contraction 0 0.3 source",
        "mouth_loss contraction 1.01 0.3 source",
        "mouth_loss contraction 0.85 0.31 source",
        "mouth_loss contraction 0.85 0.3",
        "mouth_loss contraction 0.85 0.3 source extra",
    ] {
        let bad = text.replace("mouth_loss contraction 0.85 0.3 authored-fixed-window-demonstrator", row);
        assert!(FlueFilePerformance::from_bytes(bad.as_bytes(), 37, &CancelGate::new_clock_free()).is_err());
    }
    let bad_angle = text.replace("angle-rad\t3.5e-2", "angle-rad\t2e0");
    assert!(FlueFilePerformance::from_bytes(bad_angle.as_bytes(), 37, &CancelGate::new_clock_free()).is_err());
    let dir = directory();
    let input = dir.join("bad-angle.performance");
    let output = dir.join("bad-angle.wav");
    std::fs::write(&input, bad_angle).unwrap();
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}
