//! One physical jet drives every branch, regional gas and radiating outlet.
use super::*;
use fs_couple::bernoulli_aperture::network::NetworkNode;
use fs_couple::flue_loop::{FlueVoice, network::performance::FlueOutletPerformance};

const NETWORK: &str = include_str!("../../../examples/flue-network.performance");
fn source() -> String {
    let head = NETWORK.split_once("\ngesture\n").unwrap().0
        .replace("output 24000 5", "output 1024 1");
    format!("{head}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t3e2\nevents\t0\n")
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
