//! Complete supplied networks through wind/ensemble, not only library construction.
use super::*;
const BRANCHED:&str=include_str!("../../examples/plate-valve-branched.performance");

#[test]
fn supplied_branches_walls_and_side_chambers_reach_exact_native_and_decimated_pcm() {
    let dir=directory();let input=dir.join("network.performance");
    for (i,text) in [BRANCHED.to_string(),BRANCHED.replace("audio 96000 9602","audio 48000 4801")].iter().enumerate() {
        std::fs::write(&input,text).unwrap();let (pressure,_)=reference(text);
        let (wav,clips)=encode_pcm16_wav(&pressure,48000,0.01).unwrap();
        assert!(wav[44..].iter().any(|b|*b!=0));
        for block in [1,37,512] {
            let output=dir.join(format!("network-{i}-{block}.wav"));success(run(&input,&output,block,i==0));
            assert_eq!(std::fs::read(&output).unwrap(),wav);
            let sidecar=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            for token in ["\"node_count\":5","\"section_count\":4","\"radiation_terminal_count\":1",
                "\"observed_node\":3","\"node\":2,\"kind\":\"shunt\"","\"node\":4,\"kind\":\"impedance\"",
                "one_way_mechanical_samples","total_represented_section_length_m",
                "\"radiation_feedback_added\":true","exterior output selects one outlet"] {
                assert!(sidecar.contains(token),"{token}: {sidecar}");
            }
            assert!(sidecar.contains(&format!("\"clipped_samples\":{clips}")));
        }
        let relocated=dir.join(format!("relocated-{i}.data"));std::fs::write(&relocated,text).unwrap();
        let output=dir.join(format!("relocated-{i}.wav"));success(run(&relocated,&output,37,i==0));
        assert_eq!(std::fs::read(output.with_extension("provenance.json")).unwrap(),
            std::fs::read(dir.join(format!("network-{i}-37.provenance.json"))).unwrap());
        assert!(!run(&relocated,&output,37,i==0).status.success());assert_eq!(std::fs::read(&output).unwrap(),wav);
    }
    // Explicit node observations can address side chambers, not just node 1.
    let text=BRANCHED.replace("observation network-baffled 3 0 0 0.2 8 32 1500","observation network-node 4");
    std::fs::write(&input,&text).unwrap();let output=dir.join("cavity-pressure.wav");
    success(run(&input,&output,37,true));
    assert_eq!(std::fs::read(&output).unwrap(),encode_pcm16_wav(&reference(&text).0,48000,0.01).unwrap().0);
    let sidecar=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(sidecar.contains("\"pressure_point\":\"network-node\"") && sidecar.contains("\"observed_node\":4"));
    assert!(!sidecar.contains("\"outlet_receiver\""));
}

#[test]
fn branched_ensemble_keeps_every_source_clock_and_has_no_hidden_input_scale_gain() {
    let dir=directory();let high=dir.join("high.performance");let low=dir.join("low.performance");
    let low_text=BRANCHED.replace("audio 96000 9602","audio 48000 4801")
        .replace("network-baffled 3 0 0 0.2","network-baffled 3 0.1 0 0.4");
    std::fs::write(&low,&low_text).unwrap();
    let (a,da)=reference(BRANCHED);let (b,db)=reference(&low_text);
    let mut pressure=a.clone();
    for (x,y) in pressure[da-db..].iter_mut().zip(&b) {*x+=y;}
    let expected=encode_pcm16_wav(&pressure,48000,0.02).unwrap().0;
    for (i,text) in [BRANCHED.to_string(),BRANCHED.replace("9602 0.01","9602 100")].iter().enumerate() {
        std::fs::write(&high,text).unwrap();let output=dir.join(format!("branch-mix-{i}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--valve").arg(&high).arg("--valve").arg(&low)
            .args(["--decimate","--block","37","--full-scale-pa","0.02"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(),expected);
        let sidecar=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert_eq!(sidecar.matches("\"duct_network\":").count(),2);
        assert_eq!(sidecar.matches("\"source_pcm_scale_applied\":false").count(),2);
        assert_eq!(sidecar.matches("\"propagation_delay_mechanical_samples\":").count(),2);
    }
}

#[test]
fn invalid_graphs_or_ambiguous_microphones_never_create_output_or_fall_back_to_a_tube() {
    let dir=directory();let input=dir.join("bad.performance");
    for (i,text) in [
        BRANCHED.replace("network 5 4","network 0 4"),
        BRANCHED.replace("network 5 4 1048576","network 5 4 8"),
        BRANCHED.replace("duct_node junction","duct_node inlet"),
        BRANCHED.replace("duct_section 1 4","duct_section 1 0"),
        BRANCHED.replace("wall 0.001 0.1","wall 0.001 NaN"),
        BRANCHED.replace("cavity 0.00001","cavity -1"),
        BRANCHED.replace("network-baffled 3","network-baffled 4"),
        BRANCHED.replace("observation network-baffled 3 0 0 0.2 8 32 1500","observation terminal"),
        BRANCHED.replace("8 32 1500","8 32 1600"),
    ].iter().enumerate() {
        std::fs::write(&input,text).unwrap();let output=dir.join(format!("refused-{i}.wav"));
        assert!(!run(&input,&output,37,true).status.success());
        assert!(!output.exists() && !output.with_extension("provenance.json").exists());
    }
}
