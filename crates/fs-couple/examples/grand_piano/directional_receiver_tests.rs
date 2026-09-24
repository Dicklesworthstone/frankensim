//! Ideal microphone controls act on the real field, never on the piano state.
use super::*;
use super::super::{exterior_loading,exterior_audio,playback,performance,geometry};
fn text(patterns:&str)->String {
    let base=tests::specification().replace("moving,skin","moving,soundboard_skin")
        .replace("band-hz,40,400,17","band-hz,40,300,41")
        .replace("board-band-hz,400","board-band-hz,300")
        .replace("receiver-m,0.05,0.05,1","receiver-m,0.04,0.04,0.022");
    format!("{base}receiver-m,0.04,0.04,0.022\nreceiver-evaluation,near-field\n{patterns}")
}
fn scene(patterns:&str)->super::super::Scene {
    let (board,courses,_,_)=super::super::tests::small_source_inputs();
    let controls=playback::Controls::from_texts(&courses,None,None,None).unwrap();
    super::super::prepare_controlled_body(&board,courses,None,Specification::read(&text(patterns)).unwrap(),
        &playback::Options::default(),controls,true).unwrap()
}
#[test]
fn explicit_omnis_preserve_the_scalar_path_and_incomplete_patterns_refuse() {
    let a=scene("");let b=scene("receiver-pattern,0,1,0,0,-1\nreceiver-pattern,1,1,1,0,0\n");
    let x=a.boundary.sample(&a.spec).unwrap();let y=b.boundary.sample(&b.spec).unwrap();
    assert_eq!(x.values,y.values);assert_eq!(x.delays_s,y.delays_s);
    for bad in ["receiver-pattern,0,0.5,0,0,0\n", "receiver-pattern,0,0.5,0,0,2\n",
        "receiver-pattern,0,NaN,0,0,-1\n", "receiver-pattern,2,0.5,0,0,-1\n",
        "receiver-pattern,0,1.1,0,0,-1\n", "receiver-pattern,0,0.5,0,0\n",
        "receiver-pattern,0,0.5,0,0,-1\nreceiver-pattern,0,0.5,0,0,1\n"] {
        assert!(Specification::read(&text(bad)).is_err());
    }
    let selected=text("receiver-pattern,0,0.5,0,0,-1\n");
    assert!(Specification::read(&selected.replace("receiver-evaluation,near-field\n","")).is_err());
    let extra=format!("{}receiver-pattern,1,1,0,0,1\n",tests::specification());
    assert!(Specification::read(&extra).is_err()); // explicit omni still requires a receiver
    let mut invalid=Specification::read(&selected).unwrap();invalid.near_field_receivers=false;
    assert!(ReceiverSet::for_spec(&a.boundary,&invalid).is_err());
}
#[test]
fn opposite_microphones_observe_pressure_velocity_without_changing_radiation_impedance() {
    let selected=scene("receiver-pattern,0,0.5,0,0,-1\nreceiver-pattern,1,0.5,0,0,1\n");
    let omni=Specification::read(&text("")).unwrap();
    let samples=selected.boundary.sample(&selected.spec).unwrap();
    let dipoles=Specification::read(&text("receiver-pattern,0,0,0,0,-1\nreceiver-pattern,1,0,0,0,1\n")).unwrap();
    let mut changed=false;
    for f in [0,20,40] {
        let w=samples.omega[f];
        let directional=exterior_loading::sample(&selected.boundary,&selected.spec,w).unwrap();
        let scalar=exterior_loading::sample(&selected.boundary,&omni,w).unwrap();
        let dipole=exterior_loading::sample(&selected.boundary,&dipoles,w).unwrap();
        assert_eq!(directional.impedance,scalar.impedance);assert_eq!(dipole.impedance,scalar.impedance);
        for input in 0..selected.piano.bank.board_count {
            let p=scalar.receiver_transfer[0][input];
            let front=directional.receiver_transfer[0][input];let back=directional.receiver_transfer[1][input];
            // Coincident opposite cardioids sum to the pressure field. The
            // remaining signed difference is a true vector-field observation.
            assert!((front+back-p).abs()<1e-6*p.abs().max(1e-12));
            assert!((dipole.receiver_transfer[0][input]+dipole.receiver_transfer[1][input]).abs()
                <1e-12*dipole.receiver_transfer[0][input].abs().max(1e-12));
            changed|=(front-back).abs()>1e-6*p.abs();
            for channel in 0..2 {
                let converted=samples.values[channel][input][f]*C64::new(0.,-w);
                assert!((converted-directional.receiver_transfer[channel][input]).abs()
                    <1e-8*(1.+directional.receiver_transfer[channel][input].abs()));
            }
        }
    }
    assert!(changed);assert!(selected.piano.bank.q.iter().chain(&selected.piano.bank.v).all(|v|*v==0.));
}
#[test]
fn source_hammer_directional_stereo_retains_one_clock_and_identical_reacted_mechanics() {
    let score=||performance::Performance::read("sample,event,key,value\n0,note_on,69,0.5\n1200,note_off,69,0\n",&[69],2400).unwrap();
    for feedback in [false,true] {
        let mut selected=scene("receiver-pattern,0,0.5,0,0,-1\nreceiver-pattern,1,0.5,0,0,-1\n");
        let mut manual=scene("");
        let (baked,samples,_)=super::super::bake(&mut selected,feedback).unwrap();
        let (_,reference,_)=super::super::bake(&mut manual,feedback).unwrap();
        assert_eq!(samples.delays_s,reference.delays_s);
        assert_ne!(samples.values,reference.values);
        let audio=exterior_audio::render(&mut selected.piano,score(),2400,&baked,2.).unwrap();
        let mut schedule=score();
        for frame in 0..2400 {schedule.dispatch(frame,&mut manual.piano).unwrap();manual.piano.step().unwrap();}
        assert_eq!(selected.piano.bank.q,manual.piano.bank.q);assert_eq!(selected.piano.bank.v,manual.piano.bank.v);
        assert_eq!(selected.piano.radiation_energy_j(),manual.piano.radiation_energy_j());
        assert_eq!(selected.piano.accounting.radiation_loss_j,manual.piano.accounting.radiation_loss_j);
        assert!(selected.piano.accounting.felt_loss_j>0.);assert!(audio.peak_pa>1e-14);
        assert!((selected.piano.accounting.input_work_j-selected.piano.energy_j()-selected.piano.accounting.dissipated_j()).abs()<1e-7);
        assert_eq!(u16::from_le_bytes([audio.wav[22],audio.wav[23]]),2);
        let start=audio.wav.windows(4).position(|w|w==b"data").unwrap()+8;
        for frame in audio.wav[start..].chunks_exact(4) {assert_eq!(&frame[..2],&frame[2..]);}
    }
}
#[test]
fn directional_specs_reach_both_harmonic_commands_and_bad_axes_publish_nothing() {
    use std::io::Write;
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir=std::env::temp_dir().join(format!("fs-piano-directional-{}-{stamp}",std::process::id()));
    std::fs::create_dir(&dir).unwrap();let name=|p:&str|dir.join(p).to_str().unwrap().to_owned();
    let (board,courses,_,_)=super::super::tests::small_source_inputs();
    let valid=text("receiver-pattern,0,0.5,0,0,-1\nreceiver-pattern,1,0,0,0,1\n");
    for (path,contents) in [("board.fsb",board),("scale.csv",geometry::write_scale(&courses)),
        ("mics.fspe",valid.clone()),("bad.fspe",valid.replace("0.5,0,0,-1","0.5,0,0,-2"))] {
        std::fs::OpenOptions::new().create_new(true).write(true).open(name(path)).unwrap().write_all(contents.as_bytes()).unwrap();
    }
    for command in ["response","admittance"] {
        let mut args=vec![command.into(),name("board.fsb"),name("scale.csv"),"board-skin-continuous".into(),name("mics.fspe")];
        if command=="admittance" {args.push("69".into());}
        let output=name(&format!("{command}.csv"));args.push(output.clone());
        super::super::run(&args).unwrap();let csv=std::fs::read_to_string(&output).unwrap();
        assert!(csv.contains("ideal first-order"));assert!(csv.contains("Pa-equivalent"));assert!(!csv.contains("NaN"));
        assert!(super::super::run(&args).is_err());assert_eq!(csv,std::fs::read_to_string(&output).unwrap());
    }
    let bad=vec!["render-loaded".into(),name("board.fsb"),name("scale.csv"),"board-skin-continuous".into(),name("bad.fspe"),name("refused.wav"),"0.05".into()];
    assert!(super::super::run(&bad).is_err());assert!(!std::path::Path::new(&name("refused.wav")).exists());
}
