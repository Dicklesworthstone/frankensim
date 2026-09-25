use super::*;
const INPUT: &str = include_str!("../close-drum.frm");
const FILE: &str = concat!(env!("CARGO_MANIFEST_DIR"),"/examples/percussion/close-drum.frm");

#[test]
fn independent_patterns_positions_and_omni_defaults_are_retained() {
    let spec=Spec::parse(INPUT).unwrap();assert_eq!(spec.receivers.len(),2);
    for (i,receiver) in spec.receivers.iter().enumerate() {
        let Receiver::NearField(m)=receiver else{panic!("must prepare actual near field")};
        assert_eq!(m.position_m,[0.04,0.0,if i==0{0.105}else{-0.105}]);
        assert_eq!(m.pattern.pressure_fraction(),0.5);
        assert_eq!(m.pattern.front_axis(),[0.0,0.0,if i==0{-1.0}else{1.0}]);
        assert_eq!(m.minimum_clearance_m,0.001);
    }
    let spec=Spec::parse("frankensim-microphones-v1\nnear_field,0.001,0.000001,14,2000000\nreceiver,0,0,0.1").unwrap();
    let Receiver::NearField(m)=spec.receivers[0] else{panic!()};
    assert_eq!(m.pattern,FirstOrder::default());
    let text=INPUT.replace("pattern,1,0.5,0,0,1","pattern,1,0,1,0,0");
    let spec=Spec::parse(&text).unwrap();
    let Receiver::NearField(m)=spec.receivers[1] else{panic!()};
    assert_eq!(m.pattern.pressure_fraction(),0.0);assert_eq!(m.pattern.front_axis(),[1.0,0.0,0.0]);
}

#[test]
fn incomplete_ambiguous_nonfinite_or_overbudget_microphones_refuse() {
    for text in [String::new(),INPUT.replace("frankensim-microphones-v1","wrong"),
        INPUT.replace("near_field,0.001,0.000001,14,2000000",""),
        INPUT.replace("near_field,0.001","near_field,0"),
        INPUT.replace("0.000001,14","0.1,14"),
        INPUT.replace("14,2000000","17,2000000"),
        INPUT.replace("14,2000000","14,79"),
        INPUT.replace("14,2000000","14,20000001"),
        INPUT.replace("receiver,0.04,0,0.105","receiver,NaN,0,0.105"),
        INPUT.replace("receiver,0.04,0,0.105","receiver,1e10,0,0.105"),
        INPUT.replace("pattern,1,0.5","pattern,2,0.5"),
        INPUT.replace("pattern,1,0.5","pattern,1,1.01"),
        INPUT.replace("pattern,1,0.5,0,0,1","pattern,1,0.5,0,0,2"),
        INPUT.replace("pattern,1,0.5,0,0,1","pattern,1,0.5,0,0,NaN"),
        INPUT.replace("receiver,0.04,0,-0.105", ""),
        format!("{INPUT}\nreceiver,0,0,0.2"),format!("{INPUT}\npattern,0,1,0,0,1"),
        format!("{INPUT}\nnear_field,0.001,0.000001,14,2000000"),
        format!("{INPUT}\ngain,2"),format!("{INPUT}{}"," ".repeat(MAX_BYTES))] {
        assert!(Spec::parse(&text).is_err(),"accepted invalid microphone input {text:?}");
    }
}

#[test]
fn option_and_actual_commands_reject_conflicts_before_instrument_construction() {
    let mut args=vec!["drum-mic".into(),"--microphone-spec".into(),FILE.into(),"1".into()];
    assert_eq!(option(&mut args).unwrap().unwrap().receivers.len(),2);
    assert_eq!(args,["drum-mic","1"]);
    for mut args in [vec!["--microphone-spec".into()],
        vec!["--microphone-spec".into(),"--radiation-spec".into()],
        vec!["--microphone-spec".into(),FILE.into(),"--microphone-spec".into(),FILE.into()]] {
        let before=args.clone();assert!(option(&mut args).is_err());assert_eq!(args,before);
    }
    for command in ["splash-mic","drum-mic","drum-stretch-mic","drum-modal-mic","snare-mic","snare-off-mic","hihat-mic"] {
        let spec=Spec::parse(INPUT).unwrap();assert!(spec.admit_command(command,false,false).is_ok());
        assert!(spec.admit_command(command,true,false).is_err());assert!(spec.admit_command(command,false,true).is_err());
    }
    for command in ["splash","drum","drum-wav","snare","hihat","hihat-wav"] {
        let mut args=vec![command.to_string()];
        if command.starts_with("hihat") {args.push("missing-instrument.fshh".into());}
        args.extend(["--microphone-spec".into(),FILE.into()]);
        let error=crate::run_args(args).unwrap_err().to_string();
        assert!(error.contains("requires an existing finite-point"),"{command}: {error}");
    }
    let error=crate::run_args(vec!["drum-mic".into(),"1".into(),"20".into(),"0".into(),"0".into(),"0.1".into(),
        "--microphone-spec".into(),FILE.into()]).unwrap_err().to_string();
    assert!(error.contains("supplies every receiver"),"{error}");
    let error=crate::hihat::run(vec!["hihat-mic".into(),"missing-instrument.fshh".into(),
        "--microphone-spec".into(),FILE.into(),"--microphone-right".into(),"0,0,0.1".into()]).unwrap_err().to_string();
    assert!(error.contains("supplies every receiver"),"must refuse conflict before reading missing instrument: {error}");
}

#[test]
fn file_selected_near_microphones_render_real_bem_fields_with_one_mechanical_clock() {
    use crate::{Experiment,Mechanics,acoustics::{self,Boundary,MECHANICAL_DT,SUBSTEPS}};
    use fs_bem::panel3d::SpherePanels;
    use fs_couple::render::plate::impact::{ImpactBody,BodyPotential,ImpactSystem};
    use fs_exec::CancelGate;
    // Prescribed breathing on a real closed surface; no specimen/frequency claim.
    let frames=96;
    let make=||{
        let surface=SpherePanels::icosphere(0.02,0).unwrap();
        let boundary=Boundary{triangles:surface.triangles().unwrap().to_vec(),
            weights:vec![vec![100.0;surface.areas().len()]],state_modes:vec![0]};
        let mut body=ImpactBody::free_mass(1.0,0.0,0.0001).unwrap().0;
        body.potential=BodyPotential::Linear(vec![1000.0]);
        let system=ImpactSystem::new(vec![body],vec![],vec![],vec![],
            crate::config((frames*SUBSTEPS) as u64,MECHANICAL_DT)).unwrap();
        Experiment{flexible_sticks:[None,None],mute:None,system:Mechanics::Reference(system),
            force:vec![0.0],stick_weight:1.0,second_stick:None,observer_a:vec![100.0],observer_b:vec![0.0],
            pressure:None,acoustics:Some(boundary),air:None}
    };
    let spec=Spec::parse("frankensim-microphones-v1\nnear_field,0.001,0.000001,14,2000000\nreceiver,0,0,0.025\nreceiver,0.05,0,0.04\npattern,0,0.5,0,0,-1\npattern,1,0,-1,0,0").unwrap();
    let receivers=spec.into_receivers();
    let band=acoustics::stereo::radiation_spec::Spec {band_hz:[40.0,400.0],training_intervals:16,
        max_order:8,subdivisions:0,max_panels:80,max_dense_work:100_000_000};
    let gate=CancelGate::new_clock_free();
    for loaded in [false,true] {
        let (mut actual,mut independent,bake)=if loaded {
            let (a,b)=acoustics::stereo::feedback::prepare(make(),frames,1.0,&receivers,band,&gate).unwrap();
            let (i,_)=acoustics::stereo::feedback::prepare(make(),frames,1.0,&receivers,band,&gate).unwrap();
            (a,i,Some(b))
        }else{(make(),make(),None)};
        actual.system=actual.system.into_analytic_nonlinear().unwrap();
        independent.system=independent.system.into_analytic_nonlinear().unwrap();
        for _ in 0..frames*SUBSTEPS {independent.system.step(&independent.force,&gate).unwrap();}
        let wav=match bake {
            Some(b)=>b.render(&mut actual,frames,1.0,&gate).unwrap(),
            None=>acoustics::stereo::render_receivers_with_spec(&mut actual,frames,1.0,&receivers,band,&gate).unwrap(),
        };
        assert_eq!(&wav[..4],b"RIFF");assert_eq!(&wav[22..24],&2_u16.to_le_bytes());
        assert_eq!(wav.len(),44+4*frames);
        assert!(wav[44..].iter().any(|&v|v!=0),"actual BEM observation must emit nonzero PCM");
        assert!(wav[44..].chunks_exact(4).any(|v|v[..2]!=v[2..]),"independent directional receivers must not duplicate a mono track");
        assert_eq!(actual.system.state(),independent.system.state(),"receivers cannot step mechanics twice or feed back into it");
    }
}
