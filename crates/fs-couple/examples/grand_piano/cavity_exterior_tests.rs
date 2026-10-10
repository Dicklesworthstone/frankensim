//! One geometry-derived cavity through the existing BEM/stereo playback path.
use super::*;

fn card() -> &'static str {
    "frankensim-piano-cavity-si-v1\nsource,estimated,authored sealed acoustic integration fixture\ninterface-origin-m,0,0,0\ndimensions-m,0.1,0.1,0.03\nmodes,3\ndamping-ratio,0.03\ngas,dry-air-ussa1976,293.15,101325\n"
}

fn inputs()->(String,Vec<geometry::Course>,String,Specification) {
    let (board,courses,_,_)=tests::small_source_inputs();
    // The small asymmetric physical patch couples the first x/y pressure
    // modes. The exterior is one closed enclosure, with only its top moving.
    let board=board.replace("node,4,0.05,0.05","node,4,0.043,0.054");
    let spec=exterior_geometry::tests::specification()
        .replace("band-hz,40,400,17","band-hz,40,300,41")
        .replace("board-band-hz,400","board-band-hz,300");
    let spec=Specification::read(&format!("{spec}rigid,case\nreceiver-m,0.05,0.05,1\n")).unwrap();
    let obj=exterior_geometry::tests::box_obj("case",[0.,0.,-0.03],[0.1,0.1,0.0315])
        .replace("f -4 -3 -2", "o skin\nf -4 -3 -2")
        .replace("f -8 -7 -3", "o case\nf -8 -7 -3");
    (board,courses,obj,spec)
}

#[test]
fn cavity_reaction_reaches_exterior_pressure_and_shares_the_loaded_stereo_clock() {
    for loaded in [false,true] {
        let options=playback::Options {modes:12,..playback::Options::default()};
        let make=|selected:bool| {
            let (board,courses,obj,spec)=inputs();
            let controls=playback::Controls::from_texts(&courses,None,None,None).unwrap();
            let controls=if selected {controls.with_cavity(card()).unwrap()} else {controls};
            prepare_controlled(&board,courses,&obj,spec,&options,controls).unwrap()
        };
        let mut scene=make(true);let mut manual=make(false);let mut outside_only=make(false);
        let projected=cavity::Specification::read(card()).unwrap().project(&manual.board).unwrap();
        manual.piano.configure_cavity(&projected.loaded(&manual.piano.bank).unwrap()).unwrap();
        assert!(scene.piano.has_cavity());
        assert!(scene.spec.source.contains("authored sealed acoustic integration fixture"));
        assert_eq!(scene.boundary.weights,outside_only.boundary.weights);
        let (baked,_,_)=bake(&mut scene,loaded).unwrap();
        bake(&mut manual,loaded).unwrap();bake(&mut outside_only,loaded).unwrap();
        let frames=2400;
        let score=||performance::Performance::read(
            "sample,event,key,value\n0,note_on,69,0.5\n1000,sustain,0,0.4\n1000,note_off,69,0\n",
            &[69],frames as u64).unwrap();
        let output=exterior_audio::render(&mut scene.piano,score(),frames,&baked,2.).unwrap();
        let baseline=exterior_audio::render(&mut outside_only.piano,score(),frames,&baked,2.).unwrap();
        let mut events=score();
        for sample in 0..frames {
            events.dispatch(sample as u64,&mut manual.piano).unwrap();manual.piano.step().unwrap();
        }
        assert_eq!(scene.piano.bank.q,manual.piano.bank.q);
        assert_eq!(scene.piano.bank.v,manual.piano.bank.v);
        assert_eq!(scene.piano.cavity_energy_j(),manual.piano.cavity_energy_j());
        assert_eq!(scene.piano.radiation_energy_j(),manual.piano.radiation_energy_j());
        assert_eq!(scene.piano.accounting.dissipated_j(),manual.piano.accounting.dissipated_j());
        assert_ne!(scene.piano.bank.q,outside_only.piano.bank.q);
        assert!(scene.piano.cavity_energy_j()>0. && scene.piano.accounting.cavity_loss_j>0.);
        assert!(scene.piano.accounting.felt_loss_j>0.);
        assert_eq!(scene.piano.has_radiation(),loaded);
        if loaded {assert!(scene.piano.accounting.radiation_loss_j>0.);}
        assert!(output.peak_pa>1e-14);
        assert!((output.peak_pa-baseline.peak_pa).abs()>1e-7*baseline.peak_pa,
            "cavity did not materially change modeled pressure: {} vs {}",output.peak_pa,baseline.peak_pa);
        assert!(output.report.contains("Cavity storage"));
        assert!((scene.piano.accounting.input_work_j-scene.piano.energy_j()
            -scene.piano.accounting.dissipated_j()).abs()<1e-7);
        assert_eq!(u16::from_le_bytes([output.wav[22],output.wav[23]]),2);
        let data=output.wav.windows(4).position(|w|w==b"data").unwrap()+8;
        for frame in output.wav[data..].chunks_exact(4) {assert_eq!(&frame[..2],&frame[2..]);}
    }
}

#[test]
fn exterior_cavity_requires_admitted_geometry_and_never_falls_back_after_bad_input() {
    let options=playback::Options::parse(&["--cavity".into(),"enclosure.fspc".into()]).unwrap();
    assert_eq!(options.cavity.as_deref(),Some("enclosure.fspc"));
    for args in [vec!["--cavity"],vec!["--cavity",""],vec!["--cavity","--modes","12"],
        vec!["--cavity","a","--cavity","b"]] {
        assert!(playback::Options::parse(&args.into_iter().map(str::to_owned).collect::<Vec<_>>()).is_err());
    }
    let (board,courses,obj,spec)=inputs();
    let unadmitted=playback::Controls::from_texts(&courses,None,None,None).unwrap();
    assert!(prepare_controlled(&board,courses.clone(),&obj,spec,&options,unadmitted)
        .err().unwrap().contains("cavity controls were not admitted"));
    let bad=playback::Controls::from_texts(&courses,None,None,None).unwrap()
        .with_cavity(&card().replace("interface-origin-m,0,0,0","interface-origin-m,0,0,0.01")).unwrap();
    assert!(prepare_controlled(&board,courses.clone(),&obj,inputs().3,
        &playback::Options::default(),bad).is_err());
    let (bare_board,_,bare_skin,bare_spec)=tests::small_source_inputs();
    let controls=playback::Controls::from_texts(&courses,None,None,None).unwrap().with_cavity(card()).unwrap();
    assert!(prepare_controlled(&bare_board,courses.clone(),&bare_skin,bare_spec,
        &playback::Options::default(),controls).err().unwrap().contains("moving board underside"));
    let controls=playback::Controls::from_texts(&courses,None,None,None).unwrap().with_cavity(card()).unwrap();
    assert!(prepare_controlled_body(&board,courses.clone(),None,inputs().3,
        &playback::Options::default(),controls,false).err().unwrap().contains("outer enclosure OBJ"));
    let without_motion=playback::Controls::from_texts(&courses,None,None,None).unwrap().with_cavity(card()).unwrap();
    assert!(without_motion.instrument(courses,&board::demonstration(),&playback::Options::default())
        .err().unwrap().contains("retained geometric soundboard motion"));
}
