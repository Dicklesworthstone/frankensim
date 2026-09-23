//! The real source-hammer/board/BEM/fit/feedback/stereo consumer chain.
use super::*;
use fs_math::c64::C64;

#[test]
fn loaded_bem_playback_changes_mechanics_and_preserves_one_stereo_clock() {
    let mut scene=tests::small_source_scene();
    let mut manual=tests::small_source_scene();let mut bare=tests::small_source_scene();
    let (baked,samples,report)=bake(&mut scene,true).unwrap();
    bake(&mut manual,true).unwrap();
    assert!(scene.piano.has_radiation());assert!(report.contains("Passive load:"));
    // The same velocity solve supplied the acceleration observer, with the
    // negative-time i/omega conversion, not a second independently fitted body.
    let f=samples.omega.len()/2;let w=samples.omega[f];
    let direct=exterior_loading::sample(&scene.boundary,&scene.spec,w).unwrap();
    for (channel,row) in direct.receiver_transfer.iter().enumerate() {
        for (input,&h) in row.iter().enumerate() {
            let restored=samples.values[channel][input][f]*C64::new(0.,-w);
            assert!((restored-h).abs()<1e-10*(1.+h.abs()));
        }
    }
    let score=||performance::Performance::read(
        "sample,event,key,value\n0,note_on,69,0.5\n1200,note_off,69,0\n",&[69],2400).unwrap();
    let audio=exterior_audio::render(&mut scene.piano,score(),2400,&baked,2.).unwrap();
    assert!(scene.piano.accounting.felt_loss_j>0.,"the source hammer must actually strike");
    assert!(audio.peak_pa>1e-14);assert!(audio.report.contains("Passive acoustic feedback"));
    assert!(audio.report.contains("Acoustic storage"));assert!(!audio.report.contains("no radiation backreaction"));
    assert_eq!(u16::from_le_bytes([audio.wav[22],audio.wav[23]]),2);
    let mut program=score();let mut unreacted=score();
    for n in 0..2400 {
        program.dispatch(n,&mut manual.piano).unwrap();manual.piano.step().unwrap();
        unreacted.dispatch(n,&mut bare.piano).unwrap();bare.piano.step().unwrap();
    }
    assert_eq!(scene.piano.bank.q,manual.piano.bank.q);assert_eq!(scene.piano.bank.v,manual.piano.bank.v);
    assert_eq!(scene.piano.radiation_energy_j(),manual.piano.radiation_energy_j());
    assert_ne!(scene.piano.bank.q,bare.piano.bank.q);
    assert!(scene.piano.accounting.radiation_loss_j>0.);
    let closure=scene.piano.accounting.input_work_j-scene.piano.energy_j()-scene.piano.accounting.dissipated_j();
    assert!(closure.abs()<1e-7,"{closure:e}");
    let data=audio.wav.windows(4).position(|w|w==b"data").unwrap()+8;
    for frame in audio.wav[data..].chunks_exact(4) {assert_eq!(&frame[..2],&frame[2..]);}
}

#[test]
fn failed_loaded_admission_never_attaches_a_substitute_or_advances_the_piano() {
    let mut scene=tests::small_source_scene();scene.spec.frequencies=17;
    let q=scene.piano.bank.q.clone();let energy=scene.piano.energy_j();
    assert!(bake(&mut scene,true).is_err());assert!(!scene.piano.has_radiation());
    assert_eq!(scene.piano.bank.q,q);assert_eq!(scene.piano.energy_j(),energy);
    for duration in ["NaN","-1","0"] {
        let args=["render-loaded","missing.fss","missing.csv","missing.obj","missing.fspe","unused.wav",duration].map(str::to_owned);
        assert!(run(&args).is_err());
    }
    assert!(USAGE.contains("render-loaded"));
}

#[test]
fn full_physical_controls_reach_loaded_stereo_without_a_second_mechanical_clock() {
    let options=playback::Options {substeps:8,modes:48,..playback::Options::default()};
    let make=|| {
        let (board,courses,obj,spec)=tests::small_source_inputs();
        let materials="frankensim-hammer-materials-v1\nfelt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000\nbranch,69,2000000,0.0002\n";
        let faces="frankensim-hammer-footprints-v1\nspan,69,0.008,2\n";
        let controls=playback::Controls::from_texts(&courses,Some(materials),Some(faces),Some("estimated")).unwrap();
        prepare_controlled(&board,courses,&obj,spec,&options,controls).unwrap()
    };
    let mut scene=make();let mut manual=make();
    assert_eq!(scene.piano.bank.rate,RATE*8);
    assert_eq!(scene.piano.bank.contact_strings.len(),6); // two sites on three actual strings
    assert!(scene.piano.damper_resolution().is_some());
    let (baked,_,_)=bake(&mut scene,true).unwrap();bake(&mut manual,true).unwrap();
    let score=||performance::Performance::read(
        "sample,event,key,value\n0,note_on,69,0.5\n1200,note_off,69,0\n",&[69],2400).unwrap();
    let audio=exterior_audio::render(&mut scene.piano,score(),2400,&baked,2.).unwrap();
    let mut events=score();
    for n in 0..2400 {events.dispatch(n,&mut manual.piano).unwrap();manual.piano.step().unwrap();}
    assert_eq!(scene.piano.bank.q,manual.piano.bank.q);assert_eq!(scene.piano.bank.v,manual.piano.bank.v);
    assert_eq!(scene.piano.radiation_energy_j(),manual.piano.radiation_energy_j());
    assert!(scene.piano.accounting.felt_loss_j>0.);assert!(scene.piano.accounting.damper_loss_j>0.);
    assert!(scene.piano.accounting.radiation_loss_j>0.);assert!(audio.peak_pa>1e-14);
    assert!((scene.piano.accounting.input_work_j-scene.piano.energy_j()-scene.piano.accounting.dissipated_j()).abs()<1e-7);
    let data=audio.wav.windows(4).position(|w|w==b"data").unwrap()+8;
    for frame in audio.wav[data..].chunks_exact(4) {assert_eq!(&frame[..2],&frame[2..]);}
}

#[test]
fn force_driven_csv_reaches_loaded_stereo_and_matches_direct_owner_controls() {
    let mut scene=tests::small_source_scene();let mut manual=tests::small_source_scene();
    let (baked,_,_)=bake(&mut scene,true).unwrap();bake(&mut manual,true).unwrap();
    let text="sample,event,key,value\n0,sustain,0,0.5\n0,una_corda,0,1\n24,jack_staccato,69,70\n240,sostenuto,0,1\n1200,note_off,69,0\n1440,sostenuto,0,0\n1680,sustain,0,0\n2000,una_corda,0,0\n";
    let score=playback::Score::csv(text,&[69],4800).unwrap();
    assert!(score.report.contains("peak N"));
    let audio=exterior_audio::render(&mut scene.piano,score.performance,4800,&baked,2.).unwrap();
    for n in 0..4800 {
        match n {
            0=>{manual.piano.set_sustain(0.5).unwrap();manual.piano.set_una_corda(true);},
            24=>manual.piano.jack_on(69,70.,0.007).unwrap(),
            240=>manual.piano.set_sostenuto(true),
            1200=>manual.piano.note_off(69).unwrap(),
            1440=>manual.piano.set_sostenuto(false),
            1680=>manual.piano.set_sustain(0.).unwrap(),
            2000=>manual.piano.set_una_corda(false),
            _=>{},
        }
        manual.piano.step().unwrap();
    }
    assert_eq!(scene.piano.bank.q,manual.piano.bank.q);assert_eq!(scene.piano.bank.v,manual.piano.bank.v);
    assert_eq!(scene.piano.radiation_energy_j(),manual.piano.radiation_energy_j());
    assert_eq!(scene.piano.accounting.input_work_j,manual.piano.accounting.input_work_j);
    assert!(scene.piano.accounting.input_work_j>0.);assert!(scene.piano.accounting.felt_loss_j>0.);
    assert!(scene.piano.accounting.shank_loss_j>0.);assert!(scene.piano.accounting.radiation_loss_j>0.);
    assert!(audio.peak_pa>1e-14);assert!(audio.report.contains("Passive acoustic feedback"));
    assert!((scene.piano.accounting.input_work_j-scene.piano.energy_j()-scene.piano.accounting.dissipated_j()).abs()<1e-7);
    let data=audio.wav.windows(4).position(|w|w==b"data").unwrap()+8;
    for frame in audio.wav[data..].chunks_exact(4) {assert_eq!(&frame[..2],&frame[2..]);}
}
