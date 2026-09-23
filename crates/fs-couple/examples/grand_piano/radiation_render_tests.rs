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
        "sample,event,key,value\n0,note_on,69,0.1\n1200,note_off,69,0\n",&[69],2400).unwrap();
    let audio=exterior_audio::render(&mut scene.piano,score(),2400,&baked,2.).unwrap();
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
