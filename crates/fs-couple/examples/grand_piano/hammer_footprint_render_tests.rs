//! The actual CLI constructors, score, shared stereo observers and PCM owner.
use super::*;
fn options(args:&[&str])->Result<Options,String>{
    Options::parse(&args.iter().map(|s|(*s).to_owned()).collect::<Vec<_>>())
}
fn input_file()->String{
    let id=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path=std::env::temp_dir().join(format!("fs-hammer-footprint-{}-{id}.fshp",std::process::id()));
    use std::io::Write;
    let mut file=std::fs::OpenOptions::new().write(true).create_new(true).open(&path).unwrap();
    file.write_all(b"frankensim-hammer-footprints-v1\nspan,69,0.012,4\n").unwrap();
    path.to_str().unwrap().to_owned()
}
#[test]
fn footprint_cli_is_explicit_complete_and_composes_with_source_materials(){
    let o=options(&["--preset","steinway-d","--render","piano.wav","--midi","score.mid",
        "--hammers","felt.fsh","--hammer-footprints","faces.fshp","--dampers","estimated",
        "--microphone-right","1.2,0.8,1"]).unwrap();
    assert_eq!(o.hammer_footprints.as_deref(),Some("faces.fshp"));
    for args in [vec!["--hammer-footprints","faces.fshp"],
        vec!["--render","p.wav","--hammer-footprints",""],
        vec!["--render","faces.fshp","--hammer-footprints","faces.fshp"],
        vec!["--render","p.wav","--hammer-footprints","faces.fshp","--dump-scale","faces.fshp"],
        vec!["--render","p.wav","--hammer-footprints","faces.fshp","--hammer-footprints","other.fshp"]]{
        assert!(options(&args).is_err(),"accepted {args:?}");
    }
    let path=input_file();let mut o=options(&["--preset","steinway-d","--render","p.wav",
        "--hammer-footprints",&path,"--dampers","estimated"]).unwrap();o.modes=12;
    let c=selected_scale(None,&o).unwrap()[48];let modes=board::demonstration();
    let mut piano=prepare_instrument(vec![c],&modes,&o).unwrap();
    assert_eq!(piano.hammer_contact_count(),4*c.unison);assert!(piano.damper_resolution().is_some());
    piano.jack_on(69,70.0,0.007).unwrap();for _ in 0..2400{piano.step().unwrap();}
    assert!(piano.accounting.felt_loss_j>0.0 && piano.accounting.shank_loss_j>0.0);
    // A complete supplied material card composes with the same footprint and shank.
    let material="frankensim-hammer-materials-v1\nfelt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000\nbranch,69,2000000,0.0002\n";
    let mut supplied=prepare_instrument_with_hammers(vec![c],&modes,&o,Some(material)).unwrap();
    assert_eq!(supplied.hammer_contact_count(),4*c.unison);
    supplied.note_on(69,1.0).unwrap();for _ in 0..1200{supplied.step().unwrap();}
    assert!(supplied.accounting.felt_loss_j>0.0);
    let other=geometry::demonstration_scale().unwrap()[47];
    assert!(prepare_instrument(vec![c,other],&modes,&o).is_err(),"unlisted footprint must not silently become point contact");
}

#[test]
fn finite_hammer_midi_and_csv_share_one_stereo_mechanical_timeline(){
    let path=input_file();let mut o=options(&["--render","p.wav","--hammer-footprints",&path,
        "--dampers","estimated"]).unwrap();o.modes=12;
    let c=geometry::demonstration_scale().unwrap()[48];let modes=board::demonstration();
    let bytes=b"MThd\0\0\0\x06\0\0\0\x01\x01\xe0MTrk\0\0\0\x10\
        \0\xb0\x40\x7f\0\x90\x45\x7f\x04\x90\x45\0\x08\xff\x2f\0";
    let decoded=midi::read(bytes,&[69],48_000,1500,
        midi::Mapping{maximum_velocity_m_s:2.0,..midi::Mapping::default()}).unwrap();
    let midi_score=performance::Performance::from_events(decoded.events,&[69],1500).unwrap();
    let csv=||performance::Performance::read("sample,event,key,value\n0,sustain,0,1\n0,note_on,69,2\n200,note_off,69,0\n600,sustain,0,0\n",&[69],1500).unwrap();
    let make=||prepare_instrument(vec![c],&modes,&o).unwrap();
    // Explicit kinematic surface fixture, not a measured/eigensolved panel.
    let surface:Vec<_>=[-0.25,0.25].into_iter().enumerate().map(|(side,x)|board_geometry::SurfaceSample{
        position_m:[x,0.0,0.0],area_m2:0.1,mode_shape:modes.iter().enumerate()
            .map(|(i,b)|b.volume/0.2*if(i+side)%2==0{1.0}else{-0.5}).collect(),
    }).collect();
    let positions=[[-0.4,0.1,0.8],[0.65,0.15,1.2]];
    let mut whole=audio::AudioStream::new_stereo(make(),midi_score,&surface,positions,fs_bem::helmholtz::Medium::air()).unwrap();
    let mut split=audio::AudioStream::new_stereo(make(),csv(),&surface,positions,fs_bem::helmholtz::Medium::air()).unwrap();
    let mut manual=make();let mut manual_score=csv();
    let mut a=vec![0.0;3000];let mut b=a.clone();whole.render_interleaved_block(&mut a).unwrap();
    for block in b.chunks_mut(74){split.render_interleaved_block(block).unwrap();}
    for sample in 0..1500{manual_score.dispatch(sample,&mut manual).unwrap();manual.step().unwrap();}
    assert_eq!(a,b);assert_eq!(whole.sample_position(),1500);
    assert_eq!(whole.instrument().bank.q,manual.bank.q);assert_eq!(whole.instrument().bank.v,manual.bank.v);
    assert_eq!(whole.instrument().accounting.input_work_j,manual.accounting.input_work_j);
    assert!(a.iter().all(|p|p.is_finite()) && a.iter().any(|p|p.abs()>1e-10));
    assert!(a.chunks_exact(2).any(|p|(p[0]-p[1]).abs()>1e-10));
    let(wav,_)=fs_couple::pcm_wav::encode_pcm16_wav_interleaved(&a,48_000,2,2.0).unwrap();
    assert_eq!(wav.len(),6044);assert_eq!(&wav[22..24],&2_u16.to_le_bytes());
    assert!(manual.accounting.felt_loss_j>0.0 && manual.accounting.damper_loss_j>0.0);
    assert!((manual.accounting.input_work_j-manual.energy_j()-manual.accounting.dissipated_j()).abs()<1e-7);
}
