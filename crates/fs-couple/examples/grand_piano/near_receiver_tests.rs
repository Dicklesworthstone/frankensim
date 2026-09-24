//! Real near receivers through the original structural, BEM and played paths.
//! Authored small panel fixtures, not factory piano or microphone measurements.
use super::*;
use super::super::{exterior_loading,exterior_audio,performance,playback,geometry};
use fs_math::c64::C64;

fn specification(mode:Option<&str>,close:bool) -> String {
    let mut text=tests::specification()
        .replace("moving,skin","moving,soundboard_skin")
        .replace("band-hz,40,400,17","band-hz,40,300,41")
        .replace("board-band-hz,400","board-band-hz,300");
    let receiver=if close {"receiver-m,0.04,0.04,0.022"} else {"receiver-m,0.05,0.05,1"};
    text=text.replace("receiver-m,0.05,0.05,1",receiver);
    text.push_str(&format!("{receiver}\n"));
    if let Some(mode)=mode {text.push_str(&format!("receiver-evaluation,{mode}\n"));}
    text
}
fn scene(mode:Option<&str>,close:bool) -> super::super::Scene {
    let (board,courses,_,_)=super::super::tests::small_source_inputs();
    let controls=playback::Controls::from_texts(&courses,None,None,None).unwrap();
    super::super::prepare_controlled_body(&board,courses,None,
        Specification::read(&specification(mode,close)).unwrap(),
        &playback::Options::default(),controls,true).unwrap()
}
#[test]
fn receiver_policy_is_explicit_and_legacy_omission_preserves_pressure_and_delays() {
    let original=scene(None,false);let explicit=scene(Some("centroid"),false);
    assert!(!original.spec.near_field_receivers);assert!(!explicit.spec.near_field_receivers);
    let a=original.boundary.sample(&original.spec).unwrap();let b=explicit.boundary.sample(&explicit.spec).unwrap();
    assert_eq!(a.values,b.values);assert_eq!(a.delays_s,b.delays_s);
    for mode in ["", "nearest", "near-field,1e-2"] {
        assert!(Specification::read(&specification(Some(mode),true)).is_err());
    }
    let twice=format!("{}receiver-evaluation,near-field\n",specification(Some("near-field"),true));
    assert!(Specification::read(&twice).is_err());
    let close=scene(None,true);assert!(close.boundary.sample(&close.spec).is_err());
    assert!(exterior_loading::sample(&close.boundary,&close.spec,TAU*100.).is_err());
    let close=scene(Some("near-field"),true);
    let plan=ReceiverSet::for_spec(&close.boundary,&close.spec).unwrap();
    for &delay in plan.delays_s() {assert!((delay-0.0205/340.).abs()<1e-12);}
    let p=close.spec.receivers[0];assert!(norm(sub(p,close.boundary.center))<close.boundary.radius);
    for point in [[0.05,0.05,0.],[0.05,0.05,0.0015],[0.05,0.05,0.005]] {
        assert!(ReceiverSet::new(&close.boundary,&[point],close.spec.medium,true).is_err());
    }
}
#[test]
fn pressure_acceleration_conversion_and_changed_microphones_do_not_change_the_acoustic_load() {
    let close=scene(Some("near-field"),true);let mut spec=Specification::read(&specification(Some("near-field"),false)).unwrap();
    let samples=close.boundary.sample(&close.spec).unwrap();
    for f in [0,20,40] {
        let w=samples.omega[f];
        let near=exterior_loading::sample(&close.boundary,&close.spec,w).unwrap();
        let far=exterior_loading::sample(&close.boundary,&spec,w).unwrap();
        assert_eq!(near.impedance,far.impedance);
        let mut changed=false;
        for (channel,row) in near.receiver_transfer.iter().enumerate() {
            for (input,&h) in row.iter().enumerate() {
                let observed=samples.values[channel][input][f]*C64::new(0.,-w);
                assert!((observed-h).abs()<1e-9*(1.+h.abs()));
                changed|=(h-far.receiver_transfer[channel][input]).abs()>1e-8*h.abs();
            }
        }
        assert!(changed,"actual mic position must change pressure, not impedance");
    }
    // The final assembled scene, including a rigid lid, owns clearance and
    // interior admission. Zero rigid velocity does not make it transparent.
    let asset=tests::box_obj("Lid",[0.02,0.02,0.08],[0.06,0.06,0.004]);
    let manifest="frankensim-piano-rigid-assembly-v1\nsource,estimated,receiver obstacle regression\npart,lid,lid.obj,Lid,1,0,0,0\npose,lid,0,0,0,1,0,0,0,0,0,0\n";
    let rigid=rigid::Assembly::from_text(manifest,|_|Ok(asset.clone())).unwrap();
    let boundary=rigid.attach(close.boundary).unwrap();
    spec.receivers=vec![[0.05,0.05,0.06]];
    let plan=ReceiverSet::for_spec(&boundary,&spec).unwrap();
    assert!((plan.delays_s()[0]-0.02/340.).abs()<1e-12);
    spec.receivers[0]=[0.05,0.05,0.082];
    assert!(ReceiverSet::for_spec(&boundary,&spec).is_err());
}
#[test]
fn actual_close_stereo_playback_observes_one_reacted_trajectory_without_receiver_backreaction() {
    let score=||performance::Performance::read(
        "sample,event,key,value\n0,note_on,69,0.5\n1200,note_off,69,0\n",&[69],2400).unwrap();
    for feedback in [false,true] {
        let mut close=scene(Some("near-field"),true);let mut far=scene(Some("near-field"),false);
        let (baked,samples,_)=super::super::bake(&mut close,feedback).unwrap();
        super::super::bake(&mut far,feedback).unwrap();
        assert!(samples.delays_s.iter().all(|d|*d>0. && *d<0.0001));
        let audio=exterior_audio::render(&mut close.piano,score(),2400,&baked,2.).unwrap();
        let mut manual=score();
        for n in 0..2400 {manual.dispatch(n,&mut far.piano).unwrap();far.piano.step().unwrap();}
        assert_eq!(close.piano.bank.q,far.piano.bank.q);assert_eq!(close.piano.bank.v,far.piano.bank.v);
        assert_eq!(close.piano.radiation_energy_j(),far.piano.radiation_energy_j());
        assert_eq!(close.piano.accounting.radiation_loss_j,far.piano.accounting.radiation_loss_j);
        assert!(close.piano.accounting.felt_loss_j>0.);assert!(audio.peak_pa>1e-14);
        assert!((close.piano.accounting.input_work_j-close.piano.energy_j()-close.piano.accounting.dissipated_j()).abs()<1e-7);
        assert_eq!(u16::from_le_bytes([audio.wav[22],audio.wav[23]]),2);
        let start=audio.wav.windows(4).position(|w|w==b"data").unwrap()+8;
        for frame in audio.wav[start..].chunks_exact(4) {assert_eq!(&frame[..2],&frame[2..]);}
        assert_eq!(close.piano.has_radiation(),feedback);
    }
}
#[test]
fn command_line_uses_supplied_close_positions_and_refuses_interior_output_before_publication() {
    use std::io::Write;
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir=std::env::temp_dir().join(format!("fs-piano-near-{}-{stamp}",std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let (board,courses,_,_)=super::super::tests::small_source_inputs();
    let name=|p:&str|dir.join(p).to_str().unwrap().to_owned();
    for (path,text) in [("board.fsb",board),("strings.csv",geometry::write_scale(&courses)),
        ("near.fspe",specification(Some("near-field"),true)),
        ("inside.fspe",specification(Some("near-field"),true).replace("0.04,0.04,0.022","0.04,0.04,0"))] {
        std::fs::OpenOptions::new().write(true).create_new(true).open(name(path)).unwrap().write_all(text.as_bytes()).unwrap();
    }
    let args=vec!["response".into(),name("board.fsb"),name("strings.csv"),"board-skin-continuous".into(),name("near.fspe"),name("close.csv")];
    super::super::run(&args).unwrap();
    let csv=std::fs::read_to_string(name("close.csv")).unwrap();
    assert!(csv.contains("near-field triangle-integrated"));assert!(!csv.contains("NaN"));
    assert!(csv.lines().filter(|r|!r.starts_with('#')).count()>40);
    assert!(super::super::run(&args).is_err());assert_eq!(std::fs::read_to_string(name("close.csv")).unwrap(),csv);
    let bad=vec!["render-loaded".into(),name("board.fsb"),name("strings.csv"),"board-skin-continuous".into(),
        name("inside.fspe"),name("refused.wav"),"0.05".into()];
    assert!(super::super::run(&bad).is_err());assert!(!std::path::Path::new(&name("refused.wav")).exists());
    // Keep all test-owned files; no user input or output is overwritten.
}
