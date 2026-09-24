use super::*;
use super::super::{tests as fixtures, exterior_geometry, playback, engine, performance,
    exterior_audio, prepare_controlled_body, bake, run};
use std::{io::Write, sync::atomic::{AtomicU64, Ordering}};

fn spec_text() -> String {
    format!("{}receiver-m,0.05,0.05,1\n",exterior_geometry::tests::specification()
        .replace("moving,skin","moving,soundboard_skin")
        .replace("band-hz,40,400,17","band-hz,40,300,41")
        .replace("board-band-hz,400","board-band-hz,300"))
}
fn scene(continuous: bool) -> super::super::Scene {
    let (source,courses,_,_)=fixtures::small_source_inputs();
    let controls=playback::Controls::from_texts(&courses,None,None,Some("estimated")).unwrap();
    prepare_controlled_body(&source,courses,None,Specification::read(&spec_text()).unwrap(),
        &playback::Options::default(),controls,continuous).unwrap()
}
fn balance(p:&engine::Instrument) {
    assert!((p.accounting.input_work_j-p.energy_j()-p.accounting.dissipated_j()).abs()<1e-7);
}
#[test]
fn section_derived_skin_reaches_real_source_hammers_bem_and_stereo_pcm() {
    let score=||performance::Performance::read("sample,event,key,value\n0,note_on,69,0.5\n1200,note_off,69,0\n",&[69],2400).unwrap();
    // One-way and loaded observations use the SAME source geometry and actual
    // mechanical timeline. No independently generated left/right instruments.
    for loaded in [false,true] {
        let mut actual=scene(true);let mut direct=scene(true);
        let (baked,_,_)=bake(&mut actual,loaded).unwrap();bake(&mut direct,loaded).unwrap();
        assert_eq!(actual.boundary.surface.areas().len(),16);
        assert!(actual.spec.source.contains("explicit volume-preserving continuous"));
        let out=exterior_audio::render(&mut actual.piano,score(),2400,&baked,2.).unwrap();
        assert!(out.peak_pa>1e-14);assert!(actual.piano.accounting.felt_loss_j>0.);
        assert_eq!(u16::from_le_bytes([out.wav[22],out.wav[23]]),2);
        let mut schedule=score();
        for n in 0..2400 {schedule.dispatch(n,&mut direct.piano).unwrap();direct.piano.step().unwrap();}
        assert_eq!(actual.piano.bank.q,direct.piano.bank.q);assert_eq!(actual.piano.bank.v,direct.piano.bank.v);
        assert_eq!(actual.piano.radiation_energy_j(),direct.piano.radiation_energy_j());balance(&actual.piano);
        if loaded {assert!(actual.piano.accounting.radiation_loss_j>0.);}
    }
}
#[test]
fn generated_skin_preserves_crowned_preparation_and_source_units_instead_of_flattening() {
    let (source,courses,_,_)=fixtures::small_source_inputs();
    let raised=crowned_board::elevate(&source,&[0.004;5],"explicit 4 mm reference-offset regression").unwrap();
    let controls=playback::Controls::from_texts(&courses,None,None,None).unwrap();
    let actual=prepare_controlled_body(&raised,courses,None,Specification::read(&spec_text()).unwrap(),
        &playback::Options::default(),controls,true).unwrap();
    assert!((actual.boundary.center[2]-0.004).abs()<1e-14);
    let motion=actual.board.motion.as_ref().unwrap();
    let skin=Skin::continuous_from_source(motion,&raised,0.02,2048).unwrap();
    assert!(skin.vertices.iter().all(|p|p[2]>0.));
    assert!(skin.vertices.iter().any(|p|p[2]>0.005));
}
fn file(name:&str,text:&str)->String {
    static NEXT:AtomicU64=AtomicU64::new(0);
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let p=std::env::temp_dir().join(format!("fs-section-skin-{}-{stamp}-{}-{name}",std::process::id(),NEXT.fetch_add(1,Ordering::Relaxed)));
    let mut f=std::fs::OpenOptions::new().write(true).create_new(true).open(&p).unwrap();
    f.write_all(text.as_bytes()).unwrap();p.to_str().unwrap().into()
}
#[test]
fn export_cli_writes_the_complete_skin_without_overwriting_sources_or_falling_back() {
    let (source,courses,_,_)=fixtures::small_source_inputs();
    let board=file("board.fsb",&source);
    let scale=file("scale.csv",&crate::geometry::write_scale(&courses));
    let spec=file("skin.fspe",&spec_text());let output=format!("{board}.obj");
    let args=vec!["export-skin".into(),board.clone(),scale,spec,output.clone(),"--continuous-thickness".into()];
    run(&args).unwrap();let text=std::fs::read_to_string(&output).unwrap();
    let doc=fs_io::obj::read_obj_document(&text).unwrap();assert_eq!(doc.soup.triangles.len(),16);
    assert!(text.contains("explicit volume-preserving continuous"));
    assert!(run(&args).is_err());assert_eq!(std::fs::read_to_string(&output).unwrap(),text);
    assert_eq!(std::fs::read_to_string(&board).unwrap(),source);
    let bogus=file("fake.obj",BODY);let (obj,continuous)=read_body(&bogus).unwrap();
    assert!(obj.is_some() && !continuous);
    let actual=scene(false);
    assert!(boundary(obj.as_deref(),&source,&actual.spec,actual.board.motion.as_ref().unwrap(),false).is_err());
    for name in [BODY,CONTINUOUS_BODY] {assert!(read_body(name).unwrap().0.is_none());}
    assert!(read_body("/missing/section-skin.obj").is_err());
    assert!(run(&["export-skin".into(),board.clone(),"x".into(),"y".into(),"z".into(),"--unknown".into()]).is_err());
}
