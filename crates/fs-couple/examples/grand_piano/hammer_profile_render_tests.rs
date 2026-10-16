//! A supplied face must reach the ordinary material/score/audio preparation.
use super::*;
use std::{io::Write, sync::atomic::{AtomicU64, Ordering}};

fn input(text: &str) -> String {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
        .unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-hammer-profile-render-{}-{stamp}-{}.txt",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::OpenOptions::new().create_new(true).write(true).open(&path).unwrap()
        .write_all(text.as_bytes()).unwrap();
    path.to_str().unwrap().to_owned()
}

#[test]
fn supplied_crown_and_local_felt_reach_one_nonlinear_stereo_timeline() {
    let c = geometry::demonstration_scale().unwrap()[48];
    let face = input(&format!("{}\nprofile,69\n\
        site,69,-0.0015,0.0001,{},0.25\n\
        site,69,0,0,{},0.5\n\
        site,69,0.0015,0.0001,{},0.25\n",
        linear::hammer_footprint::HEADER,
        0.9*c.felt_thickness_m, c.felt_thickness_m, 1.1*c.felt_thickness_m));
    let material = input("frankensim-hammer-materials-v1\n\
        felt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000\n\
        branch,69,2000000,0.0002\n");
    let axial = input(&format!("{}\nstretch,69,150000,0.2\n",
        linear::string_stretching::HEADER));
    let options = Options::parse(&[
        "--render", "profile.wav", "--modes", "12",
        "--hammers", &material, "--hammer-footprints", &face,
        "--string-stretching", &axial, "--dampers", "estimated",
    ].iter().map(|s| (*s).to_owned()).collect::<Vec<_>>()).unwrap();
    let modes = board::demonstration();
    let make = || prepare_instrument(vec![c], &modes, &options).unwrap();
    let score = || performance::Performance::read(
        "sample,event,key,value\n0,sustain,0,1\n0,note_on,69,1.5\n\
         200,note_off,69,0\n600,sustain,0,0\n", &[69], 1500).unwrap();
    // Declared kinematic patch, not an eigensolved or measured soundboard.
    let surface = [board_geometry::SurfaceSample {
        position_m: [0.0; 3], area_m2: 0.2,
        mode_shape: modes.iter().map(|m| m.volume/0.2).collect(),
    }];
    let microphones = [[0.0, 0.0, 1.0], [0.3, 0.0, 1.2]];
    let mut whole = audio::AudioStream::new_stereo(make(), score(), &surface,
        microphones, fs_bem::helmholtz::Medium::air()).unwrap();
    let mut split = audio::AudioStream::new_stereo(make(), score(), &surface,
        microphones, fs_bem::helmholtz::Medium::air()).unwrap();
    let mut manual = make();
    let mut manual_score = score();
    assert_eq!(manual.hammer_contact_count(), 3*c.unison);
    assert!(manual.bank.has_string_stretching() && manual.damper_resolution().is_some());
    assert_eq!(manual.bank.contact_recession_m(0), 0.0001);
    assert_eq!(manual.bank.contact_felt_thickness_m(0), Some(0.9*c.felt_thickness_m));
    let mut a = vec![0.0; 3000];
    let mut b = a.clone();
    whole.render_interleaved_block(&mut a).unwrap();
    for block in b.chunks_mut(74) { split.render_interleaved_block(block).unwrap(); }
    for sample in 0..1500 {
        manual_score.dispatch(sample, &mut manual).unwrap();
        manual.step().unwrap();
    }
    assert_eq!(a, b);
    assert_eq!(whole.sample_position(), 1500);
    assert_eq!(whole.instrument().bank.q, manual.bank.q);
    assert_eq!(whole.instrument().bank.v, manual.bank.v);
    assert_eq!(whole.instrument().accounting.input_work_j, manual.accounting.input_work_j);
    assert_eq!(whole.instrument().accounting.dissipated_j(), manual.accounting.dissipated_j());
    assert!(a.iter().all(|x| x.is_finite()) && a.iter().any(|x| x.abs()>1e-10));
    assert!(a.chunks_exact(2).any(|frame| frame[0]!=frame[1]));
    assert!(manual.accounting.felt_relaxation_loss_j>0.0 && manual.accounting.damper_loss_j>0.0);
    assert!((manual.accounting.input_work_j-manual.energy_j()-manual.accounting.dissipated_j()).abs()<1e-7);
}
