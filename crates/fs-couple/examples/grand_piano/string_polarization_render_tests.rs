//! Actual geometric preparation and the existing played engine/audio stream.
use super::*;
use std::{fmt::Write as _, io::Write as _, sync::atomic::{AtomicU64, Ordering}};

fn options(args: &[&str]) -> Result<Options, String> {
    Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
}
fn input(text: &str) -> String {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("frankensim-polarization-{}-{stamp}-{}.txt",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(&path).unwrap();
    file.write_all(text.as_bytes()).unwrap();
    path.to_str().unwrap().to_owned()
}
fn panel() -> String {
    let mut text = String::from("frankensim-board-geometry-si-v1\nsource,estimated,asymmetric test panel\nsupport,clamped\npretension,0\ndamping,0.01\n");
    for (i, p) in [[0.,0.],[1.,0.],[1.,1.],[0.,1.],[0.43,0.54]].iter().enumerate() {
        writeln!(text, "node,{i},{},{}", p[0], p[1]).unwrap();
    }
    for (i, t) in [[0,1,4],[1,2,4],[2,3,4],[3,0,4]].iter().enumerate() {
        writeln!(text, "triangle,{i},{},{},{},0.008,450,1e10,8e8,0.3,6e8,0.27", t[0], t[1], t[2]).unwrap();
    }
    text.push_str("fixed,0\nfixed,1\nfixed,2\nfixed,3\nbridge,69,0,0,0,1\n");
    text
}
fn frames() -> String {
    format!("{}\nsource,estimated,explicit test bridge height and frame\ncourse,69,0,0,0,1,0,0,0.02,0,1,0,0,0,1,0.35\n",
        string_polarization::HEADER)
}

#[test]
fn polarization_requires_a_geometric_playback_source_and_cannot_be_silently_ignored() {
    assert!(options(&[]).unwrap().string_polarization.is_none());
    let admitted = options(&["--preset", "steinway-d", "--render", "p.wav",
        "--string-polarization", "frames.fspp", "--midi", "p.mid", "--midi-half-pedal",
        "--consistent-board-mass", "--equilibrate-board-mass"]).unwrap();
    assert_eq!(admitted.string_polarization.as_deref(), Some("frames.fspp"));
    assert!(options(&["--preset", "steinway-d", "--render", "p.wav",
        "--string-polarization", "frames.fspp", "--edge-cubic-board-mass"]).unwrap().edge_cubic_board_mass);
    for args in [vec!["--string-polarization", "frames.fspp"],
        vec!["--render", "p.wav", "--string-polarization", "frames.fspp"],
        vec!["--render", "p.wav", "--board", "modes.csv", "--string-polarization", "frames.fspp"],
        vec!["--preset", "steinway-d", "--render", "p.wav", "--string-polarization", ""],
        vec!["--preset", "steinway-d", "--render", "frames.fspp", "--string-polarization", "frames.fspp"]] {
        assert!(options(&args).is_err(), "accepted {args:?}");
    }
    let c = geometry::demonstration_scale().unwrap()[48];
    assert!(prepare_instrument(vec![c], &board::demonstration(), &admitted)
        .err().unwrap().contains("no one-plane fallback"));
    assert!(prepare_geometric_board_motion(&panel(), &[69], 400., false, false, true, 0, true)
        .unwrap().motion.unwrap().is_edge_cubic());
}

#[test]
fn cubic_bridge_frames_use_the_existing_cubic_field_inside_a_structural_triangle() {
    let text = panel().replace("bridge,69,0,0,0,1", "bridge,69,0,0.2,0.3,0.5");
    let frames = frames().replace("course,69,0,0,0,1", "course,69,0,0.2,0.3,0.5");
    let course = geometry::demonstration_scale().unwrap()[48];
    let spec = string_polarization::Specification::read(&frames, &[course]).unwrap();
    for equilibrated in [false, true] {
        let prepared = prepare_geometric_board_motion(&text, &[69], 400., equilibrated, false, true, 0, true).unwrap();
        let old = prepare_geometric_board(&text, &[69], 400., equilibrated, false, true, 0).unwrap();
        assert_eq!(prepared.modes.len(), old.modes.len());
        for (a, b) in prepared.modes.iter().zip(&old.modes) {
            assert_eq!(a.frequency_hz, b.frequency_hz);
            assert_eq!(a.bridge, b.bridge);
        }
        let motion = prepared.motion.as_ref().unwrap();
        assert!(motion.is_edge_cubic());
        let projected = spec.project(&[course], &prepared.modes, Some(motion)).unwrap();
        assert!(projected.secondary().0[0].iter().any(|g| g.abs() > 1e-8));
        let linear = board_geometry::motion::MotionSurface::new(motion.mesh.clone(), motion.shapes.clone()).unwrap();
        let (p1, _) = linear.project_at(0, [0.2, 0.3, 0.5], [0., 0., 0.02], [0., 0., 1.]).unwrap();
        assert!(p1.iter().zip(&prepared.modes).any(|(a, b)| (a - b.bridge[48]).abs() > 1e-8));
        assert!(spec.project(&[course], &prepared.modes, Some(&linear)).is_err());
    }
}

#[test]
fn supplied_frames_drive_lateral_strings_through_the_same_felt_pedals_and_stereo_clock() {
    for cubic in [false, true] {
    let frame_path = input(&frames());
    let face_path = input("frankensim-hammer-footprints-v1\nspan,69,0.008,2\n");
    let hammer_path = input("frankensim-hammer-materials-v1\nfelt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000\n");
    let mut o = options(&["--preset", "steinway-d", "--board-geometry", "panel.fsb",
        "--render", "p.wav", "--modes", "12", "--string-polarization", &frame_path,
        "--hammers", &hammer_path, "--hammer-footprints", &face_path, "--dampers", "estimated"]).unwrap();
    let c = selected_scale(None, &o).unwrap()[48];
    let board = prepare_geometric_board_motion(&panel(), &[69], 400., true, !cubic, cubic, 0, true).unwrap();
    let spec = string_polarization::Specification::load(&frame_path, &[c]).unwrap();
    let projected = spec.project(&[c], &board.modes, board.motion.as_ref()).unwrap();
    assert!(projected.secondary().0[0].iter().any(|g| g.abs() > 1e-8));
    let stretch = linear::string_stretching::Specification::read(
        "frankensim-piano-string-stretching-v1\nstretch,69,100000,0.2\n", &[c]).unwrap();
    let a = prepare_instrument_with_physical_controls(vec![c], &board.modes, &o, Some(&stretch), Some(&projected)).unwrap();
    let b = prepare_instrument_with_physical_controls(vec![c], &board.modes, &o, Some(&stretch), Some(&projected)).unwrap();
    o.string_polarization = None;
    let old = prepare_instrument_with_physical_controls(vec![c], &board.modes, &o, Some(&stretch), None).unwrap();
    assert_eq!(a.hammer_contact_count(), old.hammer_contact_count());
    assert_eq!(a.bank.strings.len(), 2 * old.bank.strings.len());
    assert_eq!(a.sample_rate(), old.sample_rate());
    assert_eq!(a.bank.rate, old.bank.rate);
    assert_eq!(a.hammer_contact_count(), c.unison * 2);
    let score = || performance::Performance::read("sample,event,key,value\n0,sustain,0,1\n0,note_on,69,2\n600,note_off,69,0\n800,sustain,0,0.5\n1200,sustain,0,0\n", &[69], 1800).unwrap();
    let receivers = [[0.3,0.7,1.2], [0.7,0.2,1.0]];
    let mut a = audio::AudioStream::new_stereo(a, score(), &board.surface, receivers, fs_bem::helmholtz::Medium::air()).unwrap();
    let mut b = audio::AudioStream::new_stereo(b, score(), &board.surface, receivers, fs_bem::helmholtz::Medium::air()).unwrap();
    let mut whole = vec![0.; 3600];
    let mut split = whole.clone();
    a.render_interleaved_block(&mut whole).unwrap();
    for block in split.chunks_mut(254) { b.render_interleaved_block(block).unwrap(); }
    assert_eq!(whole, split);
    assert_eq!(a.sample_position(), 1800);
    assert_eq!(a.instrument().bank.q, b.instrument().bank.q);
    assert_eq!(a.instrument().bank.v, b.instrument().bank.v);
    assert!(whole.iter().all(|p| p.is_finite()));
    assert!(whole.iter().any(|p| p.abs() > 1e-10));
    assert!(whole.chunks_exact(2).any(|f| f[0] != f[1]));
    let p = a.instrument();
    assert!(p.accounting.felt_loss_j > 0. && p.accounting.damper_loss_j > 0.);
    assert!(p.bank.strings.iter().filter(|s| s.polarization == 1).any(|s|
        p.bank.q[s.modes.clone()].iter().any(|q| q.abs() > 1e-14)));
    assert!((p.accounting.input_work_j - p.energy_j() - p.accounting.dissipated_j()).abs() < 1e-7);
    }
}

#[test]
fn the_demonstration_material_front_door_uses_the_same_projected_geometry() {
    let o = Options::default();
    let c = geometry::demonstration_scale().unwrap()[48];
    let board = prepare_geometric_board_motion(&panel(), &[69], 400., false, false, false, 0, true).unwrap();
    let spec = string_polarization::Specification::read(&frames(), &[c]).unwrap();
    let projected = spec.project(&[c], &board.modes, board.motion.as_ref()).unwrap();
    let midplane = string_polarization::Specification::read(
        &frames().replace("0,0,0.02", "0,0,0"), &[c]).unwrap()
        .project(&[c], &board.modes, board.motion.as_ref()).unwrap();
    assert!(midplane.secondary().0[0].iter().all(|g| *g == 0.));
    let mut selected = prepare_instrument_with_physical_controls(vec![c], &board.modes, &o, None, Some(&projected)).unwrap();
    let mut old = prepare_instrument(vec![c], &board.modes, &o).unwrap();
    selected.note_on(69, 1.).unwrap(); old.note_on(69, 1.).unwrap();
    assert_eq!(selected.accounting.input_work_j, old.accounting.input_work_j);
    assert!(selected.bank.has_secondary_polarization());
    assert!(!old.bank.has_secondary_polarization());
    assert_eq!(selected.hammer_contact_count(), old.hammer_contact_count());
}
