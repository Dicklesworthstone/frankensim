//! Exercise the public render preparation, not a separate nonlinear oscillator.
use super::*;
use std::{io::Write, sync::atomic::{AtomicU64, Ordering}};

fn options(args: &[&str]) -> Result<Options, String> {
    Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
}
fn input(text: &str) -> String {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("frankensim-string-render-{}-{stamp}-{}.fspx",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(&path).unwrap();
    file.write_all(text.as_bytes()).unwrap();
    path.to_str().unwrap().to_owned()
}
fn selection(row: &str) -> String {
    input(&format!("{}\n{row}\n", linear::string_stretching::HEADER))
}
fn balanced(p: &engine::Instrument) {
    assert!((p.energy_j() + p.accounting.dissipated_j() - p.accounting.input_work_j).abs() < 1e-7);
}

#[test]
fn extension_option_composes_with_existing_playback_and_protects_all_output_paths() {
    assert!(options(&[]).unwrap().string_stretching.is_none());
    let o = options(&["--preset", "steinway-d", "--render", "p.wav", "--midi", "p.mid",
        "--string-stretching", "axial.fspx", "--dampers", "estimated",
        "--hammer-footprints", "faces.fshp", "--microphone-right", "1,1,1"]).unwrap();
    assert_eq!(o.string_stretching.as_deref(), Some("axial.fspx"));
    assert_eq!(o.midi.as_deref(), Some("p.mid"));
    for args in [vec!["--string-stretching"], vec!["--string-stretching", "axial.fspx"],
        vec!["--render", "p.wav", "--string-stretching", ""],
        vec!["--render", "p.wav", "--string-stretching", "  "],
        vec!["--render", "p.wav", "--string-stretching", "--bogus"],
        vec!["--render", "p.wav", "--string-stretching", "a", "--string-stretching", "b"]] {
        assert!(options(&args).is_err(), "accepted {args:?}");
    }
    for output in ["--render", "--dump-scale", "--dump-board", "--dump-geometry", "--dump-obj"] {
        let mut args = vec!["--preset", "steinway-d", "--string-stretching", "axial.fspx"];
        if output != "--render" { args.extend(["--render", "p.wav"]); }
        args.extend([output, "axial.fspx"]);
        assert!(options(&args).is_err(), "source overwritten by {output}");
    }
}

#[test]
fn supplied_extension_is_never_ignored_or_replaced_by_an_estimate() {
    let c = geometry::demonstration_scale().unwrap()[48];
    let mut o = options(&["--render", "p.wav", "--modes", "12"]).unwrap();
    for text in ["not a material file", "frankensim-piano-string-stretching-v1\n",
        "frankensim-piano-string-stretching-v1\nstretch,60,150000,0.2\n",
        "frankensim-piano-string-stretching-v1\nstretch,69,NaN,0.2\n",
        "frankensim-piano-string-stretching-v1\nstretch,69,150000,0.4\n"] {
        o.string_stretching = Some(input(text));
        assert!(prepare_instrument(vec![c], &board::demonstration(), &o).is_err());
    }
    // A path under an existing regular file cannot resolve on any platform.
    o.string_stretching = Some(format!("{}/missing.fspx", input("regular file")));
    assert!(prepare_instrument(vec![c], &board::demonstration(), &o).is_err());
    o.string_stretching = Some(selection("stretch,69,150000,0.2"));
    let all = geometry::demonstration_scale().unwrap();
    assert!(prepare_instrument(vec![c, all[51]], &board::demonstration(), &o).is_err());
}

#[test]
fn renderer_selection_changes_actual_string_tension_without_retuning_the_linear_bank() {
    let c = geometry::demonstration_scale().unwrap()[48];
    let mut o = options(&["--render", "p.wav", "--modes", "12"]).unwrap();
    let mut original = prepare_instrument(vec![c], &board::demonstration(), &o).unwrap();
    o.string_stretching = Some(selection("stretch,69,150000,0.2"));
    let mut selected = prepare_instrument(vec![c], &board::demonstration(), &o).unwrap();
    assert!(selected.bank.has_string_stretching());
    assert!(!original.bank.has_string_stretching());
    assert_eq!(selected.bank.modes.iter().map(|m| m.omega).collect::<Vec<_>>(),
        original.bank.modes.iter().map(|m| m.omega).collect::<Vec<_>>());
    assert_eq!(selected.hammer_contact_count(), original.hammer_contact_count());
    let rest = selected.bank.string_stretching_observation(0).unwrap().tension_n;
    let mut largest_tension = rest;
    let mut largest_energy = 0.0_f64;
    for p in [&mut original, &mut selected] { p.note_on(69, 2.0).unwrap(); }
    for tick in 0..1500 {
        if tick == 600 { for p in [&mut original, &mut selected] { p.note_off(69).unwrap(); } }
        original.step().unwrap(); selected.step().unwrap();
        let o = selected.bank.string_stretching_observation(0).unwrap();
        largest_tension = largest_tension.max(o.tension_n);
        largest_energy = largest_energy.max(o.stretching_energy_j);
        assert!(o.slope_bound <= 0.2);
        balanced(&selected);
    }
    assert!(largest_tension > rest && largest_energy > 0.0);
    assert_ne!(selected.bank.q, original.bank.q);
    assert!(selected.accounting.felt_loss_j > 0.0 && selected.accounting.damper_loss_j > 0.0);
}

#[test]
fn all_linear_file_preserves_default_and_preset_playback_bitwise() {
    let path = selection("linear,69");
    for preset in [false, true] {
        let mut o = options(&["--render", "p.wav", "--modes", "12"]).unwrap();
        if preset { o.preset = Some("steinway-d".to_owned()); }
        let c = selected_scale(None, &o).unwrap()[48];
        let mut original = prepare_instrument(vec![c], &board::demonstration(), &o).unwrap();
        o.string_stretching = Some(path.clone());
        let mut selected = prepare_instrument(vec![c], &board::demonstration(), &o).unwrap();
        assert!(!selected.bank.has_string_stretching());
        for p in [&mut original, &mut selected] { p.note_on(69, 2.0).unwrap(); }
        for tick in 0..1500 {
            if tick == 600 { for p in [&mut original, &mut selected] { p.note_off(69).unwrap(); } }
            assert_eq!(original.step().unwrap().to_bits(), selected.step().unwrap().to_bits());
        }
        assert_eq!(original.bank.q, selected.bank.q);
        assert_eq!(original.bank.v, selected.bank.v);
        assert_eq!(original.accounting.dissipated_j(), selected.accounting.dissipated_j());
        balanced(&selected);
    }
}

#[test]
fn midi_and_csv_reach_the_same_nonlinear_stereo_stream_across_block_boundaries() {
    let c = geometry::demonstration_scale().unwrap()[48];
    let modes = board::demonstration();
    let mut o = options(&["--render", "p.wav", "--modes", "12", "--dampers", "estimated"]).unwrap();
    o.string_stretching = Some(selection("stretch,69,150000,0.2"));
    o.hammer_footprints = Some(input("frankensim-hammer-footprints-v1\nspan,69,0.012,4\n"));
    let a = prepare_instrument(vec![c], &modes, &o).unwrap();
    let b = prepare_instrument(vec![c], &modes, &o).unwrap();
    assert!(a.damper_resolution().is_some());
    assert_eq!(a.hammer_contact_count(), c.unison * 4);
    let bytes = b"MThd\0\0\0\x06\0\0\0\x01\x01\xe0MTrk\0\0\0\x10\
        \0\xb0\x40\x7f\0\x90\x45\x7f\x04\x90\x45\0\x08\xff\x2f\0";
    let midi = midi::read(bytes, &[69], 48_000, 1500,
        midi::Mapping { maximum_velocity_m_s: 2.0, ..midi::Mapping::default() }).unwrap();
    let imported = performance::Performance::from_events(midi.events, &[69], 1500).unwrap();
    let csv = performance::Performance::read("sample,event,key,value\n0,sustain,0,1\n\
        0,note_on,69,2\n200,note_off,69,0\n600,sustain,0,0\n", &[69], 1500).unwrap();
    // Existing manufactured Rayleigh patch fixture, not measured board geometry.
    let surface = [board_geometry::SurfaceSample { position_m: [0.0; 3], area_m2: 0.2,
        mode_shape: modes.iter().map(|m| m.volume / 0.2).collect() }];
    let receivers = [[0.0, 0.0, 1.0], [0.3, 0.0, 1.2]];
    let mut a = audio::AudioStream::new_stereo(a, imported, &surface, receivers, fs_bem::helmholtz::Medium::air()).unwrap();
    let mut b = audio::AudioStream::new_stereo(b, csv, &surface, receivers, fs_bem::helmholtz::Medium::air()).unwrap();
    let mut whole = vec![0.0; 3000]; let mut split = whole.clone();
    a.render_interleaved_block(&mut whole).unwrap();
    for block in split.chunks_mut(274) { b.render_interleaved_block(block).unwrap(); }
    assert_eq!(whole, split);
    assert!(whole.iter().any(|x| x.abs() > 1e-10));
    assert!(whole.chunks_exact(2).any(|f| f[0] != f[1]));
    assert_eq!(a.sample_position(), 1500);
    assert_eq!(a.instrument().bank.q, b.instrument().bank.q);
    assert_eq!(a.instrument().bank.v, b.instrument().bank.v);
    assert!(a.instrument().bank.has_string_stretching());
    balanced(a.instrument());
}
