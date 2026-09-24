//! Actual assembly files, native skin, common BEM and the existing piano host.
//! Synthetic geometry fixtures are not measured Steinway cabinet dimensions.
use super::*;
use std::{path::PathBuf, sync::atomic::{AtomicU64, Ordering}};

fn directory() -> PathBuf {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-piano-rigid-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&dir).unwrap();
    std::fs::create_dir(dir.join("parts")).unwrap();
    dir
}
fn save(path: &std::path::Path, text: &str) {
    let mut file = std::fs::OpenOptions::new().create_new(true).write(true).open(path).unwrap();
    file.write_all(text.as_bytes()).unwrap();
}
fn assembly_file() -> (PathBuf, String) {
    let dir = directory();
    let obj = exterior_geometry::tests::box_obj("Lid", [20.,20.,100.], [50.,50.,4.]);
    save(&dir.join("parts/lid.obj"), &obj);
    let manifest = format!("{}\nsource,estimated,static cabinet composition regression\npart,lid,parts/lid.obj,Lid,0.001,0,0,0\npose,lid,0.02,0.02,0.1,1,0,0,30,0,0,0\n",
        exterior_geometry::rigid::HEADER);
    save(&dir.join("scene.fspr"), &manifest);
    let name = dir.join("scene.fspr").to_str().unwrap().to_owned();
    (dir, name)
}
fn acoustic_text() -> String {
    format!("{}receiver-m,0.05,0.05,1\n", exterior_geometry::tests::specification()
        .replace("moving,skin", "moving,soundboard_skin")
        .replace("band-hz,40,400,17", "band-hz,40,300,41")
        .replace("board-band-hz,400", "board-band-hz,300"))
}
fn scene(path: Option<&str>) -> Scene {
    let (board, courses, _, _) = tests::small_source_inputs();
    let options = playback::Options { rigid_assembly: path.map(str::to_owned), ..playback::Options::default() };
    let controls = playback::Controls::from_texts(&courses, None, None, None).unwrap();
    prepare_controlled_body(&board, courses, None, Specification::read(&acoustic_text()).unwrap(),
        &options, controls, true).unwrap()
}
fn args(values: &[&str]) -> Vec<String> { values.iter().map(|s| (*s).to_owned()).collect() }

#[test]
fn rigid_selection_is_explicit_and_bad_sources_refuse_before_structural_preparation() {
    let good = args(&["--rigid-assembly", "scene.fspr", "--modes", "128", "--substeps", "8"]);
    assert_eq!(playback::Options::parse(&good).unwrap().rigid_assembly.as_deref(), Some("scene.fspr"));
    assert!(playback::Options::harmonic(&good).is_ok());
    let mut lossless = good.clone(); lossless.push("--lossless-structure".into());
    assert!(!admittance_options(&lossless).unwrap().1);
    for bad in [vec!["--rigid-assembly"], vec!["--rigid-assembly", " "],
        vec!["--rigid-assembly", "a", "--rigid-assembly", "b"],
        vec!["--modes", "--rigid-assembly", "a", "24"],
        vec!["--rigid-assembly", "--modes", "24"]] {
        assert!(playback::Options::parse(&args(&bad)).is_err());
        assert!(admittance_options(&args(&bad)).is_err());
    }
    let (dir, path) = assembly_file();
    let mut options = playback::Options { rigid_assembly: Some(path), ..playback::Options::default() };
    options.rigid_assembly = Some(dir.join("absent.fspr").to_str().unwrap().to_owned());
    let (_, courses, _, _) = tests::small_source_inputs();
    let controls = playback::Controls::from_texts(&courses, None, None, None).unwrap();
    let failed = prepare_controlled_body("invalid structure", courses, None,
        Specification::read(&acoustic_text()).unwrap(), &options, controls, true);
    let error = match failed { Err(error) => error, Ok(_) => panic!("missing assembly selected a substitute") };
    assert!(error.contains("absent.fspr"), "{error}");
}

#[test]
fn relative_asset_export_retains_admitted_geometry_and_never_overwrites_input() {
    let (dir, path) = assembly_file(); let output = dir.join("posed.obj");
    let name = output.to_str().unwrap();
    let admitted = Assembly::load(&path).unwrap();
    let source = std::fs::read_to_string(dir.join("parts/lid.obj")).unwrap();
    run(&args(&["export-rigid", &path, name])).unwrap();
    let exported = std::fs::read_to_string(&output).unwrap();
    assert_eq!(exported, admitted.obj());
    assert_eq!(fs_io::obj::read_obj_document(&exported).unwrap().soup.triangles.len(), 12);
    assert!(run(&args(&["export-rigid", &path, name])).is_err());
    assert_eq!(std::fs::read_to_string(&output).unwrap(), exported);
    assert_eq!(std::fs::read_to_string(dir.join("parts/lid.obj")).unwrap(), source);
    // The admitted object holds geometry, not a filename reopened during prep.
    // This is a mutation of the fresh test-owned input only, never user assets.
    std::fs::write(dir.join("parts/lid.obj"), "deliberately invalid later input").unwrap();
    assert_eq!(admitted.obj(), exported); assert!(Assembly::load(&path).is_err());
}

#[test]
fn imported_lid_and_native_board_skin_reach_reacted_stereo_on_one_mechanical_clock() {
    let (_, path) = assembly_file();
    for feedback in [false, true] {
        let mut actual = scene(Some(&path)); let mut manual = scene(Some(&path)); let mut bare = scene(None);
        let base_count = bare.boundary.surface.areas().len();
        assert_eq!(actual.boundary.surface.areas().len(), base_count + 12);
        assert_eq!(actual.board.mass_kg, bare.board.mass_kg);
        assert_eq!(actual.board.frequency_intervals_hz, bare.board.frequency_intervals_hz);
        for (with_lid, without) in actual.boundary.weights.iter().zip(&bare.boundary.weights) {
            assert_eq!(&with_lid[..base_count], without); assert!(with_lid[base_count..].iter().all(|x| *x == 0.));
        }
        assert_eq!(&actual.boundary.surface.triangles().unwrap()[..base_count], bare.boundary.surface.triangles().unwrap());
        assert!(actual.spec.source.contains("Rigid assembly"));
        let w = std::f64::consts::TAU * 200.;
        let reacted = exterior_loading::sample(&actual.boundary, &actual.spec, w).unwrap();
        let reference = exterior_loading::sample(&bare.boundary, &bare.spec, w).unwrap();
        let difference: f64 = reacted.impedance.iter().zip(&reference.impedance).map(|(a,b)| (*a-*b).abs()).sum();
        let magnitude: f64 = reference.impedance.iter().map(|z| z.abs()).sum();
        assert!(difference > 1e-9 * magnitude, "rigid geometry must change physical fluid loading");
        let (baked, _, _) = bake(&mut actual, feedback).unwrap();
        bake(&mut manual, feedback).unwrap();
        let score = || performance::Performance::read(
            "sample,event,key,value\n0,note_on,69,0.5\n1200,note_off,69,0\n", &[69], 2400).unwrap();
        let audio = exterior_audio::render(&mut actual.piano, score(), 2400, &baked, 2.).unwrap();
        let mut schedule = score(); let mut bare_schedule = score();
        for n in 0..2400 {
            schedule.dispatch(n, &mut manual.piano).unwrap(); manual.piano.step().unwrap();
            bare_schedule.dispatch(n, &mut bare.piano).unwrap(); bare.piano.step().unwrap();
        }
        assert_eq!(actual.piano.bank.q, manual.piano.bank.q); assert_eq!(actual.piano.bank.v, manual.piano.bank.v);
        assert_eq!(actual.piano.radiation_energy_j(), manual.piano.radiation_energy_j());
        if !feedback { assert_eq!(actual.piano.bank.q, bare.piano.bank.q); }
        assert!(actual.piano.accounting.felt_loss_j > 0.); assert!(audio.peak_pa > 1e-14);
        assert!((actual.piano.accounting.input_work_j - actual.piano.energy_j() - actual.piano.accounting.dissipated_j()).abs() < 1e-7);
        assert_eq!(u16::from_le_bytes([audio.wav[22], audio.wav[23]]), 2);
        let start = audio.wav.windows(4).position(|w| w == b"data").unwrap() + 8;
        for frame in audio.wav[start..].chunks_exact(4) { assert_eq!(&frame[..2], &frame[2..]); }
    }
}

#[test]
fn response_and_admittance_cli_use_the_same_loaded_native_skin_and_imported_parts() {
    let (dir, path) = assembly_file(); let (board, courses, _, _) = tests::small_source_inputs();
    save(&dir.join("board.fsb"), &board); save(&dir.join("strings.csv"), &geometry::write_scale(&courses));
    save(&dir.join("air.fspe"), &acoustic_text());
    let p = |name: &str| dir.join(name).to_str().unwrap().to_owned();
    run(&["response".into(),p("board.fsb"),p("strings.csv"),"board-skin-continuous".into(),
        p("air.fspe"),p("response.csv"),"--rigid-assembly".into(),path.clone()]).unwrap();
    run(&["admittance".into(),p("board.fsb"),p("strings.csv"),"board-skin-continuous".into(),
        p("air.fspe"),"69".into(),p("bridge.csv"),"--rigid-assembly".into(),path]).unwrap();
    for name in ["response.csv", "bridge.csv"] {
        let csv = std::fs::read_to_string(p(name)).unwrap();
        assert!(csv.contains("Rigid assembly")); assert!(csv.contains("30 deg"));
        assert!(!csv.contains("NaN") && !csv.contains("inf,"));
        assert!(csv.lines().filter(|r| !r.starts_with('#')).count() > 40);
    }
}
