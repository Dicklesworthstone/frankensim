//! The harmonic front door must describe the same retained linear piano as
//! playback, including supplied frames and the selected source loss law.
use super::*;

fn args(text: &str) -> Vec<String> { text.split_whitespace().map(str::to_owned).collect() }

fn frame_card() -> String {
    format!("{}\nsource,estimated,harmonic frame integration fixture\ncourse,69,0,0,0,1,0,0,0.02,0,1,0,0,0,1,0.3\n",
        string_polarization::HEADER)
}

fn card_file(text: &str) -> String {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
        .unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-piano-harmonic-{}-{stamp}",std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("frames.fspp");
    std::fs::OpenOptions::new().create_new(true).write(true).open(&path).unwrap()
        .write_all(text.as_bytes()).unwrap();
    path.to_str().unwrap().to_owned()
}

#[test]
fn harmonic_physical_controls_preserve_flag_boundaries_and_explicit_loss_selection() {
    let (options,damping) = admittance_options(&args(
        "--string-polarization frames.fspp --rt0425-string-damping --modes 6 --substeps 8 --edge-cubic-board-mass"
    )).unwrap();
    assert!(damping && options.rt0425_string_damping && options.edge_cubic_board_mass);
    assert_eq!(options.string_polarization.as_deref(),Some("frames.fspp"));
    assert_eq!((options.modes,options.substeps),(6,8));
    assert!(playback::Options::harmonic(&args(
        "--rt0425-string-damping --string-polarization frames.fspp")).is_ok());
    for invalid in [
        "--string-polarization --rt0425-string-damping",
        "--modes --rt0425-string-damping 6",
        "--rt0425-string-damping --rt0425-string-damping",
        "--string-polarization a.fspp --string-polarization b.fspp",
        "--lossless-structure --rt0425-string-damping",
        "--rt0425-string-damping --lossless-structure",
        "--rt0425-string-damping --string-stretching axial.fsps",
        "--string-polarization frames.fspp --dampers estimated",
    ] { assert!(admittance_options(&args(invalid)).is_err(),"{invalid}"); }
    let (_,damping) = admittance_options(&args(
        "--string-polarization frames.fspp --lossless-structure")).unwrap();
    assert!(!damping);
}

#[test]
fn harmonic_source_loss_requires_source_scale_before_loading_files() {
    for command in ["response","admittance"] {
        let mut input = args(&format!("{command} missing.fsb other.csv missing.obj missing.fspe"));
        if command == "admittance" { input.push("69".into()); }
        input.extend(args("unused.csv --rt0425-string-damping"));
        assert!(run(&input).unwrap_err().contains("steinway-d source scale"));
    }
}

#[test]
fn harmonic_frame_admission_precedes_eigensolve_and_refuses_inconsistent_primary_motion() {
    let (board,courses,obj,spec) = tests::small_source_inputs();
    let options = playback::Options {string_polarization:Some(card_file("invalid frame card")),
        ..playback::Options::default()};
    let error = admittance_controlled("invalid board",&courses,&obj,&spec,69,&options,true).unwrap_err();
    assert!(error.contains("polarization"),"{error}");
    let options = playback::Options {string_polarization:Some(card_file(&frame_card()
        .replace("course,69,0,0,0,1,","course,69,0,1,0,0,"))),
        ..playback::Options::default()};
    let error = admittance_controlled(&board,&courses,&obj,&spec,69,&options,true).unwrap_err();
    assert!(error.contains("does not match the existing bridge"),"{error}");
}

#[test]
fn harmonic_and_played_preparation_share_vector_modes_and_loaded_bem_for_both_board_fields() {
    let card = frame_card();
    let path = card_file(&card);
    for cubic in [false,true] {
        let (board,courses,obj,spec) = tests::small_source_inputs();
        let board = board.replace("node,4,0.05,0.05","node,4,0.043,0.054");
        let options = playback::Options {string_polarization:Some(path.clone()),
            rt0425_string_damping:true,modes:6,edge_cubic_board_mass:cubic,
            equilibrate_board_mass:true,..playback::Options::default()};
        let controls = playback::Controls::load(&options,&courses).unwrap();
        let scene = prepare_controlled_body(&board,courses.clone(),Some(&obj),spec,
            &options,controls,false).unwrap();
        let projected = string_polarization::Specification::read(&card,&courses).unwrap()
            .project(&courses,&scene.board.modes,scene.board.motion.as_ref()).unwrap();
        let model = bridge_response::BridgeResponse::new_with_string_damping(&courses,&scene.board.modes,
            RATE*options.substeps as u32,0.45*f64::from(RATE),options.modes,true,
            Some(projected.secondary().0),true).unwrap();
        let played = &scene.piano.bank;
        assert!(played.has_secondary_polarization() && model.bank().has_secondary_polarization());
        assert_eq!(played.modes.len(),model.bank().modes.len());
        assert_eq!(played.strings.len(),model.bank().strings.len());
        for (a,b) in played.strings.iter().zip(&model.bank().strings) {
            assert_eq!((a.course,a.member,a.polarization,a.duplex),(b.course,b.member,b.polarization,b.duplex));
            assert_eq!(a.modes,b.modes);assert_eq!(a.bridge,b.bridge);
        }
        assert!(played.strings.iter().filter(|s|s.polarization==1)
            .any(|s|s.bridge.iter().any(|g|g.abs()>1e-8)));
        for i in 0..played.board_count {
            let mut row=vec![0.;played.board_count];row[i]=1.;
            assert_eq!(played.project_board_shape(&row).unwrap(),model.bank().project_board_shape(&row).unwrap());
        }
        let expected = exterior_loading::sweep(&scene.boundary,&model,&scene.spec,69).unwrap();
        let csv = admittance_controlled(&board,&courses,&obj,&scene.spec,69,&options,true).unwrap();
        assert!(csv.ends_with(&expected));
        assert!(csv.contains("two transverse directions=true"));
        assert!(csv.contains("RT-0425 per-key R_u and eta_u"));
        let estimated = playback::Options {rt0425_string_damping:false,..options.clone()};
        let old = admittance_controlled(&board,&courses,&obj,&scene.spec,69,&estimated,true).unwrap();
        let data = |s:&str|s.lines().filter(|l|!l.starts_with('#')).map(str::to_owned).collect::<Vec<_>>();
        assert_ne!(data(&csv),data(&old),"the source selection must reach physical mobility/power");
        let map = response(&scene).unwrap();
        assert!(map.contains("two transverse directions=true"));
        assert!(map.contains("not force-driven structural response"));
        assert!(played.q.iter().chain(&played.v).all(|x|*x==0.),"harmonic observers must not advance mechanics");
    }
}
