//! Reduced board motion, dense physical loss and the finite acoustic body must
//! reach both front doors together. This is a small authored integration case,
//! not a mesh-convergence or measured-instrument acceptance claim.
use super::*;

fn args(text:&str)->Vec<String> {text.split_whitespace().map(str::to_owned).collect()}

fn inputs()->(String,Vec<geometry::Course>,String,Specification) {
    let (board,courses,obj,_)=tests::small_source_inputs();
    // Keep the existing mass, thickness, supports and closed acoustic box.
    // Explicitly soften the authored material so all three free plate DOFs
    // fit in the source band; asymmetry couples the bridge to their motion.
    let board=board.replace("node,4,0.05,0.05","node,4,0.043,0.054")
        .replace(",1e7,8e5,0.3,6e5,0",",1e3,8e1,0.3,6e1,0");
    let spec=Specification::read(&exterior_geometry::tests::specification()
        .replace("band-hz,40,400,17","band-hz,40,300,17")
        .replace("board-band-hz,400","board-band-hz,300")).unwrap();
    (board,courses,obj,spec)
}

#[test]
fn reduced_board_dense_loss_and_motion_reach_playback_harmonic_response_and_bem() {
    let (board,courses,obj,spec)=inputs();
    let options=playback::Options::parse(&args(
        "--board-reduction 2,0,80,160 --equilibrate-board-mass --modes 6"
    )).unwrap();
    let source=board_geometry::BoardGeometry::read(&board).unwrap()
        .prepare_with_motion_mass_equilibrated(&[69],spec.board_band_hz).unwrap();
    assert_eq!((source.free_dofs,source.modes.len()),(3,3));
    assert!(source.reduction.is_none());
    let controls=playback::Controls::from_texts(&courses,None,None,None).unwrap();
    let mut scene=prepare_controlled_body(&board,courses.clone(),Some(&obj),spec,
        &options,controls,false).unwrap();
    assert!(scene.board.reduction.is_some());
    assert_eq!(scene.board.modes.len(),2);
    assert_eq!(scene.board.free_dofs,source.free_dofs);
    let c=scene.board.physical_damping.clone().expect("Ritz board needs its full physical damping");
    assert_eq!(c.len(),4);assert!(c.iter().all(|v|v.is_finite()));
    assert!(c[1]!=0.0 && c[1]==c[2],"mixed Ritz coordinates must retain cross damping");
    let motion=scene.board.motion.as_ref().unwrap();
    let build=||playback::Controls::from_texts(&courses,None,None,None).unwrap()
        .instrument_with_motion(courses.clone(),&scene.board.modes,Some(motion),&options).unwrap();
    let mut manual=build();manual.configure_bare_board_damping(&c).unwrap();
    let mut omitted=build();
    let mut model=bridge_response::BridgeResponse::new_with_string_damping(&courses,
        &scene.board.modes,RATE*options.substeps as u32,0.45*f64::from(RATE),
        options.modes,true,None,false).unwrap();
    model.configure_bare_board_damping(&c).unwrap();
    for (a,b) in scene.piano.bank.strings.iter().zip(&model.bank().strings) {
        assert_eq!(a.modes,b.modes);assert_eq!(a.bridge,b.bridge);
    }
    let boundary=Boundary::from_obj(&obj,&scene.spec,motion).unwrap().loaded(model.bank()).unwrap();
    assert_eq!(scene.boundary.weights,boundary.weights);
    assert_eq!(scene.boundary.surface.areas(),boundary.surface.areas());
    assert_eq!(boundary.weights.len(),2);
    let expected=exterior_loading::sweep(&boundary,&model,&scene.spec,69).unwrap();
    let actual=admittance_controlled(&board,&courses,&obj,&scene.spec,69,&options,true).unwrap();
    assert!(actual.ends_with(&expected),"admittance must use the same reduced C and acoustic basis");
    let rows:Vec<_>=actual.lines().filter(|l|!l.starts_with('#')).collect();
    assert_eq!(rows.len(),1+scene.spec.frequencies*(courses.len()+scene.spec.receivers.len()));
    assert!(rows[1..].iter().all(|row|row.split(',').nth(8).unwrap().parse::<f64>().unwrap()>0.0));
    let without=bridge_response::BridgeResponse::new_with_string_damping(&courses,
        &scene.board.modes,RATE*options.substeps as u32,0.45*f64::from(RATE),
        options.modes,true,None,false).unwrap();
    let selected=model.solve(80.0,69,fs_math::c64::C64::ONE,None).unwrap();
    let missing=without.solve(80.0,69,fs_math::c64::C64::ONE,None).unwrap();
    assert!((selected.bridge_velocity[0]-missing.bridge_velocity[0]).abs()
        >1e-12*selected.bridge_velocity[0].abs(),"omitting projected wood loss must change mobility");
    assert!(scene.piano.bank.q.iter().chain(&scene.piano.bank.v).all(|v|*v==0.0));
    for piano in [&mut scene.piano,&mut manual,&mut omitted] {piano.note_on(69,1.5).unwrap();}
    let mut observed=vec![0.0;scene.piano.board_trace_len()];let mut reference=observed.clone();
    let mut changed=false;
    for _ in 0..1200 {
        assert_eq!(scene.piano.step_with_board_trace(&mut observed).unwrap(),
            manual.step_with_board_trace(&mut reference).unwrap());
        assert_eq!(observed,reference);omitted.step().unwrap();
        changed|=scene.piano.bank.v.iter().zip(&omitted.bank.v).any(|(a,b)|(a-b).abs()>1e-15);
    }
    assert!(changed,"played preparation must install the physical damping replacement");
    assert_eq!(scene.piano.bank.q,manual.bank.q);assert_eq!(scene.piano.bank.v,manual.bank.v);
    assert_eq!(scene.piano.accounting.dissipated_j(),manual.accounting.dissipated_j());
    let defect=scene.piano.accounting.input_work_j-scene.piano.energy_j()
        -scene.piano.accounting.dissipated_j();
    assert!(defect.abs()<1e-7);
}

#[test]
fn exterior_reduction_accepts_the_supplied_crown_and_preserves_its_motion() {
    let (flat,_,_,_)=inputs();
    let crown=crowned_board::elevate(&flat,&[0.,0.,0.,0.,0.0015],
        "authored soft shell for reduction regression").unwrap();
    let options=playback::Options::parse(&args("--board-reduction 3,1,80,160")).unwrap();
    let board=prepare_board_motion(&crown,&[69],20000.,&options).unwrap();
    assert!(board.reduction.as_ref().unwrap().source_modes>board.modes.len());
    assert_eq!(board.modes.len(),3);
    assert_eq!(board.physical_damping.as_ref().unwrap().len(),9);
    let motion=board.motion.as_ref().unwrap();
    assert_eq!(motion.mesh.nodes[4][2],0.0015);
    assert_eq!(motion.shapes.len(),board.modes.len());
    let (primary,_)=motion.project_at(0,[0.,0.,1.],[0.;3],[0.,0.,1.]).unwrap();
    for (mode,g) in board.modes.iter().zip(primary) {
        assert!((mode.bridge[48]-g).abs()<1e-12);
    }
}

#[test]
fn exterior_reduction_options_preserve_bounds_flag_values_and_source_band() {
    let parsed=playback::Options::parse(&args(
        "--board-reduction 2,0,80,160 --equilibrate-board-mass"
    )).unwrap();
    let selection=parsed.board_reduction.as_ref().unwrap();
    assert_eq!((selection.max_modes,selection.keep_low_modes),(2,0));
    assert_eq!(selection.sample_hz,vec![80.0,160.0]);
    let (harmonic,damping)=admittance_options(&args(
        "--board-reduction 2,0,80,160 --lossless-structure"
    )).unwrap();
    assert!(!damping);assert!(harmonic.board_reduction.is_some());
    for invalid in [
        "--board-reduction", "--board-reduction --equilibrate-board-mass",
        "--board-reduction 2,0,80 --board-reduction 2,0,160",
        "--board-reduction 0,0,80", "--board-reduction 129,0,80",
        "--board-reduction 2,3,80", "--board-reduction 2,0",
        "--board-reduction 2,0,NaN", "--board-reduction 2,0,0",
        "--board-reduction 2,0,80 --dampers estimated",
    ] {assert!(admittance_options(&args(invalid)).is_err(),"accepted {invalid}");}
    let (board,_,_,spec)=inputs();
    let outside=playback::Options::parse(&args("--board-reduction 2,0,400")).unwrap();
    assert!(prepare_board_motion(&board,&[69],spec.board_band_hz,&outside).is_err());
    let error=prepare_board_motion(crowned_board::HEADER,&[69],spec.board_band_hz,&parsed)
        .unwrap_err();
    assert!(error.contains("flat"),"a crowned header must refuse flat-board mass controls before parsing: {error}");
}
