//! Frontend admission and propagation of the complete reduced material law.
use super::*;

fn options(text: &str) -> Result<Options, String> {
    Options::parse(&text.split_whitespace().map(str::to_owned).collect::<Vec<_>>())
}

#[test]
fn reduction_requires_a_played_geometric_model_and_a_representable_export() {
    assert!(options("").unwrap().board_reduction.is_none());
    let accepted = options("--preset steinway-d --render ritz-piano.wav --board-band-hz 2600 --board-reduction 96,24,800,1600,2400").unwrap();
    let reduction = accepted.board_reduction.unwrap();
    assert_eq!((reduction.max_modes,reduction.keep_low_modes),(96,24));
    assert_eq!(reduction.sample_hz,vec![800.,1600.,2400.]);
    assert_eq!(reduction.transfer_tolerance,None);
    let adaptive = options("--preset steinway-d --render ritz-piano.wav --board-band-hz 2600 --board-reduction adaptive,96,24,0.05,800,1600,2400").unwrap().board_reduction.unwrap();
    assert_eq!((adaptive.max_modes,adaptive.keep_low_modes),(96,24));
    assert_eq!(adaptive.sample_hz,reduction.sample_hz);
    assert_eq!(adaptive.transfer_tolerance,Some(0.05));
    for args in [
        "--render ritz-piano.wav --board-reduction 2,0,100",
        "--preset steinway-d --board-reduction 2,0,100",
        "--preset steinway-d --render ritz-piano.wav --dump-board ritz-board.csv --board-reduction 2,0,100",
        "--preset steinway-d --render ritz-piano.wav --board-band-hz 200 --board-reduction 2,0,300",
        "--preset steinway-d --render ritz-piano.wav --board-reduction 2,0,100 --board-reduction 3,1,200",
        "--preset steinway-d --render ritz-piano.wav --board-reduction 129,0,100",
        "--preset steinway-d --render ritz-piano.wav --board-reduction --board-band-hz 200",
        "--preset steinway-d --render ritz-piano.wav --board-band-hz 200 --board-reduction adaptive,2,0,0.05,300",
        "--preset steinway-d --render ritz-piano.wav --board-reduction adaptive,2,0,0,100",
        "--preset steinway-d --render ritz-piano.wav --board-reduction adaptive,2,0,1,100",
        "--preset steinway-d --render ritz-piano.wav --board-reduction adaptive,2,0,NaN,100",
    ] { assert!(options(args).is_err(),"{args}"); }
    let reduction = board_geometry::ritz::RitzOptions::parse("2,0,100").unwrap();
    let error = prepare_geometric_board_with_reduction(crowned_board::HEADER,&[69],300.,
        true,false,false,0,true,Some(&reduction)).unwrap_err();
    assert!(error.contains("flat geometric board"));
}

#[test]
fn played_frontend_requires_and_applies_the_full_reduced_material_operator() {
    let course = geometry::demonstration_scale().unwrap()[48];
    let modes = board::demonstration();
    let mut selected = Options::default();
    selected.modes = 12;
    selected.dampers = Some("estimated".into());
    selected.board_reduction = Some(board_geometry::ritz::RitzOptions {
        max_modes:modes.len(),keep_low_modes:0,sample_hz:vec![100.],
        selection:board_geometry::ritz::Selection::Pivoted,transfer_tolerance:None,
    });
    assert!(prepare_instrument_with_physical_controls(vec![course],&modes,&selected,None,None)
        .err().unwrap().contains("full projected material damping"));
    let n = modes.len();
    let mut c = vec![2.0;n*n];
    for i in 0..n { c[i*n+i] += 4.0; }
    let mut actual = prepare_instrument_with_board_damping(vec![course],&modes,&selected,
        None,None,Some(&c)).unwrap();
    let mut reference_options = Options::default();
    reference_options.modes = selected.modes;
    reference_options.dampers = selected.dampers.clone();
    let mut reference = prepare_instrument_with_physical_controls(vec![course],&modes,
        &reference_options,None,None).unwrap();
    let mut diagonal = prepare_instrument_with_physical_controls(vec![course],&modes,
        &reference_options,None,None).unwrap();
    reference.configure_bare_board_damping(&c).unwrap();
    for piano in [&mut actual,&mut reference,&mut diagonal] { piano.note_on(69,2.0).unwrap(); }
    for _ in 0..1500 {
        assert_eq!(actual.step().unwrap(),reference.step().unwrap());
        diagonal.step().unwrap();
    }
    assert_eq!(actual.bank.q,reference.bank.q);
    assert_eq!(actual.bank.v,reference.bank.v);
    assert_ne!(actual.bank.v,diagonal.bank.v);
    assert!(actual.accounting.felt_loss_j>0.0 && actual.accounting.modal_loss_j>0.0);
    assert!((actual.accounting.input_work_j-actual.energy_j()-actual.accounting.dissipated_j()).abs()<1e-7);
}

#[test]
fn unframed_wave_force_recovers_the_boards_axial_static_compliance() {
    // One free node has vertical displacement and two independent rocking DOFs.
    // Only key 69 has a longitudinal string and a contact on the horizontal bridge.
    let mut text = String::from("frankensim-board-geometry-si-v1\nsource,estimated,three-DOF axial panel\nsupport,clamped\npretension,0\ndamping,0.01\n");
    for (i, (x, y)) in [(0., 0.5), (0.5, 0.), (1., 0.5), (0.5, 1.), (0.5, 0.5)].iter().enumerate() {
        text.push_str(&format!("node,{i},{x},{y}\n"));
    }
    for (i, t) in [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]].iter().enumerate() {
        text.push_str(&format!("triangle,{i},{},{},{},0.08,450,1000000,80000,0.3,60000,0\n", t[0], t[1], t[2]));
    }
    text.push_str("fixed,0\nfixed,1\nfixed,2\nfixed,3\nbridge,69,0,0,0,1\nbridge,70,0,0,1,0\n");
    text.push_str("stiffener,1000000,80000,0.0001,0.0000000008333333333333333,0.000000001,0.01,450,0,4,2\n");

    let geometry = board_geometry::BoardGeometry::read(&text).unwrap();
    let attachment = geometry.bridge_attachments(&[69]).unwrap()[0];
    assert!(geometry.bridge_attachments(&[70]).is_err());
    let compliance = |board: &board_geometry::PreparedBoard| {
        let row = attachment.axial_row(board.motion.as_ref().unwrap()).unwrap();
        board.modes.iter().zip(row).map(|(mode, a)|
            (a / (std::f64::consts::TAU * mode.frequency_hz)).powi(2)).sum::<f64>()
    };
    let keys = [69, 70];
    let source = geometry.prepare_with_motion(&keys, 1000.).unwrap();
    let expected = compliance(&source);
    assert!(expected > 0.);
    for selection in ["3,0,10,30", "pod,3,0,4,30"] {
        let reduction = board_geometry::ritz::RitzOptions::parse(selection).unwrap();
        let prepare = |wave_keys: &[u8]| prepare_geometric_board_with_source_ports(&text,
            &keys, 1000., false, false, false, 0, true, Some(&reduction), None, wave_keys);
        let vertical_only = prepare(&[]).unwrap();
        assert!(compliance(&vertical_only) < 1e-10 * expected);
        let axial = prepare(&[69]).unwrap();
        assert!((compliance(&axial) - expected).abs() < 1e-9 * expected, "{selection}");
        assert!(axial.modes.len() > vertical_only.modes.len());
        assert!(axial.modes.len() < source.modes.len(), "unused transverse rocking must stay outside the basis");
        assert!(prepare(&[70]).is_err());
    }
}
