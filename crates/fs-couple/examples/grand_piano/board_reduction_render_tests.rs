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
    for args in [
        "--render ritz-piano.wav --board-reduction 2,0,100",
        "--preset steinway-d --board-reduction 2,0,100",
        "--preset steinway-d --render ritz-piano.wav --dump-board ritz-board.csv --board-reduction 2,0,100",
        "--preset steinway-d --render ritz-piano.wav --board-band-hz 200 --board-reduction 2,0,300",
        "--preset steinway-d --render ritz-piano.wav --board-reduction 2,0,100 --board-reduction 3,1,200",
        "--preset steinway-d --render ritz-piano.wav --board-reduction 129,0,100",
        "--preset steinway-d --render ritz-piano.wav --board-reduction --board-band-hz 200",
    ] { assert!(options(args).is_err(),"{args}"); }
    let reduction = board_geometry::ritz::RitzOptions::parse("2,0,100").unwrap();
    let error = prepare_geometric_board_with_reduction(crowned_board::HEADER,&[69],300.,
        false,false,false,0,true,Some(&reduction)).unwrap_err();
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
