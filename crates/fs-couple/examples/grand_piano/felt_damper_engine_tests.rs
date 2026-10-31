//! G3: a later string-domain refusal rolls back the complete played frame,
//! including already accepted damper travel, contact history and actuator work.
use super::super::{board, geometry, linear::string_stretching};
use super::*;

fn piano(strain_limit: f64) -> Instrument {
    let mut course = geometry::demonstration_scale().unwrap()[48];
    course.unison = 1;
    course.detune_cents = 0.0;
    course.duplex_length_m = 0.0;
    let mut piano =
        Instrument::new(vec![course], &board::demonstration(), 48_000, 8, 12, true).unwrap();
    let strings = string_stretching::Specification::read(
        &format!(
            "{}\nwaves,69,150000,0.2,60000,16,0.02,{strain_limit}\n",
            string_stretching::HEADER
        ),
        &[course],
    )
    .unwrap()
    .with_fixture_axial_rows(&board::demonstration());
    piano.configure_string_stretching(&strings).unwrap();
    let pads = DamperSpecification::read(
        "frankensim-piano-felt-dampers-v1\n\
         pad,69,0.12,0.14,0.00001,0.004,0,0.001,0.001,1,0,2\n\
         felt,69,60000,0.2,2.5,3.2,0.1,0.8,300000\n\
         branch,69,200000,0.002\n",
        &[course],
    )
    .unwrap();
    piano.configure_damper_spec(&pads).unwrap();
    piano
}

fn seed(piano: &mut Instrument) {
    piano.set_sustain(0.5).unwrap();
    piano.bank.v[0] = 0.003;
    // An explicit initial-velocity experiment, not an unreported strike.
    piano.accounting.input_work_j = piano.energy_j();
}

fn strain(piano: &Instrument) -> f64 {
    (piano
        .bank
        .string_stretching_observation(0)
        .unwrap()
        .maximum_tension_n
        - piano.courses[0].tension_n)
        / 150000.0
}

fn pads(piano: &Instrument) -> Vec<felt_dampers::ContactObservation> {
    let owner = piano.felt_dampers().unwrap();
    (0..owner.resolution().1)
        .map(|i| owner.contact_state(i).unwrap())
        .collect()
}

fn axial(piano: &Instrument) -> (Vec<f64>, Vec<f64>) {
    let (q, v) = piano.bank.string_longitudinal_state(0).unwrap();
    (q.to_vec(), v.to_vec())
}

#[test]
fn felt_damper_later_substep_refusal_restores_travel_material_and_axial_state() {
    let mut probe = piano(0.05);
    seed(&mut probe);
    probe.mechanics_step().unwrap();
    probe.mechanics_step().unwrap();
    let second = strain(&probe);
    assert!(probe.felt_dampers().unwrap().pad_height_m(0).unwrap() > 0.0);
    assert!(pads(&probe).iter().any(|p| p.maximum_felt_strain > 0.0));
    probe.mechanics_step().unwrap();
    let third = strain(&probe);
    assert!(third > second && second > 0.0);
    let limit = f64::midpoint(second, third);

    // Establish that this is a refusal AFTER two real accepted mechanics
    // ticks, rather than a preparation failure or an immediate invalid input.
    let mut partial = piano(limit);
    seed(&mut partial);
    partial.mechanics_step().unwrap();
    partial.mechanics_step().unwrap();
    assert!(partial.felt_dampers().unwrap().pad_height_m(0).unwrap() > 0.0);
    assert!(pads(&partial).iter().any(|p| p.maximum_felt_strain > 0.0));
    assert!(matches!(partial.mechanics_step(), Err(Error::Contact(_))));

    let mut actual = piano(limit);
    seed(&mut actual);
    let q = actual.bank.q.clone();
    let v = actual.bank.v.clone();
    let longitudinal = axial(&actual);
    let history = pads(&actual);
    let height = actual.felt_dampers().unwrap().pad_height_m(0);
    let energy = actual.energy_j();
    let accounting = actual.accounting;
    let mut trace = vec![0.0; actual.board_trace_len()];
    for _ in 0..2 {
        assert!(matches!(
            actual.step_with_board_trace(&mut trace),
            Err(Error::Contact(_))
        ));
        assert_eq!(actual.bank.q, q);
        assert_eq!(actual.bank.v, v);
        assert_eq!(axial(&actual), longitudinal);
        assert_eq!(pads(&actual), history);
        assert_eq!(actual.felt_dampers().unwrap().pad_height_m(0), height);
        assert_eq!(actual.energy_j(), energy);
        assert_eq!(actual.damper_energy_j(), 0.0);
        assert_eq!(actual.accounting.input_work_j, accounting.input_work_j);
        assert_eq!(actual.accounting.dissipated_j(), accounting.dissipated_j());
        assert_eq!(
            actual.accounting.damper_actuator_work_j,
            accounting.damper_actuator_work_j
        );
        assert_eq!(
            actual.accounting.damper_relaxation_loss_j,
            accounting.damper_relaxation_loss_j
        );
        assert_eq!(
            actual.accounting.max_balance_error_j,
            accounting.max_balance_error_j
        );
    }
}

fn pedal_piano() -> Instrument {
    let c=geometry::demonstration_scale().unwrap()[48];
    let c=Course{unison:1,detune_cents:0.0,duplex_length_m:0.0,..c};
    let mut p=Instrument::new(vec![c],&board::demonstration(),48_000,4,12,true).unwrap();
    // Declared fixture travel: a 1 mm lift at 0.1 m/s takes 10 ms.
    let pads=DamperSpecification::read(
        "frankensim-piano-felt-dampers-v1\n\
         pad,69,0.12,0.14,0.00001,0.004,0,0.001,0.1,1,0,2\n\
         felt,69,60000,0.2,2.5,3.2,0.1,0.8,300000\n\
         branch,69,200000,0.002\n",&[c]).unwrap();
    p.configure_damper_spec(&pads).unwrap();p
}

#[test]
fn sostenuto_observes_completed_felt_pad_lift_without_recapturing_later_travel() {
    let mut p=pedal_piano();
    assert!(!p.felt_dampers().unwrap().is_fully_lifted(0));
    assert!(!p.felt_dampers().unwrap().is_fully_lifted(1));
    p.set_sustain(1.0).unwrap();p.set_sostenuto(true);
    assert!(!p.hammers[0].latched,"the full-lift command has not moved the pad yet");
    for _ in 0..240 {p.step().unwrap();}
    let half=p.felt_dampers().unwrap().pad_height_m(0).unwrap();
    assert!(half>0.0 && half<0.001);
    p.set_sostenuto(false);p.set_sostenuto(true);
    assert!(!p.hammers[0].latched,"partial travel has no identified capture height");
    for _ in 0..241 {p.step().unwrap();}
    assert!(p.felt_dampers().unwrap().is_fully_lifted(0));
    p.set_sostenuto(true);
    assert!(!p.hammers[0].latched,"the existing capture cannot expand during travel");
    p.set_sostenuto(false);p.set_sostenuto(true);
    assert!(p.hammers[0].latched&&!p.hammers[0].held);
    p.set_sustain(0.0).unwrap();
    for _ in 0..64 {p.step().unwrap();}
    assert_eq!(p.felt_dampers().unwrap().pad_height_m(0),Some(0.001));
    p.set_sostenuto(false);
    for _ in 0..64 {p.step().unwrap();}
    assert!(p.felt_dampers().unwrap().pad_height_m(0).unwrap()<0.001);
    assert_eq!(p.accounting.input_work_j,0.0);
    assert_eq!(p.energy_j(),0.0);
}

#[test]
fn sustained_felt_pad_capture_preserves_ringing_motion_and_release_work() {
    // G3: sustain and sostenuto produce identical accepted geometry and audio
    // while raised, then the same finite-speed return and physical damping.
    let mut captured=pedal_piano();let mut sustained=pedal_piano();
    for p in [&mut captured,&mut sustained] {
        p.set_sustain(1.0).unwrap();
        for _ in 0..600 {p.step().unwrap();}
        p.note_on(69,1.0).unwrap();
        for _ in 0..1200 {p.step().unwrap();}
        p.note_off(69).unwrap();
    }
    let history=pads(&captured);let input=captured.accounting.input_work_j;
    captured.set_sostenuto(true);captured.set_sustain(0.0).unwrap();
    assert!(captured.hammers[0].latched);
    assert_eq!(pads(&captured),history);
    assert_eq!(captured.accounting.input_work_j,input);
    let mut peak=0.0_f64;
    for release in [false,true] {
        if release {captured.set_sostenuto(false);sustained.set_sustain(0.0).unwrap();}
        for _ in 0..768 {
            let actual=captured.step().unwrap();
            assert_eq!(actual,sustained.step().unwrap());peak=peak.max(actual.abs());
            assert_eq!(captured.bank.q,sustained.bank.q);
            assert_eq!(captured.bank.v,sustained.bank.v);
            assert_eq!(pads(&captured),pads(&sustained));
        }
    }
    assert!(peak>1e-10);
    assert!(captured.accounting.damper_loss_j>0.0);
    assert_eq!(captured.accounting.damper_loss_j,sustained.accounting.damper_loss_j);
    assert_eq!(captured.accounting.damper_actuator_work_j,sustained.accounting.damper_actuator_work_j);
    for p in [&captured,&sustained] {
        assert!((p.accounting.input_work_j-p.energy_j()-p.accounting.dissipated_j()).abs()<1e-7);
    }
}
