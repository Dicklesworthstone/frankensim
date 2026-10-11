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
    .unwrap();
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
