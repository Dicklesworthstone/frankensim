//! Played history and complete failed-frame rollback for retained axial motion.
use super::super::{board, geometry, linear::string_stretching};
use super::*;

fn piano(strain: f64, loaded: bool) -> Instrument {
    let course = geometry::demonstration_scale().unwrap()[48];
    let faces = hammer_footprint::Specification::read(
        &format!("{}\nspan,69,0.012,4", hammer_footprint::HEADER),
        &[course],
    )
    .unwrap();
    let material = string_stretching::Specification::read(
        &format!(
            "{}\nwaves,69,150000,0.2,60000,16,0.02,{strain}",
            string_stretching::HEADER
        ),
        &[course],
    )
    .unwrap()
    .with_fixture_axial_rows(&board::demonstration());
    let mut p = Instrument::new_with_footprints(
        vec![course],
        &board::demonstration(),
        48000,
        8,
        12,
        true,
        &faces,
    )
    .unwrap();
    p.configure_string_stretching(&material).unwrap();
    if loaded {
        let ports = p.bank.board_count;
        p.configure_radiation(&radiation::Model {
            ports,
            poles: vec![radiation::Pole {
                omega: std::f64::consts::TAU * 420.,
                zeta: 0.1,
                coupling: vec![35.; ports],
            }],
        })
        .unwrap();
    }
    p
}

fn axial(p: &Instrument) -> Vec<(Vec<f64>, Vec<f64>)> {
    p.bank
        .strings
        .iter()
        .enumerate()
        .filter(|(_, s)| s.polarization == 0)
        .filter_map(|(i, _)| p.bank.string_longitudinal_state(i))
        .map(|(q, v)| (q.to_vec(), v.to_vec()))
        .collect()
}

fn same(a: &Instrument, b: &Instrument) {
    assert_eq!(a.bank.q, b.bank.q);
    assert_eq!(a.bank.v, b.bank.v);
    assert_eq!(axial(a), axial(b));
    assert_eq!(a.energy_j(), b.energy_j());
    assert_eq!(a.radiation_energy_j(), b.radiation_energy_j());
    assert_eq!(a.accounting.input_work_j, b.accounting.input_work_j);
    assert_eq!(a.accounting.dissipated_j(), b.accounting.dissipated_j());
    for (a, b) in a.contacts.iter().zip(&b.contacts) {
        assert_eq!(a.overlap, b.overlap);
        assert_eq!(a.force, b.force);
        assert_eq!(a.memory, b.memory);
        assert_eq!(a.state.eps_max, b.state.eps_max);
        assert_eq!(a.state.sig_max, b.state.sig_max);
    }
    for (a, b) in a.hammers.iter().zip(&b.hammers) {
        assert_eq!(a.motion.q, b.motion.q);
        assert_eq!(a.motion.v, b.motion.v);
        assert_eq!(a.jack.elapsed_s, b.jack.elapsed_s);
        assert_eq!(a.held, b.held);
        assert_eq!(a.latched, b.latched);
    }
}

#[test]
fn longitudinal_played_history_and_retry_keep_felt_air_and_pedals_on_one_clock() {
    let mut a = piano(0.05, true);
    let mut b = piano(0.05, true);
    assert!(a.bank.longitudinal_mode_count() > a.bank.strings.len());
    for p in [&mut a, &mut b] {
        p.set_una_corda(true);
        p.note_on(69, 1.).unwrap();
    }
    let mut trace = vec![0.; a.board_trace_len()];
    let mut reference = trace.clone();
    for _ in 0..360 {
        assert_eq!(
            a.step_with_board_trace(&mut trace).unwrap(),
            b.step_with_board_trace(&mut reference).unwrap()
        );
        assert_eq!(trace, reference);
    }
    assert!(a.bank.strings.iter().enumerate().any(|(i, _)| {
        a.bank
            .string_stretching_observation(i)
            .unwrap()
            .longitudinal_kinetic_energy_j
            > 1e-15
    }));
    assert!(a.radiation_energy_j() > 0.);
    let diagonal = a.contact_h[0];
    a.contact_h[0] = f64::NAN;
    assert!(a.step_with_board_trace(&mut trace).is_err());
    a.contact_h[0] = diagonal;
    same(&a, &b);
    for p in [&mut a, &mut b] {
        p.set_sustain(0.3).unwrap();
        p.note_off(69).unwrap();
    }
    for _ in 0..192 {
        assert_eq!(
            a.step_with_board_trace(&mut trace).unwrap(),
            b.step_with_board_trace(&mut reference).unwrap()
        );
        assert_eq!(trace, reference);
    }
    same(&a, &b);
    assert!(
        a.accounting.felt_loss_j > 0.
            && a.accounting.modal_loss_j > 0.
            && a.accounting.damper_loss_j > 0.
    );
    assert!((a.energy_j() + a.accounting.dissipated_j() - a.accounting.input_work_j).abs() < 1e-7);
}

#[test]
fn longitudinal_strain_refusal_after_accepted_substeps_restores_all_axial_history() {
    fn seed(p: &mut Instrument) {
        p.set_sustain(1.).unwrap();
        p.bank.v[0] = 0.001;
        p.accounting.input_work_j = p.energy_j(); // explicit initial kinetic-energy experiment
    }
    fn bound(p: &Instrument) -> f64 {
        let rest = p.courses[0]
            .tension_at_cents(-p.courses[0].detune_cents)
            .unwrap();
        (p.bank
            .string_stretching_observation(0)
            .unwrap()
            .maximum_tension_n
            - rest)
            / 150000.
    }
    let mut probe = piano(0.05, false);
    seed(&mut probe);
    probe.mechanics_step().unwrap();
    probe.mechanics_step().unwrap();
    let second = bound(&probe);
    probe.mechanics_step().unwrap();
    let third = bound(&probe);
    assert!(third > second && second > 0.);
    let limit = f64::midpoint(second, third);
    let mut a = piano(limit, false);
    let mut b = piano(limit, false);
    seed(&mut a);
    seed(&mut b);
    let before_q = a.bank.q.clone();
    let before_v = a.bank.v.clone();
    let before_axial = axial(&a);
    let before_energy = a.energy_j();
    let before_loss = a.accounting.dissipated_j();
    b.mechanics_step().unwrap();
    b.mechanics_step().unwrap();
    assert!(
        axial(&b)
            .iter()
            .any(|(q, v)| q.iter().chain(v).any(|x| *x != 0.))
    );
    assert!(matches!(b.mechanics_step(), Err(Error::Contact(_))));
    for _ in 0..2 {
        assert!(matches!(a.step(), Err(Error::Contact(_))));
        assert_eq!(a.bank.q, before_q);
        assert_eq!(a.bank.v, before_v);
        assert_eq!(axial(&a), before_axial);
        assert_eq!(a.energy_j(), before_energy);
        assert_eq!(a.accounting.dissipated_j(), before_loss);
    }
}
