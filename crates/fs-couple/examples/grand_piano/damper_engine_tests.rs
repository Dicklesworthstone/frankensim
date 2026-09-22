//! Spatial damping through the actual hammer/string/board engine, not a PCM envelope.
use super::*;

fn instrument(spatial: bool) -> Instrument {
    let c = super::super::geometry::demonstration_scale().unwrap()[48];
    let mut p = Instrument::new(vec![c], &super::super::board::demonstration(),
        48_000, 4, 24, true).unwrap();
    if spatial {
        p.configure_dampers(&dampers::Specification::estimated(&[c]).unwrap()).unwrap();
    }
    p
}
fn balance(p: &Instrument) {
    assert!((p.accounting.input_work_j-p.energy_j()-p.accounting.dissipated_j()).abs() < 1e-7);
}

#[test]
fn pad_release_changes_the_real_board_trace_not_the_held_attack() {
    let mut point = instrument(false); let mut pad = instrument(true);
    assert_eq!(point.damper_resolution(), None);
    let (strings,cells) = pad.damper_resolution().unwrap();
    assert_eq!(strings, 3); assert!(cells >= 24);
    point.note_on(69, 2.0).unwrap(); pad.note_on(69, 2.0).unwrap();
    let mut a = vec![0.0;point.board_trace_len()]; let mut b = a.clone();
    for _ in 0..1200 {
        point.step_with_board_trace(&mut a).unwrap(); pad.step_with_board_trace(&mut b).unwrap();
        assert_eq!(a,b); assert_eq!(point.bank.q,pad.bank.q); assert_eq!(point.bank.v,pad.bank.v);
    }
    assert!(pad.accounting.felt_loss_j > 0.0);
    point.note_off(69).unwrap(); pad.note_off(69).unwrap();
    point.set_sustain(0.5).unwrap(); pad.set_sustain(0.5).unwrap();
    let mut changed = 0.0_f64;
    for _ in 0..1200 {
        point.step_with_board_trace(&mut a).unwrap(); pad.step_with_board_trace(&mut b).unwrap();
        for (x,y) in a.iter().zip(&b) { changed = changed.max((x-y).abs()); }
    }
    assert!(changed > 1e-15, "pad geometry must reach the same board trace used by pressure rendering");
    assert!(pad.accounting.damper_loss_j > 0.0); balance(&point); balance(&pad);
}

#[test]
fn sostenuto_and_full_sustain_keep_spatial_dampers_lifted_until_release() {
    let mut p = instrument(true); p.note_on(69, 2.0).unwrap();
    for _ in 0..1200 { p.step().unwrap(); }
    p.set_sostenuto(true); p.note_off(69).unwrap();
    for _ in 0..128 { p.step().unwrap(); }
    assert_eq!(p.accounting.damper_loss_j, 0.0);
    p.set_sustain(1.0).unwrap(); p.set_sostenuto(false);
    for _ in 0..128 { p.step().unwrap(); }
    assert_eq!(p.accounting.damper_loss_j, 0.0);
    p.set_sustain(0.5).unwrap();
    for _ in 0..128 { p.step().unwrap(); }
    assert!(p.accounting.damper_loss_j > 0.0); balance(&p);
    // Bad replacement is cold and transactional: preserve the installed pad,
    // motion and material histories rather than falling back to point drag.
    let mut other = p.courses[0]; other.midi = 70;
    let wrong = dampers::Specification::estimated(&[other]).unwrap();
    let resolution = p.damper_resolution(); let energy = p.energy_j();
    let q = p.bank.q.clone(); let v = p.bank.v.clone();
    assert!(p.configure_dampers(&wrong).is_err());
    assert_eq!(p.damper_resolution(),resolution); assert_eq!(p.energy_j(),energy);
    assert_eq!(p.bank.q,q); assert_eq!(p.bank.v,v);
}

#[test]
fn contact_refusal_rolls_back_a_damped_trial_and_retries_exactly() {
    let mut a = instrument(true); let mut b = instrument(true);
    a.note_on(69,2.0).unwrap(); b.note_on(69,2.0).unwrap();
    // Reach a real first contact with kinetic energy, while the hammer is still
    // active. Releasing the key now engages the new damper before contact solve.
    for _ in 0..1200 {
        a.step().unwrap(); b.step().unwrap();
        if a.contacts.iter().any(|c|c.force > 0.0) { break; }
    }
    assert!(a.hammers[0].active && a.contacts.iter().any(|c|c.force > 0.0));
    a.note_off(69).unwrap(); b.note_off(69).unwrap();
    let q = a.bank.q.clone(); let v = a.bank.v.clone(); let energy = a.energy_j();
    let history = a.contacts.clone(); let motion = a.hammers[0].motion;
    let accounting = a.accounting;
    let mut probe = v.clone();
    assert!(a.spatial_dampers.as_ref().unwrap().apply(&mut probe,
        0.5/f64::from(a.bank.rate),0.0, |_|false).unwrap() > 0.0);
    assert_ne!(probe,v, "the refused trial really has a nontrivial damping prefix");
    // Inject a contact-owner refusal AFTER the first damping half-step, not an
    // invalid pedal caught before mechanics. No production failure hook added.
    let diagonal = a.contact_h[0]; a.contact_h[0] = f64::NAN;
    assert!(matches!(a.step(),Err(Error::Contact(_)))); a.contact_h[0] = diagonal;
    assert_eq!(a.bank.q,q); assert_eq!(a.bank.v,v); assert_eq!(a.energy_j(),energy);
    assert_eq!(a.hammers[0].motion,motion);
    assert_eq!(a.accounting.damper_loss_j,accounting.damper_loss_j);
    assert_eq!(a.accounting.input_work_j,accounting.input_work_j);
    assert_eq!(a.accounting.dissipated_j(),accounting.dissipated_j());
    for (old,now) in history.iter().zip(&a.contacts) {
        assert_eq!(old.state.eps_max,now.state.eps_max); assert_eq!(old.memory,now.memory);
        assert_eq!(old.overlap,now.overlap); assert_eq!(old.force,now.force);
    }
    let mut left = vec![0.0;a.board_trace_len()]; let mut right = left.clone();
    for _ in 0..128 {
        a.step_with_board_trace(&mut left).unwrap(); b.step_with_board_trace(&mut right).unwrap();
        assert_eq!(left,right); assert_eq!(a.bank.q,b.bank.q); assert_eq!(a.bank.v,b.bank.v);
        assert_eq!(a.accounting.damper_loss_j,b.accounting.damper_loss_j);
    }
    balance(&a); balance(&b);
}
