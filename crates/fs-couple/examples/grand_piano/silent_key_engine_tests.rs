//! Silent key controls on the real hammer/string/shared-board instrument.
use super::*;

fn instrument(spatial: bool, receiver_connected: bool) -> Instrument {
    // Two equally tuned physical strings with distinct key addresses isolate
    // sympathetic transfer; MIDI numbers do not supply oscillator frequencies.
    let c = super::super::geometry::demonstration_scale().unwrap()[39];
    let c = Course { unison: 1, duplex_length_m: 0.0, detune_cents: 0.0, ..c };
    let courses = vec![c, Course { midi: 64, ..c }];
    let mut bridge = [0.0; 88];
    bridge[39] = 0.7;
    bridge[43] = if receiver_connected { 0.6 } else { 0.0 };
    let mut p = Instrument::new(courses.clone(), &[BoardMode {
        frequency_hz: 180.0, damping_ratio: 0.0, bridge, volume: 0.1,
    }], 48_000, 4, 12, false).unwrap();
    if spatial {
        p.configure_dampers(&dampers::Specification::estimated(&courses).unwrap()).unwrap();
    }
    p
}

fn assert_only_key_hold_changes(p: &mut Instrument, key: u8) {
    let ci = p.key_index(key).unwrap();
    let old = p.hammers[ci];
    let history = p.contacts.clone();
    let q = p.bank.q.clone();
    let v = p.bank.v.clone();
    let energy = p.energy_j();
    let accounting = p.accounting;
    p.silent_key_down(key).unwrap();
    let now = p.hammers[ci];
    assert!(now.held);
    assert_eq!(now.motion, old.motion);
    assert_eq!((now.active, now.latched, now.on_rest), (old.active, old.latched, old.on_rest));
    assert_eq!((now.jack.peak_n, now.jack.duration_s, now.jack.elapsed_s),
        (old.jack.peak_n, old.jack.duration_s, old.jack.elapsed_s));
    assert_eq!(p.bank.q, q);
    assert_eq!(p.bank.v, v);
    assert_eq!(p.energy_j(), energy);
    assert_eq!(p.accounting.input_work_j, accounting.input_work_j);
    assert_eq!(p.accounting.dissipated_j(), accounting.dissipated_j());
    assert_eq!(p.accounting.max_balance_error_j, accounting.max_balance_error_j);
    for (before, after) in history.iter().zip(&p.contacts) {
        assert_eq!(before.state, after.state);
        assert_eq!(before.memory, after.memory);
        assert_eq!((before.overlap, before.force, before.enabled),
            (after.overlap, after.force, after.enabled));
    }
}

#[test]
fn silent_key_hold_adds_no_work_and_preserves_rearming_and_contact_history() {
    let mut p = instrument(false, true);
    assert!(matches!(p.silent_key_down(61), Err(Error::UnknownKey)));
    assert!(p.hammers.iter().all(|h| !h.held));
    assert_only_key_hold_changes(&mut p, 60);
    assert_only_key_hold_changes(&mut p, 60);
    assert!(matches!(p.note_on(60, 2.0), Err(Error::NotRearmed)));
    for _ in 0..128 { assert_eq!(p.step().unwrap(), 0.0); }
    assert_eq!(p.energy_j(), 0.0);
    assert_eq!(p.accounting.input_work_j, 0.0);
    assert_eq!(p.accounting.dissipated_j(), 0.0);
    assert!(p.hammers.iter().all(|h| !h.active && h.on_rest));

    p.note_off(60).unwrap();
    p.note_on(60, 3.0).unwrap();
    for _ in 0..1200 {
        p.step().unwrap();
        if p.contacts.iter().any(|c| c.force > 0.0) { break; }
    }
    assert!(p.hammers[0].active);
    assert!(p.contacts.iter().any(|c| c.state.eps_max > 0.0));
    assert!(p.contacts.iter().any(|c| c.memory.0.iter().any(|v| *v != 0.0)));
    p.note_off(60).unwrap();
    // Catch the real, moving instrument silently, preserving nonzero felt
    // history and an airborne hammer instead of synthesizing a quiet strike.
    assert_only_key_hold_changes(&mut p, 60);
    assert_only_key_hold_changes(&mut p, 60);
    assert!(matches!(p.note_on(60, 2.0), Err(Error::NotRearmed)));
    p.note_off(60).unwrap();
    assert!(!p.hammers[0].held);
}

#[test]
fn silent_sostenuto_receives_real_shared_board_energy_then_half_pedal_damps_it() {
    for spatial in [false, true] {
        let mut p = instrument(spatial, true);
        let mut disconnected = instrument(spatial, false);
        let receiver_rest = p.hammers[1].motion;
        for piano in [&mut p, &mut disconnected] {
            piano.silent_key_down(64).unwrap();
            assert_eq!(piano.accounting.input_work_j, 0.0);
            piano.set_sostenuto(true);
            piano.note_off(64).unwrap();
            assert!(!piano.hammers[1].held && piano.hammers[1].latched);
            piano.note_on(60, 2.0).unwrap();
            assert!(!piano.hammers[0].latched, "later strikes must not join the capture");
        }
        let input = p.accounting.input_work_j;
        assert_eq!(input, disconnected.accounting.input_work_j);
        let mut peak = 0.0_f64;
        for _ in 0..1600 {
            p.step().unwrap();
            disconnected.step().unwrap();
            let receiver = &p.bank.strings[1];
            let bridge = receiver.bridge.iter().zip(&p.bank.q[p.bank.modes.len()..])
                .map(|(g, q)| g*q).sum::<f64>();
            for k in receiver.modes.clone() {
                // Fixed-interface deformation, excluding mere endpoint lift.
                peak = peak.max((p.bank.q[k]-p.bank.modes[k].beta*bridge).abs());
            }
        }
        assert!(peak > 1e-10, "shared-board transfer was absent: {peak:e}");
        for k in disconnected.bank.strings[1].modes.clone() {
            assert_eq!(disconnected.bank.q[k], 0.0);
            assert_eq!(disconnected.bank.v[k], 0.0);
        }
        assert!(p.accounting.felt_loss_j > 0.0, "the source must be a real hammer strike");
        assert_eq!(p.hammers[1].motion, receiver_rest);
        assert!(!p.hammers[1].active && p.hammers[1].on_rest);
        for (i, c) in p.contacts.iter().enumerate() {
            if p.bank.strings[p.bank.contact_strings[i]].course == 1 {
                assert_eq!(c.force, 0.0);
                assert_eq!(c.state.eps_max, 0.0);
                assert_eq!(c.memory, relaxation::Memory::default());
            }
        }
        assert_eq!(p.accounting.damper_loss_j, 0.0);
        p.set_sustain(1.0).unwrap();
        p.set_sostenuto(false);
        assert!(!p.hammers[1].held && !p.hammers[1].latched);
        for _ in 0..128 { p.step().unwrap(); }
        assert_eq!(p.accounting.damper_loss_j, 0.0);
        p.set_sustain(0.5).unwrap();
        for _ in 0..256 { p.step().unwrap(); }
        assert!(p.hammers[0].held);
        assert!(p.accounting.damper_loss_j > 0.0, "only the silent receiver's pads are engaged");
        assert_eq!(p.accounting.input_work_j, input);
        let defect = input-p.energy_j()-p.accounting.dissipated_j();
        assert!(defect.abs() < 1e-7, "unaccounted work {defect:e} J");
    }
}
