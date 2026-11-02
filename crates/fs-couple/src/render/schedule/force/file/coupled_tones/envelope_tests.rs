use super::*;

const EXAMPLE: &[u8] =
    include_bytes!("../../../../../../../../examples/coupled_enveloped_tones.performance");

#[test]
fn coupled_envelopes_keep_callback_phase_and_leave_physical_ringdown_after_release() {
    let text = std::str::from_utf8(EXAMPLE).unwrap();
    let prefix = text.split("envelopes ").next().unwrap();
    let receiver = &text[text.find("coupled_broadband ").unwrap()..];
    let plain = format!("{prefix}{receiver}");
    let mut whole = CoupledTonePerformance::from_bytes(EXAMPLE, 4800).unwrap();
    let mut split = CoupledTonePerformance::from_bytes(EXAMPLE, 37).unwrap();
    let mut carrier = CoupledTonePerformance::from_bytes(plain.as_bytes(), 4800).unwrap();
    assert_eq!(whole.envelope_count(), 3);
    let mut a = vec![0.0; 4800];
    let mut b = vec![0.0; 4800];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(37) {
        split.block(chunk).unwrap();
    }
    carrier.block(&mut vec![0.0; 4800]).unwrap();
    assert_eq!(a, b);
    assert_eq!(whole.cosine, split.cosine);
    assert_eq!(whole.sine, split.sine);
    // Carrier propagation is independent of gain, including the silent gap
    // and re-attack, and retains original source tone order through coupling.
    assert_eq!(whole.cosine, carrier.cosine);
    assert_eq!(whole.sine, carrier.sine);
    assert_eq!(whole.held, carrier.held);
    assert_eq!(whole.network().last_frame(), split.network().last_frame());
    assert_eq!(
        whole.info().input_hash,
        fs_blake3::hash_domain(MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN, EXAMPLE)
    );
    assert!(a.iter().any(|p| p.abs() > 1e-7));

    // Release every held input as well. After sample3200 all tone gains are
    // zero, but retained structural motion and receiver history must continue.
    let released = text.replace("force 613 1 0 0.04", "force 613 1 0 0");
    let mut ring = CoupledTonePerformance::from_bytes(released.as_bytes(), 3200).unwrap();
    ring.block(&mut vec![0.0; 3200]).unwrap();
    assert_eq!(ring.held, vec![0.0, 0.0]);
    for tone in 0..3 {
        assert_eq!(
            ring.envelopes
                .as_ref()
                .unwrap()
                .endpoints(tone, 3200)
                .unwrap(),
            [0.0, 0.0]
        );
    }
    let before = ring.network().total_energy_j().unwrap();
    let mut tail = vec![0.0; 1600];
    ring.block(&mut tail).unwrap();
    let after = ring.network().total_energy_j().unwrap();
    assert!(before > after && after > 0.0);
    assert!(tail.iter().any(|p| p.abs() > 1e-9));
    let frame = ring.network().last_frame().unwrap();
    assert_eq!(frame.external_work_j, 0.0);
    assert!(frame.total_dissipation_j() > 0.0);
    assert_eq!(frame.numerical_dissipation_j, 0.0);
}

#[test]
fn invalid_knots_refuse_and_failed_coupled_sample_keeps_accepted_envelope_cursor() {
    let valid = std::str::from_utf8(EXAMPLE).unwrap();
    for invalid in [
        valid.replacen("gain 0 0", "gain 1 0", 1),
        valid.replace("envelope 1 4", "envelope 0 4"),
        valid.replace("gain 4800 0", "gain 4801 0"),
        valid.replace("gain 240 1", "gain 240 NaN"),
        valid.replace("gain 960 0", "gain 240 0"),
    ] {
        assert!(CoupledTonePerformance::from_bytes(invalid.as_bytes(), 240).is_err());
    }
    let failing = valid.replace("events 3\n", "events 4\nforce 240 0 0 1e20\n");
    let mut performance = CoupledTonePerformance::from_bytes(failing.as_bytes(), 240).unwrap();
    performance.block(&mut vec![0.0; 240]).unwrap();
    let phase = (performance.cosine.clone(), performance.sine.clone());
    let held = performance.held.clone();
    let event = performance.next_event;
    let gains = performance
        .envelopes
        .as_ref()
        .unwrap()
        .endpoints(0, 240)
        .unwrap();
    let frame = performance.network().last_frame().cloned();
    let states: Vec<_> = performance
        .network()
        .components()
        .iter()
        .map(|model| model.states().to_vec())
        .collect();
    assert!(performance.block(&mut [0.0]).is_err());
    assert_eq!(performance.samples_rendered(), 240);
    assert_eq!(performance.network().samples_rendered(), 240);
    assert_eq!(performance.network().last_frame(), frame.as_ref());
    assert_eq!(performance.next_event, event);
    assert_eq!(performance.held, held);
    assert_eq!(
        (performance.cosine.clone(), performance.sine.clone()),
        phase
    );
    assert_eq!(
        performance
            .envelopes
            .as_ref()
            .unwrap()
            .endpoints(0, 240)
            .unwrap(),
        gains
    );
    for (model, state) in performance.network().components().iter().zip(states) {
        assert_eq!(model.states(), state);
    }
}
