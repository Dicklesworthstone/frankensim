use super::*;
use fs_math::c64::C64;

const EXAMPLE: &[u8] =
    include_bytes!("../../../../../../../../examples/coupled_broadband_tones.performance");

#[test]
fn complete_receiver_file_keeps_signed_modal_rows_delay_and_callback_phase() {
    let mut whole = CoupledTonePerformance::from_bytes(EXAMPLE, 4800).unwrap();
    let mut split = CoupledTonePerformance::from_bytes(EXAMPLE, 37).unwrap();
    let report = whole.broadband_report().unwrap();
    assert_eq!(report.source_id, "example:signed-coupled-one-pole");
    assert_eq!(report.delay_samples, 3.0);
    assert_eq!(report.sample_rate_hz, 48000);
    assert!(report.held_out_max_rel_error <= report.controls.max_held_out_rel_error);
    assert!(!report.constant_narrow_band_compatibility);
    assert_eq!(whole.continuous_force_tones(), 3);
    assert_eq!(whole.compiled_controls(), 3);
    let broadband::Source::Broadband(receiver) = &whole.source else {
        panic!("receiver not constructed")
    };
    // Compare every signed, flattened modal row with the declared delayed
    // rational operator. This checks the file's row ordering and convention;
    // the receiver owner's tests cover its independent time/filter oracle.
    for omega in [175.0, 750.0, 1800.0, 4411.424404170788, 7800.0, 11000.0] {
        let phase: f64 = omega * 3.0 / 48000.0;
        let delay = C64::new(phase.cos(), phase.sin());
        for (mode, (direct, gain)) in [(0.25, 1.5), (-0.4, -0.6)].into_iter().enumerate() {
            let expected = (C64::from_re(direct)
                + C64::from_re(gain * 900.0) / C64::new(900.0, -omega))
                * delay;
            let actual = receiver.transfer_at(mode, omega).unwrap();
            let tolerance = 2e-6 * expected.re.hypot(expected.im);
            assert!((actual.re - expected.re).hypot(actual.im - expected.im) <= tolerance);
        }
    }
    let mut a = vec![0.0; 4800];
    let mut b = vec![0.0; 4800];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(37) {
        split.block(chunk).unwrap();
    }
    assert_eq!(a, b);
    assert_eq!(&a[..3], &[0.0; 3]);
    assert!(a[3..].iter().any(|p| p.abs() > 1e-6));
    assert_eq!(whole.cosine, split.cosine);
    assert_eq!(whole.sine, split.sine);
    assert_eq!(whole.held, vec![0.0, -0.7 * 0.04]);
    assert_eq!(whole.next_event, whole.events.len());
    assert_eq!(whole.network().last_frame(), split.network().last_frame());
    assert_eq!(
        whole.info().input_hash,
        fs_blake3::hash_domain(MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN, EXAMPLE)
    );
}

fn without_receiver(text: &str) -> String {
    text.split("coupled_broadband ")
        .next()
        .unwrap()
        .lines()
        .map(|line| {
            let mut fields: Vec<_> = line.split_ascii_whitespace().collect();
            if fields.first() == Some(&"mode") {
                fields[3] = "0";
                fields[4] = "0";
            }
            fields.join(" ") + "\n"
        })
        .collect()
}

#[test]
fn receiver_retains_every_accepted_source_state_event_and_work_heat_record() {
    let text = std::str::from_utf8(EXAMPLE).unwrap();
    // A supplied legacy quadrature must be observation-only metadata when an
    // explicit causal receiver owns pressure; it cannot affect the mechanics.
    let with_legacy = text
        .replace("0.01 2 0 0 0", "0.01 2 0.3 0 0")
        .replace("0.02 -1 0 0 0", "0.02 -1 -0.2 0 0");
    let mut observed = CoupledTonePerformance::from_bytes(with_legacy.as_bytes(), 1).unwrap();
    let neutral = without_receiver(text);
    let mut physical = CoupledTonePerformance::from_bytes(neutral.as_bytes(), 1).unwrap();
    assert!(physical.broadband_report().is_none());
    let initial = physical.network().total_energy_j().unwrap();
    assert!(initial > 0.0); // Complete spring preload remains in both paths.
    let mut work = 0.0;
    let mut heat = 0.0;
    for sample in 1..=4800 {
        observed.block(&mut [0.0]).unwrap();
        physical.block(&mut [0.0]).unwrap();
        assert_eq!(
            observed.network().last_frame(),
            physical.network().last_frame()
        );
        assert_eq!(observed.held, physical.held);
        assert_eq!(observed.cosine, physical.cosine);
        assert_eq!(observed.sine, physical.sine);
        for (a, b) in observed
            .network()
            .components()
            .iter()
            .zip(physical.network().components())
        {
            assert_eq!(a.states(), b.states());
            assert!(
                a.modes()
                    .iter()
                    .all(|mode| mode.pressure_per_modal_velocity == C64::ZERO)
            );
        }
        let frame = observed.network().last_frame().unwrap();
        assert_eq!(frame.sample, sample);
        assert_eq!(frame.numerical_dissipation_j, 0.0);
        work += frame.external_work_j;
        heat += frame.total_dissipation_j();
    }
    let final_energy = observed.network().total_energy_j().unwrap();
    assert!(heat > 0.0);
    assert!((final_energy - initial + heat - work).abs() < 1e-11);
    assert_eq!(observed.samples_rendered(), 4800);
}

fn change_first_row(text: &str, kind: &str, edit: impl FnOnce(&mut Vec<String>)) -> String {
    let mut edit = Some(edit);
    text.lines()
        .map(|line| {
            let mut fields: Vec<_> = line.split_ascii_whitespace().map(str::to_owned).collect();
            if fields.first().is_some_and(|key| key == kind)
                && let Some(edit) = edit.take()
            {
                edit(&mut fields);
            }
            fields.join(" ") + "\n"
        })
        .collect()
}

#[test]
fn incomplete_signed_receiver_clock_fit_and_preparation_budgets_refuse() {
    let valid = std::str::from_utf8(EXAMPLE).unwrap();
    let invalid = [
        change_first_row(valid, "frequency", |fields| {
            fields.pop();
        }),
        change_first_row(valid, "frequency", |fields| fields.push("0".into())),
        change_first_row(valid, "withheld", |fields| fields[2] = "100".into()),
        change_first_row(valid, "withheld", |fields| fields[1] = "200000".into()),
        change_first_row(valid, "coupled_broadband", |fields| fields[1] = "0".into()),
        change_first_row(valid, "coupled_broadband", |fields| fields[2] = "1".into()),
        change_first_row(valid, "coupled_broadband", |fields| fields[3] = "1".into()),
        change_first_row(valid, "transfer", |fields| fields[2] = "-1".into()),
        format!("{valid}unexpected receiver tail\n"),
    ];
    for bad in invalid {
        assert!(
            CoupledTonePerformance::from_bytes(bad.as_bytes(), 32).is_err(),
            "admitted {bad}"
        );
    }
    let admitted = CoupledTonePerformance::from_bytes(EXAMPLE, 32).unwrap();
    assert!(admitted.validate_sample_rate(44100).is_err());
    assert!(admitted.validate_sample_count(4801).is_err());
}
