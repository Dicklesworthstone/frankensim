use super::*;

const HEAD: &str = "frankensim-reed-performance-v1\naudio 48000 4800 20000\n\
    ambient 293.15 101325 0\n";
const TAIL: &str = "listener_m 1\ntermination flanged\nsegments 1\ncylinder 0.0073 0.5\n\
    compile_limits 100000 1024\nschedule\nfrankensim-gesture-schedule-v1\n\
    control_rate_hz\t700\ntracks\t1\ntrack\tblow\tblowing-pressure\n\
    initial\tpressure-pa\t0e0\nevents\t1\nevent\t0e0\t1e-2\tpressure-pa\t2.8e3\n";
const SIX: &str = "reed 0.0004 0.013 6000 3.8e-6 1400 0.077";

fn source(reed_row: &str) -> String {
    format!("{HEAD}{reed_row}\n{TAIL}")
}

fn refusal(reed_row: &str) -> &'static str {
    match ReedPerformance::from_bytes(source(reed_row).as_bytes(), 64) {
        Err(ReedPerformanceError::Input { what, .. }) => what,
        other => panic!(
            "expected an input refusal for {reed_row:?}, got {:?}",
            other.err()
        ),
    }
}

fn pressure(reed_row: &str) -> Vec<f64> {
    let mut performance = ReedPerformance::from_bytes(source(reed_row).as_bytes(), 64)
        .expect("a valid reed row must load");
    let mut out = vec![0.0; 4800];
    for chunk in out.chunks_mut(64) {
        performance
            .block(chunk)
            .expect("the short window must render");
    }
    out
}

#[test]
fn a_four_field_lay_loads_and_the_six_field_row_is_untouched() {
    let plain = ReedPerformance::from_bytes(source(SIX).as_bytes(), 64).expect("six fields");
    let with_lay =
        ReedPerformance::from_bytes(source(&format!("{SIX} 7.9e6 2 2e-4 6.47")).as_bytes(), 64)
            .expect("six fields and a lay");
    assert!(plain.info().massive_reed && with_lay.info().massive_reed);
    assert_ne!(plain.info().input_hash, with_lay.info().input_hash);
}

#[test]
fn a_malformed_lay_is_refused_at_its_field() {
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 2e-4")),
        "missing field",
        "a lay needs all four fields"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 2e-4 6.47 1")),
        "unexpected extra field"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 0.5 2e-4 6.47")),
        "lay exponent must be at least 1 and its loss nonnegative"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 2e-4 -1")),
        "lay exponent must be at least 1 and its loss nonnegative"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 4e-4 6.47")),
        "lay onset must lie between zero and the reed's rest opening"
    );
    assert_eq!(
        refusal(&format!("{SIX} 7.9e6 2 -1e-4 6.47")),
        "lay onset must lie between zero and the reed's rest opening"
    );
    assert_eq!(
        refusal(&format!("{SIX} 0 2 2e-4 6.47")),
        "physical dimension/scale must be positive"
    );
    assert_eq!(
        refusal("reed 0.0004 0.013 6000 0 0 0.35 7.9e6 2 2e-4 6.47"),
        "a measured lay needs a massive reed"
    );
}

#[test]
fn a_lay_that_meets_the_reed_before_the_slit_shuts_changes_the_sound() {
    // At 2.8 kPa on a 6 kPa reed the static opening is about 0.21 mm, so an onset
    // of 0.3 mm keeps the reed in contact for the whole window while the default
    // lay (onset 0) is untouched: the two renders cannot agree.
    let plain = pressure(SIX);
    let curved = pressure(&format!("{SIX} 7.9e6 2 3e-4 6.47"));
    assert!(plain.iter().chain(&curved).all(|p| p.is_finite()));
    let largest = plain
        .iter()
        .zip(&curved)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        largest > 1e-6,
        "the lay left the sound unchanged: {largest}"
    );
    // The default render is deterministic, so the difference above is the lay.
    assert_eq!(pressure(SIX), plain);
}

#[test]
fn tongue_contact_holds_massive_reed_closed_until_release_and_preserves_callback_invariance() {
    let untongued = source(SIX);
    let tongued = source(&format!(
        "{SIX}\ntongue_contact 0.000026 0.0000933333 0.35 0.02 0.005"
    ));
    let render_blocked = |text: &str, block: usize| -> Vec<u64> {
        let mut perf = ReedPerformance::from_bytes(text.as_bytes(), block).unwrap();
        let mut out = vec![0.0; 4800];
        for chunk in out.chunks_mut(block) {
            perf.block(chunk).unwrap();
        }
        out.into_iter().map(f64::to_bits).collect()
    };
    let out_1 = render_blocked(&tongued, 1);
    let out_37 = render_blocked(&tongued, 37);
    let out_512 = render_blocked(&tongued, 512);
    assert_eq!(out_1, out_37);
    assert_eq!(out_1, out_512);
    assert_ne!(out_1, render_blocked(&untongued, 64));
    assert_eq!(
        refusal("reed 0.0004 0.013 6000 0 0 0.35\ntongue_contact 0.000026 0.0000933333 0.35 0.02 0.005"),
        "tongue_contact requires a massive reed"
    );
    assert_eq!(
        refusal(&format!("{SIX}\ntongue_contact 0.0002 0.0000933333 0.35 0.02 0.005")),
        "tongue_contact area exceeds the vibrating reed blade area"
    );
    assert_eq!(
        refusal(&format!("{SIX}\ntongue_contact 0.000026 0.0000933333 0.35 -0.01 0.005")),
        "tongue_contact requires a nonnegative release_start_s"
    );
    assert_eq!(
        refusal(&format!("{SIX}\ntongue_contact 0.000026 0.0000933333 0.35 0.02 0")),
        "physical dimension/scale must be positive"
    );
}

/// The tongue's force history was lifted out of `append_controls` into `TongueContactInput::force_series` so the DG-valve
/// file's typed reed compiles the same numbers. The reference below is the loop as it stood before the lift, verbatim; the
/// series and the controls built from it must equal it value for value and sample for sample.
fn assert_series_is_the_original_loop(
    schedule_text: &str,
    reed: BeatingReed,
    audio_rate_hz: u32,
    last_tick: u64,
    tongue: TongueContactInput,
) {
    let schedule = GestureSchedule::from_canonical_bytes(schedule_text.as_bytes())
        .expect("a valid gesture schedule");
    let face = crate::reed_bore::reed_pressure_face(reed);
    let control_rate_hz = schedule.control_rate_hz;
    let release_tick = ((tongue.release_start_s * f64::from(control_rate_hz)).round() as u64).min(last_tick);
    let release_pressure_pa = schedule.sample("blow", release_tick).unwrap();
    let pre_release_close_n = tongue.modal_factor * tongue.force_n;
    let held_downward_n = pre_release_close_n + face * release_pressure_pa;
    let mut previous = None;
    let mut reference = Vec::new();
    for tick in 0..=last_tick {
        let t = tick as f64 / f64::from(control_rate_hz);
        let force_n = if t <= tongue.release_start_s {
            let p_m = schedule.sample("blow", tick).unwrap();
            (held_downward_n - face * p_m).max(0.0)
        } else {
            let elapsed = t - tongue.release_start_s;
            if elapsed < tongue.release_s {
                (pre_release_close_n * (1.0 - elapsed / tongue.release_s)).max(0.0)
            } else {
                0.0
            }
        };
        let bits = force_n.to_bits();
        if previous != Some(bits) {
            reference.push((tick, force_n));
            previous = Some(bits);
        }
    }
    let series = tongue.force_series(&schedule, "blow", face, last_tick).unwrap();
    assert_eq!(series.len(), reference.len());
    for (a, b) in series.iter().zip(&reference) {
        assert_eq!((a.0, a.1.to_bits()), (b.0, b.1.to_bits()));
    }
    assert!(series.len() > 2 && series.last().unwrap().1 == 0.0, "a tongue lets go");
    let mut controls = Vec::new();
    tongue
        .append_controls(&mut controls, &schedule, "blow", reed, audio_rate_hz, last_tick)
        .unwrap();
    assert_eq!(controls.len(), reference.len());
    for (control, (tick, force_n)) in controls.iter().zip(&reference) {
        let sample = (u128::from(*tick) * u128::from(audio_rate_hz)).div_ceil(u128::from(control_rate_hz)) as u64;
        assert_eq!(control.sample, sample);
        assert_eq!(control.delta, ControlDelta::SetReedClosingForce { voice: 0, force_n: *force_n });
    }
}

#[test]
fn the_shared_tongue_force_series_is_the_original_loop_and_the_controls_are_its_mapping() {
    let text = "frankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\ntrack\tblow\tblowing-pressure\n\
         initial\tpressure-pa\t0e0\nevents\t1\nevent\t0e0\t1e-2\tpressure-pa\t2.8e3\n";
    let reed = BeatingReed {
        rest_opening_m: 0.0004,
        width_m: 0.013,
        closing_pressure_pa: 6000.0,
        mass_kg: 3.8e-6,
        stiffness_n_m: 1400.0,
        damping_ratio: 0.077,
        blowing_pressure_pa: 0.0,
        attack_s: 0.0,
    };
    let last_tick = 4799 * 700 / 48_000;
    for (area, force, start, release) in
        [(0.000026, 0.35, 0.02, 0.005), (0.000026, 0.2, 0.0, 0.01), (0.00001, 0.5, 0.05, 0.002)]
    {
        let tongue = TongueContactInput::new(area, 0.0000933333, force, start, release).unwrap();
        assert_series_is_the_original_loop(text, reed, 48_000, last_tick, tongue);
    }
}

/// The oboe F4 page's own tongue card: `reed 0.000357 0.0056 17720 2.2e-5 3464 0.114`, `tongue_contact 1.12e-5 6.979e-5
/// 0.2 0.1845 0.005`, audio 96 kHz, 254 400 frames, a 48 kHz control rate and the card's blowing-pressure schedule (a 20 ms
/// ramp to 6.865 kPa at 0.15 s, release at 2.45 s): the control list the reed-performance route built before the lift is
/// the list it builds now, sample for sample and bit for bit.
#[test]
fn the_f4_tongue_card_builds_the_same_controls_as_before_the_series_was_shared() {
    let text = "frankensim-gesture-schedule-v1\ncontrol_rate_hz\t48000\ntracks\t1\ntrack\tblow\tblowing-pressure\n\
         initial\tpressure-pa\t0e0\nevents\t2\nevent\t1.5e-1\t2e-2\tpressure-pa\t6.865e3\n\
         event\t2.45e0\t6e-2\tpressure-pa\t0e0\n";
    let reed = BeatingReed {
        rest_opening_m: 0.000357,
        width_m: 0.0056,
        closing_pressure_pa: 17720.0,
        mass_kg: 2.2e-5,
        stiffness_n_m: 3464.0,
        damping_ratio: 0.114,
        blowing_pressure_pa: 0.0,
        attack_s: 0.0,
    };
    let tongue = TongueContactInput::new(1.12e-5, 6.979e-5, 0.2, 0.1845, 0.005).unwrap();
    let last_tick = 254_399_u64 * 48_000 / 96_000;
    assert_eq!(last_tick, 127_199);
    assert_series_is_the_original_loop(text, reed, 96_000, last_tick, tongue);
}
