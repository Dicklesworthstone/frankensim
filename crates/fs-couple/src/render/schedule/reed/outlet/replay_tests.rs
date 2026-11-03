//! The baffled clarinet observer no longer replays one window late (14773/14776).
//!
//! Before this change every clarinet site card's outlet kernel kept the
//! transfer's pre-arrival spread at the END of its 16384-tap window. At the
//! fundamental that wrapped part was as loud as the direct path (G3 bell
//! +0.1 dB, D5 +1.0 dB), so after the breath-off the observer replayed the
//! note for one more window (0.341 s) and then cut to zero.
use super::*;
use crate::driving_point::regional::outlet::{OUTLET_REPLAY_GAIN_BOUND, replay_gain};

const G3: &str = include_str!(
    "../../../../../../../data/instruments/descriptions/clarinet/performances/BbClar.mf.G3.performance"
);
const D5: &str = include_str!(
    "../../../../../../../data/instruments/descriptions/clarinet/performances/BbClar.mf.D5.performance"
);

/// Exact impulse response of a primed observer: the cascade is linear, so an
/// impulse minus the same history fed zeros cancels the prime.
fn impulse_response(observer: &OutletPressure, len: usize) -> Vec<f64> {
    let (mut pulse, mut quiet) = (observer.clone(), observer.clone());
    (0..len)
        .map(|i| pulse.step(if i == 0 { 1.0 } else { 0.0 }).unwrap() - quiet.step(0.0).unwrap())
        .collect()
}

/// The circular window (arrival at tap 0) recovered from the realized observer:
/// nothing precedes `axial + exterior.0`; fold the response modulo the grid from
/// there and undo the latency rotation. The short Rayleigh stage is a common
/// factor at every bin, so the replay ratios are those of the kernel itself.
fn circular_window(observer: &OutletPressure) -> Vec<f64> {
    let info = observer.info();
    let n = info.fft_samples;
    let start = info.axial_delay_samples + info.exterior_delay_samples.0;
    let h = impulse_response(observer, start + n + info.exterior_delay_samples.1 + 2);
    let mut rotated = vec![0.0; n];
    for (i, v) in h.iter().enumerate().skip(start) {
        rotated[(i - start) % n] += v;
    }
    let latency = info.pre_arrival_latency_samples;
    let mut window = rotated[latency..].to_vec();
    window.extend_from_slice(&rotated[..latency]);
    window
}

/// |replayed| / |direct| at one frequency, with the same split as the gate.
fn gain_at(window: &[f64], latency: usize, frequency_hz: f64, rate: f64) -> f64 {
    let n = window.len();
    let (mut direct, mut late) = ((0.0, 0.0), (0.0, 0.0));
    for (j, &v) in window.iter().enumerate() {
        let lag = if j >= n - latency {
            j as f64 - n as f64
        } else {
            j as f64
        };
        let phase = -core::f64::consts::TAU * frequency_hz * lag / rate;
        let term = (v * phase.cos(), v * phase.sin());
        let part = if j >= n / 2 && j < n - latency {
            &mut late
        } else {
            &mut direct
        };
        part.0 += term.0;
        part.1 += term.1;
    }
    late.0.hypot(late.1) / direct.0.hypot(direct.1)
}

fn db(gain: f64) -> f64 {
    20.0 * gain.log10()
}

#[test]
fn clarinet_outlet_observers_admit_no_audible_one_window_replay() {
    for (label, text, f0) in [("G3", G3, 195.74), ("D5", D5, 589.60)] {
        let performance = ReedPerformance::from_bytes(text.as_bytes(), 512)
            .unwrap_or_else(|error| panic!("{label}: {error}"));
        let terminal = performance
            .terminal
            .as_ref()
            .expect("the card selects the baffled outlet");
        let bell = terminal.pressure.info();
        let rate = f64::from(bell.sample_rate_hz);
        let band = (
            0.25 * rate / bell.axial_delay_samples as f64,
            bell.receiver.maximum_frequency_hz,
        );
        for (index, observer) in core::iter::once(&terminal.pressure)
            .chain(&terminal.holes)
            .enumerate()
        {
            let info = observer.info();
            // The receiver gate admits both cards on their requested grid: no
            // regrid, so no added latency and the same output.
            assert_eq!(info.fft_samples, 16_384, "{label} #{index}");
            assert_eq!(info.pre_arrival_latency_samples, info.fft_samples / 4);
            assert!(
                info.replay_gain.0 <= OUTLET_REPLAY_GAIN_BOUND,
                "{label} #{index}: {:?}",
                info.replay_gain
            );
            let window = circular_window(observer);
            let mut line = format!(
                "{label} mouth {index} (0 = bell): admitted worst {:.1} dB at {:.1} Hz;",
                db(info.replay_gain.0),
                info.replay_gain.1
            );
            for latency in [0, 64, 256, 1024, 4096] {
                let (worst, at) = replay_gain(&window, latency, bell.sample_rate_hz, band);
                let harmonics = (1..=5)
                    .map(|h| format!("{:.1}", db(gain_at(&window, latency, h as f64 * f0, rate))))
                    .collect::<Vec<_>>()
                    .join(" ");
                line += &format!(
                    " L={latency}: band worst {:.1} dB at {at:.0} Hz, f0..5f0 [{harmonics}];",
                    db(worst)
                );
            }
            println!("{line}");
            // The production latency holds the bound at the note's own harmonics too.
            for h in 1..=5 {
                let gain = gain_at(
                    &window,
                    info.pre_arrival_latency_samples,
                    h as f64 * f0,
                    rate,
                );
                assert!(
                    gain <= OUTLET_REPLAY_GAIN_BOUND,
                    "{label} #{index} harmonic {h}: {:.1} dB",
                    db(gain)
                );
            }
            // Planted negative: the old whole-delay removal (no latency) fails on the G3 bell.
            if label == "G3" && index == 0 {
                let (old, at) = replay_gain(&window, 0, bell.sample_rate_hz, band);
                assert!(
                    old > OUTLET_REPLAY_GAIN_BOUND,
                    "old removal on the G3 bell: {:.1} dB at {at} Hz",
                    db(old)
                );
                assert!(
                    gain_at(&window, 0, f0, rate) > 0.1,
                    "old removal replays G3's f0 near the direct level"
                );
            }
        }
    }
}

/// The known hold: G3 with the breath-off moved to 1.00 s held -2 dB (re the
/// loudest 10 ms frame) until 1.34 s, then fell to -37 dB by 1.41 s. Now the
/// level must be at least 20 dB down over 1.30-1.40 s (the observer latency
/// is 85 ms of that). Renders 1.5 s of the full card: minutes in release.
#[test]
#[ignore = "renders 1.5 s of the G3 site card; run by name"]
fn clarinet_g3_release_no_longer_holds_one_window() {
    let text = G3.replace(
        "event\t2.45e0\t6e-2\tpressure-pa\t0e0\n",
        "event\t1e0\t6e-2\tpressure-pa\t0e0\n",
    );
    assert_ne!(text, G3);
    let mut performance = ReedPerformance::from_bytes(text.as_bytes(), 512).unwrap();
    let mut out = vec![0.0; 72_000];
    for part in out.chunks_mut(480) {
        performance.pressure_block(part).unwrap();
    }
    let rms: Vec<f64> = out
        .chunks(480)
        .map(|c| (c.iter().map(|v| v * v).sum::<f64>() / 480.0).sqrt())
        .collect();
    let top = rms.iter().copied().fold(0.0, f64::max);
    let level = |frame: usize| 20.0 * (rms[frame] / top).log10();
    println!(
        "G3 release@1.00: 10 ms dB 0.90-1.50 s: {}",
        (90..150)
            .map(|i| format!("{:.0}", level(i)))
            .collect::<Vec<_>>()
            .join(" ")
    );
    println!(
        "G3 sustained 0.50-0.90 s: mean {:.2} dB re the loudest frame",
        (50..90).map(level).sum::<f64>() / 40.0
    );
    for frame in 130..140 {
        assert!(
            level(frame) <= -20.0,
            "frame {frame}: {:.1} dB, still holding",
            level(frame)
        );
    }
}
