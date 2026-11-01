//! The realized bank against the image it is realized from, on the two heads as
//! they stand in the repository. One smoothed flow step drives the same-sample
//! reflection. Its first difference is a pulse, so the ratio of transforms is
//! the input impedance the loop plays. The first ten peaks of that curve must
//! sit within 2 cents and 0.3 dB of the image's own. A second test reads the
//! bank's peaks in air and heliox at equal c dt against the images' ratios.

use super::*;
use crate::mm_line::{CupState, cup_junction_boundary};
use fs_fft::{C64 as FftC64, Fft};
use fs_math::det;

const TROMBONE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../data/instruments/descriptions/trombone/tenor-trombone.head"
));
const TRUMPET: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../data/instruments/descriptions/trumpet/trumpet-open.head"
));
// The parser wants a gesture. Nothing here plays it.
const SCHEDULE: &str = "schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\n\
    tracks\t2\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t0\n\
    track\tlip\tpre-stress\ninitial\tmultiplier\t1e0\nevents\t0\n";
const FLOW_M3_S: f64 = 1e-7;
const RAMP_S: f64 = 1e-4;
// The pass rests on this window: at 2 s the trombone reads 0.08 cents and 0.01
// dB. Cut to 0.5 s it fails at 774 cents and 40.6 dB, because the cut stops the
// response before it decays and leaves a false peak at 73 Hz. Don't shorten it.
const SETTLE_S: f64 = 2.0;
const TRANSFORM: usize = 1 << 17;
const IMAGE_FIRST_HZ: f64 = 20.0;
const IMAGE_STEP_HZ: f64 = 0.25;
const PEAKS: usize = 10;

fn step_flow(time_s: f64) -> f64 {
    if time_s <= 0.0 {
        0.0
    } else if time_s < RAMP_S {
        FLOW_M3_S * 0.5 * (1.0 - det::cos(core::f64::consts::PI * time_s / RAMP_S))
    } else {
        FLOW_M3_S
    }
}

// The first difference, zero-padded and transformed.
fn pulse_spectrum(series: &[f64]) -> Vec<FftC64> {
    let mut buf = vec![FftC64::new(0.0, 0.0); TRANSFORM];
    let mut previous = 0.0;
    for (slot, value) in buf.iter_mut().zip(series) {
        *slot = FftC64::new(value - previous, 0.0);
        previous = *value;
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); TRANSFORM];
    Fft::new(TRANSFORM).forward(&mut buf, &mut scratch);
    buf
}

/// First peaks above 45 Hz, each hz and dB from a parabola through log magnitude.
fn peaks(first_hz: f64, step_hz: f64, magnitude: &[f64]) -> Vec<(f64, f64)> {
    let mut found = Vec::new();
    for k in 2..magnitude.len() - 2 {
        let hz = first_hz + k as f64 * step_hz;
        if hz > 45.0 && magnitude[k] > magnitude[k - 1] && magnitude[k] >= magnitude[k + 1] {
            let (la, lb, lc) = (
                magnitude[k - 1].ln(),
                magnitude[k].ln(),
                magnitude[k + 1].ln(),
            );
            let offset = 0.5 * (la - lc) / (la - 2.0 * lb + lc);
            found.push((
                hz + offset * step_hz,
                20.0 * (lb - 0.25 * (la - lc) * offset) / core::f64::consts::LN_10,
            ));
            if found.len() == PEAKS {
                break;
            }
        }
    }
    found
}

fn check(tag: &str, rate_hz: u32, head: &str) {
    let text =
        format!("frankensim-lip-performance-v1\naudio {rate_hz} {rate_hz}\n{head}{SCHEDULE}");
    let description = LipDescription::from_bytes(text.as_bytes()).unwrap();
    let gas = *description.ambient_gas();
    let config = MmLineConfig {
        sample_rate_hz: description.rate_hz,
        n_modes: description.modes,
        extra_slices: description.extra_slices,
        loss: LossModel::Bessel,
    };
    let load = MmLoad::Analytic(description.termination);
    let mut voice = BrassVoice::new(
        core::slice::from_ref(&description.duct),
        &["held"],
        &gas,
        &load,
        &config,
        description.lip.clone(),
        description.cup_volume_m3,
    )
    .unwrap();
    let zc0 = voice.port_impedance();
    let dt = 1.0 / f64::from(rate_hz);
    let compliance = voice.cup_compliance;
    let samples = (SETTLE_S / dt) as usize;

    let mut cup = CupState::default();
    let mut flow = Vec::with_capacity(samples);
    let mut pressure = Vec::with_capacity(samples);
    let mut worst_read = 0.0f64;
    for n in 0..samples {
        let volume_flow = step_flow(n as f64 * dt);
        let boundary = voice.bank.boundary_response().unwrap();
        let (outgoing, reflected, bore) =
            cup_junction_boundary(volume_flow, boundary, &mut cup, compliance, zc0, dt);
        let actual = voice.bank.push(outgoing).unwrap();
        worst_read = worst_read.max((actual - reflected).abs() / (1.0 + actual.abs()));
        flow.push(volume_flow);
        pressure.push(bore);
    }
    assert!(worst_read < 1e-9, "{tag}: the affine read left the push");

    let df = f64::from(rate_hz) / TRANSFORM as f64;
    let top = (1500.0 / df) as usize;
    let (pressure, flow) = (pulse_spectrum(&pressure), pulse_spectrum(&flow));
    let played: Vec<f64> = (0..top)
        .map(|k| pressure[k].re.hypot(pressure[k].im) / flow[k].re.hypot(flow[k].im))
        .collect();
    let played = peaks(0.0, df, &played);

    let image_hz: Vec<f64> = (0..5921)
        .map(|k| IMAGE_FIRST_HZ + f64::from(k) * IMAGE_STEP_HZ)
        .collect();
    let image: Vec<f64> = description
        .bore_scan(&image_hz)
        .unwrap()
        .rows
        .iter()
        .map(|row| row.impedance_multimodal[0].hypot(row.impedance_multimodal[1]))
        .collect();
    let image = peaks(IMAGE_FIRST_HZ, IMAGE_STEP_HZ, &image);

    assert_eq!(image.len(), PEAKS, "{tag}: the image has too few peaks");
    assert_eq!(played.len(), PEAKS, "{tag}: the bank has too few peaks");
    let (mut worst_cents, mut worst_db) = (0.0f64, 0.0f64);
    for (image, played) in image.iter().zip(&played) {
        worst_cents = worst_cents.max((1200.0 * (played.0 / image.0).log2()).abs());
        worst_db = worst_db.max((played.1 - image.1).abs());
    }
    // The card receipt's own read, taps against the spectrum they were
    // realized from on the bank's bins, printed beside the played one.
    let receipt = &voice.receipts()[0];
    println!(
        "bank limit {tag}: worst {worst_cents:.2} cents, {worst_db:.2} dB; receipt {:.2} cents, {:.2} dB over {} peaks",
        receipt.image_worst_cents, receipt.image_worst_db, receipt.image_peaks
    );
    assert!(
        worst_cents <= 2.0 && worst_db <= 0.3,
        "{tag}: worst {worst_cents:.2} cents, {worst_db:.2} dB"
    );
}

// This costs about 14 minutes single-threaded on a loaded host, mostly the
// trombone's 65536 taps. It stays a plain test on purpose.
#[test]
fn bank_same_sample_limit() {
    check("trumpet", 32_000, TRUMPET);
    check("trombone", 19_200, TROMBONE);
}

/// The first ten |Z| peaks of the bank's own taps, Z = zc0 (1 + R) / (1 - R),
/// and of the image the bank is realized from (the multimodal column of
/// `bore_scan`), up to `band_hz`; and the gas's sound speed. Read from the taps,
/// the bank sits within 0.09 c and 0.01 dB of the same bank played by the loop.
fn bank_and_image(
    rate_hz: u32,
    head: &str,
    band_hz: f64,
) -> (Vec<(f64, f64)>, Vec<(f64, f64)>, f64) {
    let text =
        format!("frankensim-lip-performance-v1\naudio {rate_hz} {rate_hz}\n{head}{SCHEDULE}");
    let description = LipDescription::from_bytes(text.as_bytes()).unwrap();
    let gas = *description.ambient_gas();
    let config = MmLineConfig {
        sample_rate_hz: description.rate_hz,
        n_modes: description.modes,
        extra_slices: description.extra_slices,
        loss: LossModel::Bessel,
    };
    let load = MmLoad::Analytic(description.termination);
    let voice = BrassVoice::new(
        core::slice::from_ref(&description.duct),
        &["held"],
        &gas,
        &load,
        &config,
        description.lip.clone(),
        description.cup_volume_m3,
    )
    .unwrap();
    let zc0 = voice.port_impedance();
    let mut buf = vec![FftC64::new(0.0, 0.0); TRANSFORM];
    for (slot, tap) in buf.iter_mut().zip(voice.bank.combo_taps(0)) {
        *slot = FftC64::new(*tap, 0.0);
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); TRANSFORM];
    Fft::new(TRANSFORM).forward(&mut buf, &mut scratch);
    let df = f64::from(rate_hz) / TRANSFORM as f64;
    let bank: Vec<f64> = buf[..(band_hz / df) as usize]
        .iter()
        .map(|r| zc0 * (1.0 + r.re).hypot(r.im) / (1.0 - r.re).hypot(r.im))
        .collect();
    let image_hz: Vec<f64> = (0..=((band_hz - IMAGE_FIRST_HZ) / IMAGE_STEP_HZ) as usize)
        .map(|k| IMAGE_FIRST_HZ + k as f64 * IMAGE_STEP_HZ)
        .collect();
    let image: Vec<f64> = description
        .bore_scan(&image_hz)
        .unwrap()
        .rows
        .iter()
        .map(|row| row.impedance_multimodal[0].hypot(row.impedance_multimodal[1]))
        .collect();
    (
        peaks(0.0, df, &bank),
        peaks(IMAGE_FIRST_HZ, IMAGE_STEP_HZ, &image),
        gas.sound_speed,
    )
}

/// Each realized gas/air peak ratio against the predicted one, in cents: Ok
/// when all ten are within 1 c.
fn ratios(gas: &[(f64, f64)], air: &[(f64, f64)], predicted: &[f64]) -> Result<Vec<f64>, Vec<f64>> {
    let cents: Vec<f64> = gas
        .iter()
        .zip(air)
        .zip(predicted)
        .map(|((g, a), p)| 1200.0 * (g.0 / a.0 / p).log2())
        .collect();
    if cents.len() == PEAKS && cents.iter().all(|c| c.abs() <= 1.0) {
        Ok(cents)
    } else {
        Err(cents)
    }
}

// Ruling 10707: the declared gas reaches the rendered bank. Air at 32 kHz and
// heliox at 59,154 Hz have the same c dt (10.73 mm), so the same bank error. The
// prediction is the image each bank is realized from, the multimodal column of
// `bore_scan` with each gas's own c; the scalar column moves the ratio by at most
// 0.16 c on the first six peaks. The 1 c limit is a regression limit, set from
// the measured 0.13 c between the two banks' per-peak errors, so a pass is not
// independent evidence. The planted negative, the air bank on both sides, misses
// by about 1,000 c. It took 496 s on ts1 at bcb39833c (worst -0.14 c), past the
// ruled five minutes, so it runs by name.
#[test]
#[ignore = "496 s; run: cargo test --release -p fs-couple --lib -- --ignored the_bore_peaks_move_with_the_gas_in_the_rendered_bank"]
fn the_bore_peaks_move_with_the_gas_in_the_rendered_bank() {
    const AIR_ROW: &str = "ambient moist_air 293.15 101325 0";
    const HELIOX_ROW: &str = "ambient mixture 293.15 101325 mole 2 He 0.8 O2 0.2";
    assert_eq!(TRUMPET.matches(AIR_ROW).count(), 1);
    let (air, air_image, c_air) = bank_and_image(32_000, TRUMPET, 1500.0);
    let heliox_head = TRUMPET.replace(AIR_ROW, HELIOX_ROW);
    let (heliox, heliox_image, c_heliox) = bank_and_image(59_154, &heliox_head, 2800.0);
    assert!(
        (59_154.0 - 32_000.0 * c_heliox / c_air).abs() <= 0.5,
        "equal c dt wants {} Hz",
        32_000.0 * c_heliox / c_air
    );
    for (name, found) in [
        ("air bank", &air),
        ("air image", &air_image),
        ("heliox bank", &heliox),
        ("heliox image", &heliox_image),
    ] {
        assert_eq!(found.len(), PEAKS, "{name}: {found:?}");
    }
    let predicted: Vec<f64> = heliox_image
        .iter()
        .zip(&air_image)
        .map(|(h, a)| h.0 / a.0)
        .collect();
    let cents = ratios(&heliox, &air, &predicted)
        .unwrap_or_else(|c| panic!("heliox/air peak ratios miss the image's by {c:?} cents"));
    println!("bank gas ratio: predicted {predicted:?}; realized against predicted {cents:?} cents");
    let negative = ratios(&air, &air, &predicted);
    assert!(
        negative.is_err(),
        "the air bank on both sides passes: {negative:?}"
    );
}
