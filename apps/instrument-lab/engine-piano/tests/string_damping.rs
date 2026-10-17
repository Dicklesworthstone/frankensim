//! Intrinsic string loss law (`Config::rt0425_string_damping`). The default stays the
//! reference's authored law, bit for bit. The published RT-0425 law, which the browser
//! plays, stops the bass's upper partials decaying ~5x faster than a real piano's.
use fs_live_piano::{Config, Description, LivePiano};

const RATE: u32 = 48_000;

fn render(rt0425: Option<bool>, midi: u8, seconds: f64) -> Vec<f64> {
    let mut cfg = Config::live(RATE);
    if let Some(on) = rt0425 {
        cfg.rt0425_string_damping = on;
    }
    let mut p = LivePiano::new(&Description::embedded().unwrap(), cfg).unwrap();
    p.note_on(midi, 4.4).unwrap();
    (0..(seconds * f64::from(RATE)) as usize).map(|_| p.step()).collect()
}

/// Energy of `x[start..start+len]` in [lo, hi) Hz (Hann window, direct DFT on 2 Hz bins).
fn band(x: &[f64], start: usize, len: usize, lo: f64, hi: f64) -> f64 {
    let w: Vec<f64> = (0..len)
        .map(|n| (0.5 - 0.5 * (std::f64::consts::TAU * n as f64 / len as f64).cos()) * x[start + n])
        .collect();
    let mut energy = 0.0;
    let mut f = lo;
    while f < hi {
        let step = std::f64::consts::TAU * f / f64::from(RATE);
        let (mut re, mut im) = (0.0, 0.0);
        for (n, v) in w.iter().enumerate() {
            re += v * (step * n as f64).cos();
            im -= v * (step * n as f64).sin();
        }
        energy += re * re + im * im;
        f += 2.0;
    }
    energy
}

/// dB lost by A0's 1.5-2.2 kHz band between the windows centred 0.3 s and 0.9 s.
fn high_band_drop(rt0425: bool) -> f64 {
    let x = render(Some(rt0425), 21, 1.0);
    let early = band(&x, 4_800 + 480, 9_600, 1_500.0, 2_200.0);
    let late = band(&x, 33_600 + 480, 9_600, 1_500.0, 2_200.0);
    10.0 * (early / late).log10()
}

#[test]
fn the_authored_law_remains_the_default_bit_for_bit() {
    let default = render(None, 45, 0.25);
    let authored = render(Some(false), 45, 0.25);
    assert!(default.iter().zip(&authored).all(|(a, b)| a.to_bits() == b.to_bits()));
    let published = render(Some(true), 45, 0.25);
    assert!(default.iter().zip(&published).any(|(a, b)| a != b), "the published law changed nothing");
}

#[test]
fn published_losses_let_the_bass_upper_partials_ring() {
    let authored = high_band_drop(false);
    let published = high_band_drop(true);
    println!("A0 1.5-2.2 kHz drop over 0.6 s: authored {authored:.1} dB, RT-0425 {published:.1} dB");
    // A real A0 (Salamander v16) loses ~8 dB here (12-13 dB/s).
    assert!(published > 0.5, "the upper partials must still decay ({published:.2} dB)");
    assert!(authored > 2.5 * published, "authored {authored:.1} dB vs RT-0425 {published:.1} dB");
}
