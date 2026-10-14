//! Render one drum-kit piece to a 48 kHz float WAV.
//!
//! usage: kit PIECE OUT.wav [SPEED_M_S] [SECONDS]
//! PIECE: kick | snare | snare-off | tom-12 | tom-14 | hihat-closed | hihat-open | crash | ride
use fs_vkplate::kit::{self, Piece};
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        eprintln!("usage: kit PIECE OUT.wav [SPEED_M_S] [SECONDS]");
        std::process::exit(2);
    }
    let speed = a.get(3).map(|s| s.parse::<f64>().expect("speed")).filter(|v| *v > 0.0);
    let secs = a.get(4).map(|s| s.parse::<f64>().expect("seconds"));
    let v = |d: f64| speed.unwrap_or(d);
    // FS_KIT_TUNE scales a drum's head tensions (tuning experiments)
    let tune: f64 = std::env::var("FS_KIT_TUNE").ok().and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let d = |s: kit::DrumSpec| kit::retune(s, tune);
    let (piece, dur): (Result<Piece, String>, f64) = match a[1].as_str() {
        "kick" => (kit::drum(&d(kit::kick()), &[(0.005, v(3.0))]), 1.2),
        "snare" => (kit::drum(&d(kit::snare()), &[(0.005, v(4.5))]), 1.2),
        "snare-off" => (kit::drum(&d(kit::snare_wires_off()), &[(0.005, v(4.5))]), 1.2),
        "tom-12" => (kit::drum(&d(kit::tom(12.0)), &[(0.005, v(3.5))]), 2.0),
        "tom-14" => (kit::drum(&d(kit::tom(14.0)), &[(0.005, v(3.5))]), 2.5),
        "hihat-closed" => (kit::hihat(0.0, 80.0, &[(0.005, v(2.5))]), 0.8),
        "hihat-open" => (kit::hihat(0.008, 0.0, &[(0.005, v(2.5))]), 2.5),
        "crash" => (kit::cymbal(&kit::crash(), &[(0.005, v(6.0), 0.95, 0.008)]), 5.0),
        "ride" => (kit::cymbal(&kit::ride(), &[(0.005, v(2.5), 0.6, 0.004)]), 4.0),
        other => {
            eprintln!("unknown piece {other}");
            std::process::exit(2);
        }
    };
    let mut p = piece.unwrap_or_else(|e| {
        eprintln!("refused: {e}");
        std::process::exit(1)
    });
    let dur = secs.unwrap_or(dur);
    for n in &p.notes {
        eprintln!("{n}");
    }
    let (ratio, who) = p.instrument.stability();
    eprintln!("linear CFL use {ratio:.3} ({who}); internal rate {} kHz", p.instrument.rate / 1000.0);
    if ratio >= 0.95 {
        eprintln!("refused: leapfrog would be unstable");
        std::process::exit(1);
    }
    let t0 = Instant::now();
    let out = kit::render(&mut p, dur, |f| eprint!("\r{:5.1} %", 100.0 * f));
    eprintln!();
    let peak = out.iter().flatten().fold(0.0f64, |m, &x| m.max(x.abs()));
    let wall = t0.elapsed().as_secs_f64();
    for c in &p.instrument.contacts {
        if c.peak_force > 0.0 && c.a.len() + c.b.len() > 2 {
            eprintln!("strike peak force {:.1} N", c.peak_force);
        }
    }
    eprintln!("{dur} s rendered in {wall:.1} s ({:.1}x real time); peak {peak:.3} Pa", wall / dur);
    assert!(peak.is_finite(), "non-finite output");
    std::fs::write(&a[2], fs_vkplate::wav_f32(&out, kit::OUT_RATE as u32)).expect("write");
}
