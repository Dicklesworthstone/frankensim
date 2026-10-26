//! Render one drum-kit piece to a 48 kHz float WAV.
//!
//! usage: kit PIECE OUT.wav [SPEED_M_S] [SECONDS] [--preset-heads]
//! PIECE: kick | snare | snare-off | tom-12 | tom-14 | hihat-closed | hihat-open | splash | crash | ride
//!
//! The drums' heads and stick come from sourced material cards, and no piece keeps a value fitted
//! to a recording (see `drum_cards`); `--preset-heads` renders the recording-fitted presets
//! instead, and says so.
use fs_vkplate::kit::{self, Piece};
use std::time::Instant;

mod drum_cards;

fn main() {
    let preset_heads = std::env::args().any(|a| a == drum_cards::PRESET_HEADS_FLAG);
    let a: Vec<String> = std::env::args().filter(|a| a != drum_cards::PRESET_HEADS_FLAG).collect();
    if a.len() < 3 {
        eprintln!("usage: kit PIECE OUT.wav [SPEED_M_S] [SECONDS] [--preset-heads]");
        std::process::exit(2);
    }
    let speed = a.get(3).map(|s| s.parse::<f64>().expect("speed")).filter(|v| *v > 0.0);
    let secs = a.get(4).map(|s| s.parse::<f64>().expect("seconds"));
    let v = |d: f64| speed.unwrap_or(d);
    // FS_KIT_TUNE scales a drum's head tensions (tuning experiments)
    let tune: f64 = std::env::var("FS_KIT_TUNE").ok().and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let drum = |s: kit::DrumSpec, kick: bool, v: f64| drum_cards::drum(&s, tune, &[(0.005, v)], kick, preset_heads);
    // FS_KIT_GRID scales a cymbal's grid spacing (refinement studies; 1 is the kit's own grid)
    let grid: f64 = std::env::var("FS_KIT_GRID").ok().and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let cy = |s: kit::CymbalSpec| kit::CymbalSpec { grid: s.grid * grid, ..s };
    let cymbal = |s: kit::CymbalSpec, hit: (f64, f64, f64, f64)| drum_cards::cymbal(&cy(s), &[hit], preset_heads);
    let hihat = |gap: f64, pedal_n: f64, v: f64| drum_cards::hihat(gap, pedal_n, &[(0.005, v)], preset_heads);
    let (piece, dur): (Result<Piece, String>, f64) = match a[1].as_str() {
        "kick" => (drum(kit::kick(), true, v(3.0)), 1.2),
        "snare" => (drum(kit::snare(), false, v(4.5)), 1.2),
        "snare-off" => (drum(kit::snare_wires_off(), false, v(4.5)), 1.2),
        "tom-12" => (drum(kit::tom(12.0), false, v(3.5)), 2.0),
        "tom-14" => (drum(kit::tom(14.0), false, v(3.5)), 2.5),
        "hihat-closed" => (hihat(0.0, 80.0, v(2.5)), 0.8),
        "hihat-open" => (drum_cards::open_hihat(&[(0.005, v(2.5))], preset_heads), 2.5),
        "splash" => (cymbal(kit::splash(), (0.005, v(2.0), 0.78, 0.005)), 2.6),
        "crash" => (cymbal(kit::crash(), (0.005, v(6.0), 0.95, 0.008)), 5.0),
        "ride" => (cymbal(drum_cards::ride(preset_heads), (0.005, v(2.5), 0.6, 0.004)), 4.0),
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
    // the striker's contact (a drum's "striker", a cymbal's or hi-hat's "stick"), not the wires or the rims; every
    // striker now lands in point contact, so a count of contact nodes no longer tells them apart
    let striker = |b: usize| matches!(p.instrument.bodies[b].name.as_str(), "striker" | "stick");
    for c in &p.instrument.contacts {
        if c.peak_force > 0.0 && c.b.iter().any(|&(b, _, _)| striker(b)) {
            eprintln!("strike peak force {:.1} N", c.peak_force);
        }
    }
    eprintln!("{dur} s rendered in {wall:.1} s ({:.1}x real time); peak {peak:.3} Pa", wall / dur);
    assert!(peak.is_finite(), "non-finite output");
    std::fs::write(&a[2], fs_vkplate::wav_f32(&out, kit::OUT_RATE as u32)).expect("write");
}
