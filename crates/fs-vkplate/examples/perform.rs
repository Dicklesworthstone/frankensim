//! Play one kit piece through a list of strokes: `perform PIECE HITS.csv OUT.wav SECONDS [--preset-heads]`.
//! HITS.csv holds `time_s,speed_m_s` lines (a header line is allowed); every stroke goes through
//! the same simulation, so each one lands on whatever the piece is still doing from the last.
//! The drums' heads and stick come from sourced material cards, and no piece keeps a value fitted
//! to a recording (see `drum_cards`); `--preset-heads` renders the recording-fitted presets
//! instead, and says so.
use fs_vkplate::kit::{self, Piece};
use std::time::Instant;

mod drum_cards;

fn main() {
    let preset_heads = std::env::args().any(|a| a == drum_cards::PRESET_HEADS_FLAG);
    let a: Vec<String> = std::env::args().filter(|a| a != drum_cards::PRESET_HEADS_FLAG).collect();
    if a.len() != 5 {
        eprintln!("usage: perform PIECE HITS.csv OUT.wav SECONDS [--preset-heads]");
        std::process::exit(2);
    }
    let text = std::fs::read_to_string(&a[2]).expect("read hits");
    let hits: Vec<(f64, f64)> = text
        .lines()
        .filter_map(|l| {
            let mut it = l.split(',').map(|s| s.trim().parse::<f64>());
            match (it.next(), it.next()) {
                (Some(Ok(t)), Some(Ok(v))) if t >= 0.0 && v > 0.0 => Some((t, v)),
                _ => None,
            }
        })
        .collect();
    let secs: f64 = a[4].parse().expect("seconds");
    let cym = |r: f64, tip: f64| hits.iter().map(|&(t, v)| (t, v, r, tip)).collect::<Vec<_>>();
    let piece: Result<Piece, String> = match a[1].as_str() {
        "kick" => drum_cards::drum(&kit::kick(), 1.0, &hits, true, preset_heads),
        "snare" => drum_cards::drum(&kit::snare(), 1.0, &hits, false, preset_heads),
        "snare-off" => drum_cards::drum(&kit::snare_wires_off(), 1.0, &hits, false, preset_heads),
        "tom-12" => drum_cards::drum(&kit::tom(12.0), 1.0, &hits, false, preset_heads),
        "tom-14" => drum_cards::drum(&kit::tom(14.0), 1.0, &hits, false, preset_heads),
        "hihat-closed" => drum_cards::hihat(0.0, 80.0, &hits, preset_heads),
        "hihat-open" => drum_cards::open_hihat(&hits, preset_heads),
        "crash" => drum_cards::cymbal(&kit::crash(), &cym(0.95, 0.008), preset_heads),
        "ride" => drum_cards::cymbal(&drum_cards::ride(preset_heads), &cym(0.6, 0.004), preset_heads),
        other => {
            eprintln!("unknown piece {other}");
            std::process::exit(2);
        }
    };
    let mut p = piece.unwrap_or_else(|e| {
        eprintln!("refused: {e}");
        std::process::exit(1)
    });
    for n in &p.notes {
        eprintln!("{n}");
    }
    let (ratio, who) = p.instrument.stability();
    if ratio >= 0.95 {
        eprintln!("refused: leapfrog would be unstable ({ratio:.3}, {who})");
        std::process::exit(1);
    }
    eprintln!("{} strokes over {secs} s", hits.len());
    let t0 = Instant::now();
    let out = kit::render(&mut p, secs, |f| eprint!("\r{:5.1} %", 100.0 * f));
    eprintln!();
    let peak = out.iter().flatten().fold(0.0f64, |m, &x| m.max(x.abs()));
    let wall = t0.elapsed().as_secs_f64();
    eprintln!("{secs} s rendered in {wall:.1} s ({:.1}x real time); peak {peak:.3} Pa", wall / secs);
    assert!(peak.is_finite(), "non-finite output");
    std::fs::write(&a[3], fs_vkplate::wav_f32(&out, kit::OUT_RATE as u32)).expect("write");
}
