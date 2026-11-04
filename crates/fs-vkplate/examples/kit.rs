//! Render one drum-kit piece to a 48 kHz float WAV.
//!
//! usage: kit PIECE OUT.wav [SPEED_M_S] [SECONDS] [--preset-heads] [--observer TABLE.tsv]
//! PIECE: kick | snare | snare-off | tom-12 | tom-14 | hihat-closed | hihat-open | splash | crash | crash-shaft | ride
//!
//! `--observer` hears a drum's overhead as one closed body on its floor (pexfp.378), from a table
//! baked for that drum by fs-bem's `examples/drum_observer.rs` (crates/fs-vkplate/data/closed-body);
//! it is refused for any other piece, and for a table baked for another drum, clock or gas.
//!
//! `--head-dipole`, with `--observer`, adds the heads' compact dipole radiation (`RadDipole`), its
//! length read from that table. It is an estimate; it is refused for a ported drum and on `--preset-heads`.
//!
//! The drums' heads and stick come from sourced material cards, and no piece keeps a value fitted
//! to a recording (see `drum_cards`); `--preset-heads` renders the recording-fitted presets
//! instead, and says so. `crash-shaft` is the crash struck at its edge with the hickory card's
//! shaft, the wood crushing across its grain; it has no default speed (render a ladder of
//! SPEED_M_S).
use fs_vkplate::Law;
use fs_vkplate::kit::{self, Piece};
use std::time::Instant;

mod drum_cards;

const OBSERVER_FLAG: &str = "--observer";
/// ALL-IN (MB 14741), opt-in: the drum's enclosed air as rigid-cylinder acoustic modes below F_MAX_HZ.
const CAVITY_MODES_FLAG: &str = "--cavity-modes";
/// ALL-IN (MB 14741), opt-in: `--wire-contact viscoelastic`, the snare wires' contact loss from the film card (.372).
const WIRE_CONTACT_FLAG: &str = "--wire-contact";
/// ALL-IN (MB 15788), opt-in with --cavity-modes: the air modes' viscothermal wall loss (Pierce Eq. (10-4.16)).
const CAVITY_WALL_LOSS_FLAG: &str = "--cavity-wall-loss";
/// ALL-IN, opt-in with --observer: the heads' compact dipole radiation, its length from that table.
const HEAD_DIPOLE_FLAG: &str = "--head-dipole";

fn main() {
    let preset_heads = std::env::args().any(|a| a == drum_cards::PRESET_HEADS_FLAG);
    let all: Vec<String> = std::env::args().collect();
    let observer = all.iter().position(|a| a == OBSERVER_FLAG).map(|i| {
        all.get(i + 1).cloned().unwrap_or_else(|| {
            eprintln!("usage: {OBSERVER_FLAG} TABLE.tsv");
            std::process::exit(2)
        })
    });
    let cavity_modes = all.iter().position(|a| a == CAVITY_MODES_FLAG).map(|i| {
        all.get(i + 1).and_then(|s| s.parse::<f64>().ok()).unwrap_or_else(|| {
            eprintln!("usage: {CAVITY_MODES_FLAG} F_MAX_HZ");
            std::process::exit(2)
        })
    });
    let wire_contact = all.iter().position(|a| a == WIRE_CONTACT_FLAG).map(|i| match all.get(i + 1).map(String::as_str) {
        Some("viscoelastic") => (),
        _ => {
            eprintln!("usage: {WIRE_CONTACT_FLAG} viscoelastic");
            std::process::exit(2)
        }
    });
    let wall_loss = all.iter().any(|a| a == CAVITY_WALL_LOSS_FLAG);
    if wall_loss && cavity_modes.is_none() {
        eprintln!("refused: {CAVITY_WALL_LOSS_FLAG} is the loss of the air modes; give {CAVITY_MODES_FLAG} F_MAX_HZ too");
        std::process::exit(1)
    }
    let head_dipole = all.iter().any(|a| a == HEAD_DIPOLE_FLAG).then(|| {
        let Some(path) = &observer else {
            eprintln!("refused: {HEAD_DIPOLE_FLAG} takes its length from the piece's closed-body table; give {OBSERVER_FLAG} TABLE.tsv too");
            std::process::exit(1)
        };
        std::fs::read_to_string(path)
            .map_err(|e| format!("cannot read {path}: {e}"))
            .and_then(|text| kit::head_dipole_from_table(&text, cavity_modes))
            .unwrap_or_else(|e| {
                eprintln!("refused: {e}");
                std::process::exit(1)
            })
    });
    let valued = |i: usize| i > 0 && (all[i - 1] == OBSERVER_FLAG || all[i - 1] == CAVITY_MODES_FLAG || all[i - 1] == WIRE_CONTACT_FLAG);
    let a: Vec<String> = all
        .iter()
        .enumerate()
        .filter(|(i, a)| {
            *a != drum_cards::PRESET_HEADS_FLAG && *a != OBSERVER_FLAG && *a != CAVITY_MODES_FLAG && *a != WIRE_CONTACT_FLAG && *a != CAVITY_WALL_LOSS_FLAG && *a != HEAD_DIPOLE_FLAG && !valued(*i)
        })
        .map(|(_, a)| a.clone())
        .collect();
    if a.len() < 3 {
        eprintln!("usage: kit PIECE OUT.wav [SPEED_M_S] [SECONDS] [--preset-heads] [--observer TABLE.tsv]");
        std::process::exit(2);
    }
    let speed = a.get(3).map(|s| s.parse::<f64>().expect("speed")).filter(|v| *v > 0.0);
    let secs = a.get(4).map(|s| s.parse::<f64>().expect("seconds"));
    let v = |d: f64| speed.unwrap_or(d);
    // FS_KIT_TUNE scales a drum's head tensions (tuning experiments)
    let tune: f64 = std::env::var("FS_KIT_TUNE").ok().and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let drum = |s: kit::DrumSpec, kick: bool, v: f64| drum_cards::drum(&kit::DrumSpec { head_dipole, ..s }, tune, &[(0.005, v)], kick, preset_heads);
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
        "crash-shaft" => (drum_cards::shaft_cymbal(&cy(kit::crash()), speed), 5.0),
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
    if let Some(path) = &observer {
        let spec = match a[1].as_str() {
            "kick" => kit::kick(),
            "snare" => kit::snare(),
            "snare-off" => kit::snare_wires_off(),
            "tom-12" => kit::tom(12.0),
            "tom-14" => kit::tom(14.0),
            other => {
                eprintln!("refused: {OBSERVER_FLAG} hears a drum as one closed body; {other} is not a drum");
                std::process::exit(1)
            }
        };
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
            eprintln!("refused: cannot read {path}: {e}");
            std::process::exit(1)
        });
        if let Err(e) = kit::attach_closed_body(&mut p, &spec, &text) {
            eprintln!("refused: {e}");
            std::process::exit(1)
        }
    }
    if let Some(f_max) = cavity_modes {
        let spec = match a[1].as_str() {
            "kick" => kit::kick(),
            "snare" => kit::snare(),
            "snare-off" => kit::snare_wires_off(),
            "tom-12" => kit::tom(12.0),
            "tom-14" => kit::tom(14.0),
            other => {
                eprintln!("refused: {CAVITY_MODES_FLAG} is a drum's enclosed air; {other} is not a drum");
                std::process::exit(1)
            }
        };
        match kit::attach_cavity_modes_with(&mut p, &spec, f_max, wall_loss) {
            Ok(receipt) => p.notes.push(receipt),
            Err(e) => {
                eprintln!("refused: {e}");
                std::process::exit(1)
            }
        }
    }
    if wire_contact.is_some() {
        if a[1] != "snare" {
            eprintln!("refused: {WIRE_CONTACT_FLAG} is the snare wires' contact; {} has none", a[1]);
            std::process::exit(1)
        }
        match kit::attach_viscoelastic_wires(&mut p, &kit::snare()) {
            Ok(receipt) => p.notes.push(receipt),
            Err(e) => {
                eprintln!("refused: {e}");
                std::process::exit(1)
            }
        }
    }
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
            if let Some(crush) = &c.crush {
                let Law::Power { k, .. } = c.law;
                eprintln!(
                    "  crushed to {:.1} um (crushing begins at {:.2} um), permanent set {:.1} um",
                    crush.deepest * 1e6,
                    crush.onset(k) * 1e6,
                    crush.set(k) * 1e6
                );
            }
        }
    }
    for m in &p.instrument.mics {
        if let Some(o) = m.closed_body() {
            eprintln!("{}", o.report());
        }
    }
    if let Some(dipole) = &p.instrument.rad_dipole {
        eprintln!("head dipole radiated {:.4e} J", dipole.radiated_j);
    }
    eprintln!("{dur} s rendered in {wall:.1} s ({:.1}x real time); peak {peak:.3} Pa", wall / dur);
    assert!(peak.is_finite(), "non-finite output");
    std::fs::write(&a[2], fs_vkplate::wav_f32(&out, kit::OUT_RATE as u32)).expect("write");
}

#[cfg(test)]
#[path = "kit_tests/identity.rs"]
mod identity_tests;
#[cfg(test)]
#[path = "kit_tests/closed_body_tables.rs"]
mod closed_body_tables_tests;
