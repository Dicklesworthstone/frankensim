//! Render one drum-kit piece to a 48 kHz float WAV.
//!
//! usage: kit PIECE OUT.wav [SPEED_M_S] [SECONDS] [--preset-heads] [--observer TABLE.tsv] [--port RADIUS_M [--port-flow-loss K]]
//! PIECE: kick | snare | snare-off | tom-12 | tom-14 | hihat-closed | hihat-open | splash | crash | crash-shaft | ride
//!
//! `--observer` hears a drum's overhead as one closed body on its floor (pexfp.378), from a table
//! baked for that drum by fs-bem's `examples/drum_observer.rs` (crates/fs-vkplate/data/closed-body);
//! it is refused for any other piece, and for a table baked for another drum, clock or gas.
//!
//! `--head-dipole`, with `--observer`, adds the heads' compact dipole radiation (`RadDipole`), its
//! length read from that table. It is an estimate; it is refused for a ported drum and on `--preset-heads`.
//! `--port RADIUS_M` declares a front-head opening; optional `--port-flow-loss K`
//! supplies its quasi-steady pressure-loss coefficient (no value is inferred).
//! `--port-radiation` replaces its frozen radiation resistance with the causal
//! uniform-piston response in an infinite rigid baffle, including exterior memory.
//! `--mic-first-order ALPHA FX FY FZ` declares an ideal pressure/velocity pattern
//! and its unit front axis (toward the source). ALPHA=0.5 is cardioid, 0 figure eight.
//! This is not measured capsule data; pressure-only closed-body tables refuse it.
//! `--head-mode-damping skrodzka2006`, on snare-off only, adds Skrodzka, Hojan and Proksza's
//! (2006, Table 1) measured damping ratios of another snare's batter modes to this batter's
//! matching shapes (`kit::SKRODZKA_2006_TABLE_1`). It is an upper bound: the measured ratios
//! already contain loss paths this model carries. It is refused on `--preset-heads`.
//! `--shell-loss b20`, on crash, ride, hihat-closed and hihat-open only (MB 12552), replaces each
//! shell's two declared T60 points with B20 bronze's printed loss at the same frequencies:
//! thermoelastic damping from fs-material's slab law plus an anelastic loss factor borrowed from
//! cast bells (`shell_loss`); radiation loss is missing. A hi-hat's two shells each take the loss
//! at their own free edge. `b20-low` and `b20-high` are the band's ends and never post. It is
//! refused on `--preset-heads`.
//! The card route's crash and ride are CuSn8 (`kit::CUSN8`, the recorded Paiste cymbals' alloy,
//! MB 13551), so a row is offered only on a shell of its own alloy: `b20` rows on the B20
//! hi-hats, `cusn8` rows on the crash and ride. `--shell-loss cusn8` is CuSn8's printed
//! thermoelastic damping alone; `cusn8-eta` adds the cast bells' borrowed anelastic factor;
//! `cusn8-low` is a band end (Kemper's lower conductivity) and never posts.
//! `--shell-radiation baffled-rayleigh`, on crash and ride with `--shell-loss` only (MB 14574 (d), a
//! correction), adds the radiation loss of the played shell's own rest modes
//! (`kit::attach_shell_radiation`): each mode up to the grid's resolution limit takes the baffled
//! Rayleigh resistance of its shape on the projected plane, from both faces, times the microphone's
//! own first-order dipole factor, and the modes of one sector whose overlap exceeds 0.05 share their
//! full resistance block (MB 14957); the modes above take none. Without a printed shell loss it is
//! refused (declared T60 points already lump radiation), and on the drums (their heads have their
//! own route), the hi-hats (two shells in contact) and `--preset-heads`.
//! `--kick-pillow pinta-fleece-foam02 --pillow-thickness M --pillow-backing rigid|air
//! [--pillow-radius M]`, on the kick only (MB 22635), replaces the declared muffle dashpot with the
//! pillow as a porous layer on each head (`fs_vkplate::pillow`): the layer's surface resistance from
//! the measured JCAL card at that thickness and backing, as a dashpot plus Maxwell branches fitted
//! to the card's curve alone, on the dashpot's patch (centre (0, -160) mm; radius 70 mm unless
//! `--pillow-radius`). The fleece stands in for pillow fill, a borrowed material; the thickness,
//! backing and patch are declared. It is refused on `--preset-heads`.
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
const PORT_FLAG: &str = "--port";
const PORT_FLOW_LOSS_FLAG: &str = "--port-flow-loss";
const PORT_RADIATION_FLAG: &str = "--port-radiation";
const MIC_FIRST_ORDER_FLAG: &str = "--mic-first-order";
/// Opt-in (MB 7247), snare-off only: `--head-mode-damping skrodzka2006`.
const HEAD_MODE_DAMPING_FLAG: &str = "--head-mode-damping";
/// Opt-in (MB 12552, 13551), cymbals only: `--shell-loss b20|b20-low|b20-high` on the hi-hats,
/// `--shell-loss cusn8|cusn8-eta|cusn8-low` on the crash and ride.
const SHELL_LOSS_FLAG: &str = "--shell-loss";
const SHELL_LOSS_ROWS: &str = "b20|b20-low|b20-high|cusn8|cusn8-eta|cusn8-low";
/// Opt-in (MB 14574 (d), a correction), crash and ride with `--shell-loss` only:
/// `--shell-radiation baffled-rayleigh`, the played shell's radiation loss on its own rest modes.
const SHELL_RADIATION_FLAG: &str = "--shell-radiation";
const SHELL_RADIATION_ROWS: &str = "baffled-rayleigh";
/// Opt-in (MB 22635), kick only: `--kick-pillow CARD --pillow-thickness M --pillow-backing rigid|air
/// [--pillow-radius M]`, the pillow as a porous layer on each head in place of the muffle dashpot.
const KICK_PILLOW_FLAG: &str = "--kick-pillow";
const PILLOW_THICKNESS_FLAG: &str = "--pillow-thickness";
const PILLOW_BACKING_FLAG: &str = "--pillow-backing";
const PILLOW_RADIUS_FLAG: &str = "--pillow-radius";

fn main() {
    let preset_heads = std::env::args().any(|a| a == drum_cards::PRESET_HEADS_FLAG);
    let all: Vec<String> = std::env::args().collect();
    let mic_index = all.iter().position(|a| a == MIC_FIRST_ORDER_FLAG);
    let mic_pattern = mic_index.map(|i| {
        let values: [f64; 4] = std::array::from_fn(|j| {
            all.get(i + 1 + j).and_then(|s| s.parse().ok()).filter(|x: &f64| x.is_finite()).unwrap_or_else(|| {
                eprintln!("usage: {MIC_FIRST_ORDER_FLAG} ALPHA FX FY FZ (finite values, explicit unit front axis)");
                std::process::exit(2)
            })
        });
        let axis = [values[1], values[2], values[3]];
        fs_vkplate::Mic::first_order([0.0; 3], values[0], axis).unwrap_or_else(|e| {
            eprintln!("refused: {e}");
            std::process::exit(1)
        });
        (values[0], axis)
    });
    let port_value = |flag: &str| all.iter().position(|a| a == flag).map(|i| {
        all.get(i + 1).and_then(|s| s.parse::<f64>().ok()).filter(|x| x.is_finite()).unwrap_or_else(|| {
            eprintln!("usage: {flag} FINITE_VALUE");
            std::process::exit(2)
        })
    });
    let port_radius = port_value(PORT_FLAG);
    let port_flow_loss = port_value(PORT_FLOW_LOSS_FLAG);
    let port_radiation = all.iter().any(|a| a == PORT_RADIATION_FLAG);
    if port_radius.is_some_and(|r| r <= 0.0) || port_flow_loss.is_some_and(|k| k < 0.0)
        || ((port_flow_loss.is_some() || port_radiation) && port_radius.is_none()) {
        eprintln!("refused: give {PORT_FLAG} RADIUS_M > 0, optional {PORT_FLOW_LOSS_FLAG} K >= 0 and optional {PORT_RADIATION_FLAG}");
        std::process::exit(1)
    }
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
    let head_mode_damping = all.iter().position(|a| a == HEAD_MODE_DAMPING_FLAG).map(|i| match all.get(i + 1).map(String::as_str) {
        Some("skrodzka2006") => (),
        _ => {
            eprintln!("usage: {HEAD_MODE_DAMPING_FLAG} skrodzka2006");
            std::process::exit(2);
        }
    });
    let shell_loss = all.iter().position(|a| a == SHELL_LOSS_FLAG).map(|i| {
        all.get(i + 1).and_then(|s| fs_vkplate::shell_loss::row(s)).unwrap_or_else(|| {
            eprintln!("usage: {SHELL_LOSS_FLAG} {SHELL_LOSS_ROWS}");
            std::process::exit(2)
        })
    });
    let shell_radiation = all.iter().position(|a| a == SHELL_RADIATION_FLAG).map(|i| match all.get(i + 1).map(String::as_str) {
        Some(SHELL_RADIATION_ROWS) => (),
        _ => {
            eprintln!("usage: {SHELL_RADIATION_FLAG} {SHELL_RADIATION_ROWS}");
            std::process::exit(2)
        }
    });
    let kick_pillow = all.iter().position(|a| a == KICK_PILLOW_FLAG).map(|i| {
        all.get(i + 1).and_then(|s| fs_vkplate::pillow::PillowCard::parse(s)).unwrap_or_else(|| {
            eprintln!("usage: {KICK_PILLOW_FLAG} {}", fs_vkplate::pillow::PillowCard::ROWS);
            std::process::exit(2)
        })
    });
    let pillow_thickness = port_value(PILLOW_THICKNESS_FLAG);
    let pillow_radius = port_value(PILLOW_RADIUS_FLAG);
    let pillow_backing = all.iter().position(|a| a == PILLOW_BACKING_FLAG).map(|i| {
        all.get(i + 1).and_then(|s| fs_vkplate::pillow::Backing::parse(s)).unwrap_or_else(|| {
            eprintln!("usage: {PILLOW_BACKING_FLAG} {}", fs_vkplate::pillow::Backing::ROWS);
            std::process::exit(2)
        })
    });
    if (pillow_thickness.is_some() || pillow_radius.is_some() || pillow_backing.is_some()) && kick_pillow.is_none() {
        eprintln!("refused: {PILLOW_THICKNESS_FLAG}, {PILLOW_BACKING_FLAG} and {PILLOW_RADIUS_FLAG} describe the pillow; give {KICK_PILLOW_FLAG} CARD too");
        std::process::exit(1)
    }
    if kick_pillow.is_some() && (pillow_thickness.is_none_or(|d| d <= 0.0) || pillow_backing.is_none() || pillow_radius.is_some_and(|r| r <= 0.0)) {
        eprintln!("refused: {KICK_PILLOW_FLAG} needs {PILLOW_THICKNESS_FLAG} M > 0 and {PILLOW_BACKING_FLAG} {} (optional {PILLOW_RADIUS_FLAG} M > 0)", fs_vkplate::pillow::Backing::ROWS);
        std::process::exit(1)
    }
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
    let valued = |i: usize| (i > 0 && (all[i - 1] == OBSERVER_FLAG || all[i - 1] == CAVITY_MODES_FLAG || all[i - 1] == WIRE_CONTACT_FLAG || all[i - 1] == PORT_FLAG || all[i - 1] == PORT_FLOW_LOSS_FLAG || all[i - 1] == HEAD_MODE_DAMPING_FLAG || all[i - 1] == SHELL_LOSS_FLAG || all[i - 1] == SHELL_RADIATION_FLAG || all[i - 1] == KICK_PILLOW_FLAG || all[i - 1] == PILLOW_THICKNESS_FLAG || all[i - 1] == PILLOW_BACKING_FLAG || all[i - 1] == PILLOW_RADIUS_FLAG))
        || mic_index.is_some_and(|j| i > j && i <= j + 4);
    let a: Vec<String> = all
        .iter()
        .enumerate()
        .filter(|(i, a)| {
            *a != drum_cards::PRESET_HEADS_FLAG && *a != OBSERVER_FLAG && *a != CAVITY_MODES_FLAG && *a != WIRE_CONTACT_FLAG && *a != CAVITY_WALL_LOSS_FLAG && *a != HEAD_DIPOLE_FLAG && *a != PORT_FLAG && *a != PORT_FLOW_LOSS_FLAG && *a != PORT_RADIATION_FLAG && *a != MIC_FIRST_ORDER_FLAG && *a != HEAD_MODE_DAMPING_FLAG && *a != SHELL_LOSS_FLAG && *a != SHELL_RADIATION_FLAG && *a != KICK_PILLOW_FLAG && *a != PILLOW_THICKNESS_FLAG && *a != PILLOW_BACKING_FLAG && *a != PILLOW_RADIUS_FLAG && !valued(*i)
        })
        .map(|(_, a)| a.clone())
        .collect();
    if a.len() < 3 {
        eprintln!("usage: kit PIECE OUT.wav [SPEED_M_S] [SECONDS] [--preset-heads] [--observer TABLE.tsv] [--port RADIUS_M [--port-flow-loss K] [--port-radiation]] [--mic-first-order ALPHA FX FY FZ] [--head-mode-damping skrodzka2006] [--shell-loss {SHELL_LOSS_ROWS} [--shell-radiation {SHELL_RADIATION_ROWS}]] [{KICK_PILLOW_FLAG} {} {PILLOW_THICKNESS_FLAG} M {PILLOW_BACKING_FLAG} {} [{PILLOW_RADIUS_FLAG} M]]", fs_vkplate::pillow::PillowCard::ROWS, fs_vkplate::pillow::Backing::ROWS);
        std::process::exit(2);
    }
    if kick_pillow.is_some() && (a[1] != "kick" || preset_heads) {
        eprintln!("refused: {KICK_PILLOW_FLAG} replaces the kick's declared muffle dashpot, offered on the kick without {} only", drum_cards::PRESET_HEADS_FLAG);
        std::process::exit(1)
    }
    if port_radius.is_some() && (!matches!(a[1].as_str(), "kick" | "snare" | "snare-off" | "tom-12" | "tom-14") || observer.is_some()) {
        eprintln!("refused: {PORT_FLAG} needs a drum and the vent-capable Rayleigh observer; closed-body tables contain no opening");
        std::process::exit(1)
    }
    if shell_loss.is_some() && (!matches!(a[1].as_str(), "crash" | "ride" | "hihat-closed" | "hihat-open") || preset_heads) {
        eprintln!("refused: {SHELL_LOSS_FLAG} replaces a cymbal shell's declared T60 pair, offered on crash, ride and the hi-hats without {} only", drum_cards::PRESET_HEADS_FLAG);
        std::process::exit(1)
    }
    if shell_radiation.is_some() && (shell_loss.is_none() || !matches!(a[1].as_str(), "crash" | "ride")) {
        eprintln!("refused: {SHELL_RADIATION_FLAG} adds the radiation a printed shell loss lacks, offered on crash and ride with {SHELL_LOSS_FLAG} only (declared T60 points already lump radiation; the drums' heads have their own route; the hi-hats are two shells in contact)");
        std::process::exit(1)
    }
    let speed = a.get(3).map(|s| s.parse::<f64>().expect("speed")).filter(|v| *v > 0.0);
    let secs = a.get(4).map(|s| s.parse::<f64>().expect("seconds"));
    let v = |d: f64| speed.unwrap_or(d);
    // FS_KIT_TUNE scales a drum's head tensions (tuning experiments)
    let tune: f64 = std::env::var("FS_KIT_TUNE").ok().and_then(|s| s.parse().ok()).unwrap_or(1.0);
    // FS_KIT_STROKE_TIMES lists a drum's strokes in seconds, comma separated, each at the same speed (a roll); unset, the one stroke at 5 ms
    let stroke_times: Vec<f64> = std::env::var("FS_KIT_STROKE_TIMES")
        .map_or_else(|_| vec![0.005], |s| s.split(',').map(|t| t.parse().expect("stroke time")).collect());
    // --kick-pillow: on the declared dashpot's patch (its radius unless --pillow-radius); card_drum then drops the dashpot
    let pillow = |s: &kit::DrumSpec| {
        kick_pillow.map(|card| {
            let Some((x, y, r, _)) = s.muffle else {
                eprintln!("refused: {KICK_PILLOW_FLAG} lies on the preset's muffle patch, and {} declares none", s.name);
                std::process::exit(1)
            };
            let (Some(thickness_m), Some(backing)) = (pillow_thickness, pillow_backing) else { unreachable!("checked above") };
            fs_vkplate::pillow::Pillow { card, thickness_m, backing, patch: (x, y, pillow_radius.unwrap_or(r)) }
        })
    };
    let configure_drum = |s: kit::DrumSpec| kit::DrumSpec {
        head_dipole,
        port: port_radius.or(s.port),
        port_flow_loss: port_flow_loss.unwrap_or(s.port_flow_loss),
        pillow: pillow(&s).or(s.pillow),
        ..s
    };
    let drum = |s: kit::DrumSpec, kick: bool, v: f64| {
        let hits: Vec<(f64, f64)> = stroke_times.iter().map(|&t| (t, v)).collect();
        drum_cards::drum(&configure_drum(s), tune, &hits, kick, preset_heads)
    };
    // FS_KIT_GRID scales a cymbal's grid spacing (refinement studies; 1 is the kit's own grid)
    let grid: f64 = std::env::var("FS_KIT_GRID").ok().and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let cy = |s: kit::CymbalSpec| kit::CymbalSpec { grid: s.grid * grid, ..s };
    let cymbal = |s: kit::CymbalSpec, hit: (f64, f64, f64, f64)| drum_cards::cymbal(&cy(s), &[hit], preset_heads);
    // --shell-loss: the played cymbal's shell with its T60 pair from print (fs_vkplate::shell_loss), refused on a
    // shell of another alloy than the row's; every other shell, and every shell without the flag, is returned as it came
    let played = a[1].as_str();
    let shell = |s: kit::CymbalSpec, on: bool| match &shell_loss {
        Some(loss) if on => fs_vkplate::shell_loss::with_printed_shell_loss(&s, loss)
            .map(|(spec, line)| (spec, Some(line)))
            .unwrap_or_else(|e| {
                eprintln!("refused: {e}");
                std::process::exit(1)
            }),
        _ => (s, None),
    };
    let (crash_shell, crash_line) = shell(drum_cards::crash(preset_heads), played == "crash");
    let (ride_shell, ride_line) = shell(drum_cards::ride(preset_heads), played == "ride");
    let hihat = |gap: f64, pedal_n: f64, v: f64| drum_cards::hihat(gap, pedal_n, &[(0.005, v)], preset_heads);
    // --shell-loss on a hi-hat: each of its two card shells takes the loss at its own free edge, with its own receipt
    let (hat_shells, hat_lines) = {
        let on = matches!(played, "hihat-closed" | "hihat-open");
        let (top, bottom) = kit::card_hihat_spec();
        let (top, top_line) = shell(top, on);
        let (bottom, bottom_line) = shell(bottom, on);
        ((top, bottom), top_line.into_iter().chain(bottom_line).collect::<Vec<String>>())
    };
    let (piece, dur): (Result<Piece, String>, f64) = match a[1].as_str() {
        "kick" => (drum(kit::kick(), true, v(3.0)), 1.2),
        "snare" => (drum(kit::snare(), false, v(4.5)), 1.2),
        "snare-off" => (drum(kit::snare_wires_off(), false, v(4.5)), 1.2),
        "tom-12" => (drum(kit::tom(12.0), false, v(3.5)), 2.0),
        "tom-14" => (drum(kit::tom(14.0), false, v(3.5)), 2.5),
        "hihat-closed" if shell_loss.is_some() => (drum_cards::hihat_with_shells(0.0, 80.0, &[(0.005, v(2.5))], hat_shells), 0.8),
        "hihat-closed" => (hihat(0.0, 80.0, v(2.5)), 0.8),
        "hihat-open" if shell_loss.is_some() => (drum_cards::open_hihat_with_shells(&[(0.005, v(2.5))], hat_shells), 2.5),
        "hihat-open" => (drum_cards::open_hihat(&[(0.005, v(2.5))], preset_heads), 2.5),
        "splash" => (cymbal(kit::splash(), (0.005, v(2.0), 0.78, 0.005)), 2.6),
        "crash" => (cymbal(crash_shell, (0.005, v(6.0), 0.95, 0.008)), 5.0),
        "crash-shaft" => (drum_cards::shaft_cymbal(&cy(kit::card_crash()), speed), 5.0),
        "ride" => (cymbal(ride_shell, (0.005, v(2.5), 0.6, 0.004)), 4.0),
        other => {
            eprintln!("unknown piece {other}");
            std::process::exit(2);
        }
    };
    let mut p = piece.unwrap_or_else(|e| {
        eprintln!("refused: {e}");
        std::process::exit(1)
    });
    if shell_loss.is_some() {
        match a[1].as_str() {
            "crash" => p.notes.extend(crash_line),
            "ride" => p.notes.extend(ride_line),
            _ => p.notes.extend(hat_lines),
        }
    }
    // --shell-radiation: the played shell, rebuilt as the match built it (card route, FS_KIT_GRID applied)
    let mut radiating: Option<&'static str> = None;
    if shell_radiation.is_some() {
        let spec = cy(if played == "crash" { drum_cards::crash(false) } else { drum_cards::ride(false) });
        match kit::attach_shell_radiation(&mut p, &spec) {
            Ok(lines) => p.notes.extend(lines),
            Err(e) => {
                eprintln!("refused: {e}");
                std::process::exit(1)
            }
        }
        radiating = Some(spec.name);
    }
    if let Some((alpha, axis)) = mic_pattern {
        for mic in 0..p.instrument.mics.len() {
            p.instrument.set_microphone_pattern(mic, alpha, axis).unwrap_or_else(|e| {
                eprintln!("refused: {e}");
                std::process::exit(1)
            });
        }
        p.notes.push(format!("ideal pressure/particle-velocity microphone: alpha={alpha}, front={axis:?}; finite-aperture and near-field response, no measured capsule/electronics model"));
    }
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
        match kit::attach_cavity_modes_with(&mut p, &configure_drum(spec), f_max, wall_loss) {
            Ok(receipt) => p.notes.push(receipt),
            Err(e) => {
                eprintln!("refused: {e}");
                std::process::exit(1)
            }
        }
    }
    if port_radiation {
        match kit::attach_baffled_vent_radiation(&mut p) {
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
    if head_mode_damping.is_some() {
        if a[1] != "snare-off" || preset_heads {
            eprintln!("refused: {HEAD_MODE_DAMPING_FLAG} skrodzka2006 is a snare batter's measured damping, offered on snare-off's card heads only");
            std::process::exit(1)
        }
        match kit::attach_skrodzka_2006_batter_damping(&mut p, &kit::snare_wires_off()) {
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
        for o in m.closed_bodies() {
            eprintln!("{}", o.report());
        }
    }
    if let Some(dipole) = &p.instrument.rad_dipole {
        eprintln!("head dipole radiated {:.4e} J", dipole.radiated_j);
    }
    for (index, cavity) in p.instrument.cavities.iter().enumerate() {
        if let Some(port) = &cavity.port {
            eprintln!("vent {index}: {:.4e} J removed by linear and separated-flow pressure loss (K={})",
                port.dissipated_j, port.flow_loss_coefficient());
            if port.has_baffled_radiation() {
                eprintln!("vent {index}: {:.4e} J radiated, {:.4e} J stored in the exterior radiation memory",
                    port.radiated_j(), port.radiation_energy_j());
            }
        }
    }
    if let Some(body) = radiating.and_then(|name| p.instrument.bodies.iter().find(|b| b.name == name)) {
        eprintln!("shell radiation removed {:.4e} J", body.shape_damping_dissipated_j());
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
