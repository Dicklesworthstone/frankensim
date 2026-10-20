//! Native speed of the live kernel: real-time factor versus sounding keys.
//! usage: bench [board_band_hz] [course_mode_budget] [seconds] [--string-waves]
//! `--string-waves` stretches the plain steel strings with longitudinal waves
//! (`stretch::plain_steel_waves_spec`, 10 kHz axial band).
use fs_live_piano::{Config, Description, LivePiano};
use std::time::Instant;

const KEYS: [u8; 16] = [36, 43, 48, 52, 55, 60, 64, 67, 72, 76, 79, 84, 31, 40, 57, 62];

fn main() {
    let all: Vec<String> = std::env::args().collect();
    let waves = all.iter().any(|a| a == "--string-waves");
    if let Some(bad) = all.iter().skip(1).find(|a| a.starts_with("--") && *a != "--string-waves") {
        eprintln!("unknown option {bad}");
        std::process::exit(2);
    }
    let args: Vec<String> = all.into_iter().filter(|a| !a.starts_with("--")).collect();
    let band: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let budget: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(200);
    let seconds: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2.0);
    let rate = 48_000u32;
    let t = Instant::now();
    // band 0 = the embedded 32-coordinate board the browser uses.
    let desc = if band > 0.0 { Description::model_d(band) } else { Description::embedded() }.expect("description");
    let board_s = t.elapsed().as_secs_f64();
    let mut cfg = Config::live(rate);
    if band > 0.0 { cfg.board_band_hz = band; }
    cfg.course_mode_budget = budget;
    cfg.max_courses = 88;
    cfg.max_dynamic_modes = usize::MAX / 4;
    if waves {
        cfg.string_stretching = Some(fs_live_piano::stretch::plain_steel_waves_spec(&desc.courses, 10_000.0, 0.01, 0.05));
    }
    let t = Instant::now();
    let piano = LivePiano::new(&desc, cfg.clone()).expect("piano");
    println!("{}", fs_live_piano::describe(&piano));
    println!("cold: board {board_s:.3} s, tables+receiver {:.3} s", t.elapsed().as_secs_f64());
    drop(piano);
    let steps = (seconds * f64::from(rate)) as usize;
    println!("mode,keys,strings,dyn_modes,note_on_ms_max,wall_s,realtime_factor,us_per_step,peak_pa,input_mj,stored_mj,dissipated_mj,residual_j");
    for (label, pedal) in [("pedal", true), ("damped-release", false)] {
        for n in [0usize, 1, 2, 4, 6, 8, 10, 12, 16] {
            let mut p = LivePiano::new(&desc, cfg.clone()).expect("piano");
            p.set_sustain(if pedal { 1.0 } else { 0.0 }).unwrap();
            let mut on_max = 0.0_f64;
            for &k in &KEYS[..n] {
                let t = Instant::now();
                p.note_on(k, 2.5).expect("note on");
                on_max = on_max.max(t.elapsed().as_secs_f64());
            }
            let t = Instant::now();
            p.apply_pending_notes();
            on_max = on_max.max(t.elapsed().as_secs_f64());
            let (strings, modes) = (p.dynamic_strings(), p.dynamic_modes());
            let t = Instant::now();
            let mut peak = 0.0_f64;
            for i in 0..steps {
                if !pedal && i == steps / 4 {
                    for &k in &KEYS[..n] { p.note_off(k); }
                }
                peak = peak.max(p.step().abs());
            }
            let wall = t.elapsed().as_secs_f64();
            let l = p.close_ledger();
            println!("{label},{n},{strings},{modes},{:.3},{wall:.4},{:.2},{:.2},{peak:.4},{:.4},{:.4},{:.4},{:.3e}",
                on_max * 1e3, seconds / wall, wall / steps as f64 * 1e6, l.input_j * 1e3, l.stored_j * 1e3,
                l.dissipated_j() * 1e3, l.residual_j);
            if let Some(e) = p.take_error() { println!("  last error: {e}"); }
        }
    }
}
