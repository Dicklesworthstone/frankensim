//! Robustness of the live configuration (native):
//! 1. every key at pp..beyond-ff hammer speeds: contact refusals, finite
//!    output, energy books;
//! 2. a long random performance with pedal changes and more keys than the mode
//!    budget admits: forced retirements, refusals, energy books, finite output.
//! usage: stress [random_seconds]
use fs_live_piano::{Config, Description, LivePiano};

fn main() {
    let seconds: f64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(60.0);
    let rate = 48_000u32;
    let cfg = Config::live(rate);
    let desc = Description::embedded().expect("embedded description");

    // ---- 1. every key, five hammer speeds -------------------------------
    let mut p = LivePiano::new(&desc, cfg.clone()).expect("piano");
    let speeds = [0.2, 1.0, 2.5, 4.5, 8.0];
    let mut worst_residual = 0.0_f64;
    let mut nonfinite = 0u64;
    let mut errors: Vec<String> = Vec::new();
    let mut peak_by_speed = [0.0_f64; 5];
    let mut silent: Vec<(u8, f64)> = Vec::new();
    for midi in 21u8..=108 {
        for (si, &speed) in speeds.iter().enumerate() {
            p.note_on(midi, speed).expect("note on");
            let mut peak = 0.0_f64;
            for i in 0..9_600 {
                if i == 4_800 { p.note_off(midi); }
                let x = p.step();
                if !x.is_finite() { nonfinite += 1; }
                peak = peak.max(x.abs());
            }
            if let Some(e) = p.take_error() { errors.push(format!("key {midi} at {speed} m/s: {e}")); }
            peak_by_speed[si] = peak_by_speed[si].max(peak);
            if peak < 1.0e-6 { silent.push((midi, speed)); }
        }
        let l = p.close_ledger();
        worst_residual = worst_residual.max(l.residual_j.abs());
    }
    let l = p.close_ledger();
    let c = p.counters;
    println!("all-keys sweep: 88 keys x {:?} m/s; strikes {}; refused {}; contact failures {}; nonfinite samples {}; forced retirements {}; retirements {}; rebuilds {}",
        speeds, c.note_ons, c.refused_note_ons, c.contact_failures, nonfinite, c.forced_retirements, c.retirements, c.rebuilds);
    println!("  peak pressure by speed [Pa]: {:?}; overall peak {:.3} Pa; strikes with peak < 1e-6 Pa: {:?}", peak_by_speed.map(|x| (x * 1e3).round() / 1e3), c.peak_pa, silent);
    println!("  energy: input {:.6} J, stored {:.3e} J, dissipated {:.6} J, retired {:.3e} J, residual {:.3e} J (worst during sweep {:.3e} J)",
        l.input_j, l.stored_j, l.dissipated_j(), l.retired_j, l.residual_j, worst_residual);
    for e in errors.iter().take(12) { println!("  error: {e}"); }
    if errors.len() > 12 { println!("  ... {} more errors", errors.len() - 12); }

    // ---- 2. random performance -----------------------------------------
    let mut p = LivePiano::new(&desc, cfg.clone()).expect("piano");
    let mut seed = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move || {
        seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let steps = (seconds * f64::from(rate)) as usize;
    let mut held: Vec<u8> = Vec::new();
    let (mut nonfinite, mut max_keys, mut max_modes, mut peak) = (0u64, 0usize, 0usize, 0.0_f64);
    let mut errors: Vec<String> = Vec::new();
    let mut worst_residual = 0.0_f64;
    let start = std::time::Instant::now();
    let mut worst_block = 0.0_f64;
    let mut block_start = std::time::Instant::now();
    for i in 0..steps {
        // About 14 strikes per second on average, in bursts (chords), wide range.
        if next() < 14.0 / f64::from(rate) * 0.5 {
            let chord = 1 + (next() * 5.0) as usize;
            for _ in 0..chord {
                let midi = 24 + (next() * 78.0) as u8;
                let speed = 0.3 + next() * next() * 5.0;
                if p.note_on(midi, speed).is_ok() { held.push(midi); }
            }
        }
        if !held.is_empty() && next() < 10.0 / f64::from(rate) {
            let index = (next() * held.len() as f64) as usize;
            p.note_off(held.swap_remove(index));
        }
        if next() < 0.6 / f64::from(rate) {
            let travel = [0.0, 1.0, 0.5, 1.0][(next() * 4.0) as usize];
            p.set_sustain(travel).expect("pedal");
        }
        let x = p.step();
        if !x.is_finite() { nonfinite += 1; }
        peak = peak.max(x.abs());
        max_keys = max_keys.max(p.sounding_keys());
        max_modes = max_modes.max(p.dynamic_modes());
        if let Some(e) = p.take_error() { if errors.len() < 50 { errors.push(format!("sample {i}: {e}")); } }
        if i % 128 == 127 {
            worst_block = worst_block.max(block_start.elapsed().as_secs_f64());
            block_start = std::time::Instant::now();
        }
        if i % 48_000 == 0 { worst_residual = worst_residual.max(p.close_ledger().residual_j.abs()); }
    }
    let wall = start.elapsed().as_secs_f64();
    let l = p.close_ledger();
    let c = p.counters;
    println!("random performance: {seconds} s in {wall:.2} s (real-time factor {:.2}); worst 128-sample block {:.3} ms (2.667 ms of audio)", seconds / wall, worst_block * 1e3);
    println!("  strikes {}; refused {}; contact failures {}; forced retirements {}; retirements {}; rebuilds {}; nonfinite samples {}",
        c.note_ons, c.refused_note_ons, c.contact_failures, c.forced_retirements, c.retirements, c.rebuilds, nonfinite);
    println!("  max sounding keys {max_keys}; max dynamic modes {max_modes} (budget {}); peak {peak:.3} Pa; samples >= full scale ({} Pa): {}",
        cfg.max_dynamic_modes, cfg.full_scale_pa, c.clipped_samples);
    println!("  energy: input {:.6} J, stored {:.6} J, dissipated {:.6} J, retired {:.6e} J, residual {:.3e} J (worst sampled {:.3e} J)",
        l.input_j, l.stored_j, l.dissipated_j(), l.retired_j, l.residual_j, worst_residual);
    for e in errors.iter().take(12) { println!("  error: {e}"); }
    if errors.len() > 12 { println!("  ... {} more errors", errors.len() - 12); }
}
