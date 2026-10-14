//! Native measurements and reference comparison for the live drum.
//!
//!   bench info   [params_json]                      frequencies, build time
//!   bench trace  <steps> <speed> <x> <y>            500 kHz CSV in the reference's columns
//!   bench rt     [params_json]                      real-time factor per substeps
//!   bench conv   [params_json]                      pressure deviation per substeps vs 16
//!   bench rearm  [params_json]                      launch/catch recompile cost
//!   bench wav    <out.f32> <seconds> [params_json]  raw float32 mono of the batter microphone
//!
//! Wall-clock timing lives ONLY here (std::time is not available in the wasm path).
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_exec::CancelGate;
use fs_live_drum::engine::{Voice, compile_system};
use fs_live_drum::instrument::{Instrument, LAUNCH_GAP_M};
use fs_live_drum::{DrumParams, LiveDrumCore};
use std::io::Write;
use std::time::Instant;

type Error = Box<dyn std::error::Error>;

fn params(arg: Option<&String>) -> Result<DrumParams, Error> {
    Ok(DrumParams::default().overlay_json(arg.map_or("{}", String::as_str))?)
}

fn info(arg: Option<&String>) -> Result<(), Error> {
    let p = params(arg)?;
    let start = Instant::now();
    let inst = Instrument::build(p)?;
    let build = start.elapsed();
    let start = Instant::now();
    let voice = Voice::new(&inst, p.substeps)?;
    let compile = start.elapsed();
    for (i, head) in inst.heads.iter().enumerate() {
        println!("head {i}: modes={} film_mass_kg={} frequencies_hz={:?} max_residual={:e}",
            head.mode_count(), head.film.mass_kg, head.frequencies_hz(), head.max_residual);
    }
    println!("coupled_frequencies_hz={:?}", inst.coupled_frequencies_hz);
    println!("mesh nodes={} triangles={}", inst.nodes.len(), inst.tris.len());
    println!("stick full_mass_kg={:.8} effective_mass_kg={:.8} hertz_K={:e}",
        inst.stick.full_mass_kg, inst.stick.effective_mass_kg, inst.stick.hertz_stiffness);
    println!("build_ms={:.3} voice_compile_ms={:.3}", build.as_secs_f64() * 1e3, compile.as_secs_f64() * 1e3);
    println!("filter_delay_frames={:?} travel_delay_frames={:?} radiation_multiply_adds={}",
        voice.filter_delay_frames(), voice.travel_delay_frames(), voice.radiation_multiply_adds());
    Ok(())
}

/// The reference's mechanics CSV (main.rs run_args), at its 500 kHz clock.
fn trace(args: &[String]) -> Result<(), Error> {
    let steps: u64 = args.first().ok_or("steps")?.parse()?;
    let speed: f64 = args.get(1).ok_or("speed")?.parse()?;
    let x: f64 = args.get(2).ok_or("x")?.parse()?;
    let y: f64 = args.get(3).ok_or("y")?.parse()?;
    let inst = Instrument::build(DrumParams::default())?;
    let gate = CancelGate::new_clock_free();
    let (n0, n1) = (inst.heads[0].mode_count(), inst.heads[1].mode_count());
    let row = inst.batter_shapes_at([x, y])?;
    let root = inst.stick.effective_mass_kg.sqrt();
    let tip = ModalAcousticState {
        displacement_m_sqrt_kg: -LAUNCH_GAP_M * root,
        velocity_m_sqrt_kg_per_s: speed * root,
    };
    let rest = vec![ModalAcousticState::default(); n0 + n1];
    let mut system = compile_system(&inst, 500_000, &rest, Some((tip, &row)), &gate)?;
    let bottom = {
        let head = &inst.heads[1];
        fs_couple::render::plate::impact::linear::wire::film_shapes(&head.film, &head.modes, &[[x, y]])?.remove(0)
    };
    let areas = inst.compression_areas();
    let k = inst.params.bulk_modulus_pa() / inst.params.volume_m3();
    let force = vec![0.0; system.mode_count()];
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    writeln!(out, "time_s,point_a_displacement_m,point_a_velocity_m_s,point_b_displacement_m,cavity_internal_pa,total_energy_j,felt_crush_j,loss_j,balance_j")?;
    for _ in 0..steps {
        let f = system.step(&force, &gate)?;
        let s = system.state();
        let a: f64 = row.iter().enumerate().map(|(i, b)| b * s[2 * (1 + i)]).sum();
        let av: f64 = row.iter().enumerate().map(|(i, b)| b * s[2 * (1 + i) + 1]).sum();
        let b: f64 = bottom.iter().enumerate().map(|(i, b)| b * s[2 * (1 + n0 + i)]).sum();
        let pressure = k * areas.iter().enumerate().map(|(i, a)| a * s[2 * (1 + i)]).sum::<f64>();
        writeln!(out, "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            f.time_s, a, av, b, pressure, f.stored_energy_j, 0.0, f.dissipated_energy_j, f.balance_residual_j)?;
    }
    out.flush()?;
    Ok(())
}

fn rt(arg: Option<&String>) -> Result<(), Error> {
    let p = params(arg)?;
    let inst = Instrument::build(p)?;
    println!("modes: batter={} resonant={}", inst.heads[0].mode_count(), inst.heads[1].mode_count());
    println!("substeps,scenario,audio_s,wall_s,realtime_factor,steps_per_s,ns_per_step,stick_step_fraction,contact_step_fraction,peak_pa,balance_residual_j");
    for substeps in [16, 8, 4, 2, 1] {
        if p.band_high_hz / (48_000.0 * substeps as f64) >= 0.45 { continue; }
        // one stroke then ring; a 20 Hz roll; a stroke every 5 ms (stick ~always present)
        for (name, frames, interval) in [("single", 48_000, 0), ("roll20hz", 48_000, 2_400), ("dense5ms", 24_000, 240)] {
            let start = Instant::now();
            let r = Voice::bench(&inst, substeps, frames, interval, [0.06, 0.01], 2.0)?;
            let wall = start.elapsed().as_secs_f64();
            let audio = frames as f64 / 48_000.0;
            println!("{substeps},{name},{audio:.3},{wall:.4},{:.3},{:.0},{:.1},{:.4},{:.4},{:.4},{:e}",
                audio / wall, r.mechanical_steps as f64 / wall, wall * 1e9 / r.mechanical_steps as f64,
                r.steps_with_stick as f64 / r.mechanical_steps as f64,
                r.steps_in_contact as f64 / r.mechanical_steps as f64, r.peak_pa, r.ledger.balance_residual_j);
        }
    }
    Ok(())
}

fn render(p: DrumParams, frames: usize, strikes: &[(usize, [f64; 2], f64)]) -> Result<(Vec<f64>, Vec<f64>, fs_live_drum::Ledger), Error> {
    let mut drum = LiveDrumCore::new(&p.to_json())?;
    let mut a = Vec::with_capacity(frames);
    let mut b = Vec::with_capacity(frames);
    for frame in 0..frames {
        for (at, xy, speed) in strikes {
            if *at == frame { drum.strike(xy[0], xy[1], *speed)?; }
        }
        let (x, y) = drum.render_pa(1);
        a.push(x[0]);
        b.push(y[0]);
        if let Some(r) = drum.take_refusal() { return Err(r.into()); }
    }
    let ledger = drum.snapshot().ledger;
    Ok((a, b, ledger))
}

fn rms(x: &[f64]) -> f64 { (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt() }

fn conv(arg: Option<&String>) -> Result<(), Error> {
    let base = params(arg)?;
    let frames = 24_000;
    let strikes = [(0usize, [0.06, 0.01], 2.0)];
    let reference = render(DrumParams { substeps: 16, ..base }, frames, &strikes)?;
    println!("substeps,mechanical_rate_hz,peak_pa_batter,peak_pa_resonant,rms_pa_batter,relative_rms_deviation_vs_16,max_abs_deviation_pa,balance_residual_j");
    for substeps in [16, 8, 4, 2, 1] {
        if base.band_high_hz / (48_000.0 * substeps as f64) >= 0.45 { continue; }
        let (a, b, ledger) = render(DrumParams { substeps, ..base }, frames, &strikes)?;
        // Decimator latency differs by a few frames between ratios (44..47
        // frames, 0 at ratio 1): align by the declared filter delay before
        // differencing so the number measures mechanics, not a known lag.
        let delay = |s: usize| -> usize { match s { 1 => 0, 2 => 40, 4 => 44, 8 => 46, _ => 47 } };
        let (shift_ref, shift) = (delay(16), delay(substeps));
        let n = frames - 64;
        let difference: Vec<f64> = (0..n).map(|i| a[i + shift] - reference.0[i + shift_ref]).collect();
        let peak = |x: &[f64]| x.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        println!("{substeps},{},{:.5},{:.5},{:.5},{:.3e},{:.3e},{:.3e}", 48_000 * substeps, peak(&a), peak(&b), rms(&a),
            rms(&difference) / rms(&reference.0[shift_ref..shift_ref + n]), peak(&difference), ledger.balance_residual_j);
    }
    Ok(())
}

fn rearm(arg: Option<&String>) -> Result<(), Error> {
    let p = params(arg)?;
    let inst = Instrument::build(p)?;
    let gate = CancelGate::new_clock_free();
    let (n0, n1) = (inst.heads[0].mode_count(), inst.heads[1].mode_count());
    let row = inst.batter_shapes_at([0.06, 0.01])?;
    let rest = vec![ModalAcousticState::default(); n0 + n1];
    let root = inst.stick.effective_mass_kg.sqrt();
    let tip = ModalAcousticState { displacement_m_sqrt_kg: -LAUNCH_GAP_M * root, velocity_m_sqrt_kg_per_s: 2.0 * root };
    let rate = p.mechanical_rate_hz();
    let n = 2000;
    let start = Instant::now();
    for _ in 0..n { std::hint::black_box(compile_system(&inst, rate, &rest, Some((tip, &row)), &gate)?); }
    let with = start.elapsed().as_secs_f64() / f64::from(n);
    let start = Instant::now();
    for _ in 0..n { std::hint::black_box(compile_system(&inst, rate, &rest, None, &gate)?); }
    let without = start.elapsed().as_secs_f64() / f64::from(n);
    let start = Instant::now();
    for _ in 0..n { std::hint::black_box(inst.batter_shapes_at([0.06, 0.01])?); }
    let shapes = start.elapsed().as_secs_f64() / f64::from(n);
    println!("modes={} launch_recompile_us={:.2} catch_recompile_us={:.2} strike_row_us={:.2}",
        n0 + n1, with * 1e6, without * 1e6, shapes * 1e6);
    Ok(())
}

fn wav(args: &[String]) -> Result<(), Error> {
    let path = args.first().ok_or("output path")?;
    let seconds: f64 = args.get(1).ok_or("seconds")?.parse()?;
    let p = params(args.get(2))?;
    let frames = (seconds * 48_000.0) as usize;
    let (a, _, ledger) = render(p, frames, &[(0, [0.06, 0.01], 2.0)])?;
    let mut file = std::fs::File::create(path)?;
    for v in &a { file.write_all(&((v / p.full_scale_pa) as f32).to_le_bytes())?; }
    eprintln!("wrote {frames} float32 frames; peak_pa={:.4}; ledger={ledger:?}", a.iter().fold(0.0_f64, |m, v| m.max(v.abs())));
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("info") => info(args.get(1)),
        Some("trace") => trace(&args[1..]),
        Some("rt") => rt(args.get(1)),
        Some("conv") => conv(args.get(1)),
        Some("rearm") => rearm(args.get(1)),
        Some("wav") => wav(&args[1..]),
        _ => Err("usage: bench info|trace|rt|conv|rearm|wav (see the file header)".into()),
    };
    if let Err(e) = result {
        eprintln!("bench refused: {e}");
        std::process::exit(1);
    }
}
