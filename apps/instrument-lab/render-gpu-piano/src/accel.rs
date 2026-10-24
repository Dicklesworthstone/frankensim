//! The grand_piano `Accelerator`: contact-free stretches on a backend, hammers on the CPU.
//!
//! A stretch runs from the present sample to the next control event (at most `SEGMENT` steps), as long
//! as the planner (`Instrument::contact_free_steps`) sees no hammer near its string. The backend steps
//! the bank and reports every step's contact positions; the CPU then replays the hammers' free flight
//! against those positions in complete output frames. A refusal restores the whole frame: the backend
//! reruns only the preceding complete frames, and the refused frame goes to the CPU. Backend limits,
//! traces and statistics count mechanics steps; `Accelerator::advance` takes and returns output frames.
//! Each stretch must close the instrument's energy ledger.
use crate::emulate::Emulator;
use crate::grand_piano::audio::Accelerator;
use crate::grand_piano::engine::Instrument;
use crate::plan::Program;
use std::time::Instant;

/// At most this many steps per stretch (the traces of one stretch are held in memory).
pub const SEGMENT: usize = 2048;
/// Shorter stretches go to the CPU.
const MIN_SEGMENT: usize = 16;
/// The planner keeps each hammer this far below its string's amplitude bound (times the margin).
const PLAN_CLEARANCE_M: f64 = 2.0e-4;
const PLAN_MARGIN: f64 = 1.5;
/// The replay refuses a step whose hammer comes this close to the string's actual position.
const EXACT_CLEARANCE_M: f64 = 1.0e-6;

pub trait Backend {
    fn name(&self) -> String;
    /// Launches rerun because the GPU did not schedule every block (Metal; see metal.rs).
    fn stalls(&self) -> usize {
        0
    }
    /// Load the program's damper operators (after `Program::set_dampers`).
    fn load_dampers(&mut self, p: &Program) -> Result<(), String>;
    /// Step (q, v) `steps` times; per step the board velocity into `trace` (r values) and the contact
    /// positions into `positions` (one per contact); returns (modal loss, damper loss).
    fn run(
        &mut self,
        p: &Program,
        q: &mut [f64],
        v: &mut [f64],
        steps: usize,
        trace: &mut [f64],
        positions: &mut [f64],
    ) -> Result<(f64, f64), String>;
}

pub struct Cpu(Emulator);
impl Backend for Cpu {
    fn name(&self) -> String {
        "CPU (one core, phase emulator)".into()
    }
    fn load_dampers(&mut self, _: &Program) -> Result<(), String> {
        Ok(())
    }
    fn run(
        &mut self,
        p: &Program,
        q: &mut [f64],
        v: &mut [f64],
        steps: usize,
        trace: &mut [f64],
        positions: &mut [f64],
    ) -> Result<(f64, f64), String> {
        Ok(self.0.run(p, q, v, steps, trace, positions))
    }
}

#[derive(Default)]
struct Stats {
    stretches: usize,
    steps: usize,
    rollbacks: usize,
    backend_s: f64,
    rebuilds: usize,
    rebuild_s: f64,
    hammer_s: f64,
    checked: usize,
    max_state: f64,
    max_trace: f64,
    max_position_m: f64,
    max_loss: f64,
}

pub struct Piano {
    backend: Box<dyn Backend>,
    required_backend: Option<String>,
    program: Program,
    damper_key: Option<(Vec<bool>, u64)>,
    positions: Vec<f64>,
    check: usize,
    stats: Stats,
    started: Instant,
}

/// Explicit GPU selection is a requirement and must not fall back to the CPU.
/// The factory's Option API cannot return initialization errors to its caller.
fn unavailable(choice: &str, reason: &str) -> Option<Box<dyn Accelerator>> {
    if matches!(choice, "cuda" | "metal") {
        eprintln!("render-gpu-piano: required {choice} backend unavailable: {reason}");
        std::process::exit(2);
    }
    eprintln!("render-gpu-piano: {reason}; rendering on the CPU");
    None
}

/// The accelerator for `instrument`, or None (with the reason on stderr) to render on the CPU. The
/// backend is chosen by FS_PIANO_GPU: auto (default: the GPU, else `cores`), metal, cuda, cores (the
/// phase emulator on the CPU cores), cpu (the same on one core) or off. FS_PIANO_GPU_CHECK=N compares
/// the first N steps of every stretch with the CPU's own step and limits stretches to N.
pub fn factory(instrument: &mut Instrument) -> Option<Box<dyn Accelerator>> {
    let choice = std::env::var("FS_PIANO_GPU").unwrap_or_else(|_| "auto".into());
    if choice == "off" {
        return None;
    }
    let Some((bank, dampers)) = instrument.gpu_images() else {
        return unavailable(&choice,
            "this configuration uses physics the GPU path does not model (felt or point dampers, radiation, cavity, string stretching or hammer footprints)");
    };
    let mut program = match Program::new(&bank, dampers.as_ref()) {
        Ok(p) => p,
        Err(e) => return unavailable(&choice, &e),
    };
    let (lifted, sustain) = instrument.damper_controls();
    if let Err(e) = program.set_dampers(dampers.as_ref(), &lifted, sustain) {
        return unavailable(&choice, &e);
    }
    let mut backend = match select(&choice, &program) {
        Ok(b) => b,
        Err(e) => {
            if choice == "auto" {
                eprintln!("render-gpu-piano: {e}; the phase emulator steps on the CPU cores");
                Box::new(crate::parallel::Cores(crate::parallel::Parallel::new(
                    &program,
                    crate::parallel::threads(),
                )))
            } else {
                eprintln!("render-gpu-piano: {e}");
                std::process::exit(2)
            }
        }
    };
    let check = std::env::var("FS_PIANO_GPU_CHECK")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let pads = program.np();
    let stations: usize = program.pads.iter().map(|p| p.stations()).sum();
    println!(
        "render-gpu-piano: contact-free stretches on {}; {} string modes in {} bridge groups, {} board coordinates, {} damper pads ({} stations), {} contacts{}.",
        backend.name(),
        program.n,
        program.ng,
        program.r,
        pads,
        stations,
        program.contacts.len(),
        if check > 0 {
            format!("; checking {check} steps per stretch against the CPU step")
        } else {
            String::new()
        }
    );
    if std::env::var("FS_PIANO_GPU_STATS").is_ok() {
        let mut groups: Vec<(usize, usize, usize, usize)> = (0..program.ng)
            .map(|g| {
                let (s, e) = program.group_modes[g];
                let (pb, pe) = program.group_pads[g];
                let work: usize = program.pads[pb..pe]
                    .iter()
                    .map(|p| p.stations() * p.len())
                    .sum();
                let st: usize = program.pads[pb..pe].iter().map(|p| p.stations()).sum();
                (e - s, pe - pb, st, work)
            })
            .collect();
        groups.sort_by_key(|x| std::cmp::Reverse(x.3));
        for (i, g) in groups.iter().take(12).enumerate() {
            eprintln!(
                "group #{i}: {} modes, {} pads, {} stations, {} station-modes",
                g.0, g.1, g.2, g.3
            );
        }
        let maxst = program.pads.iter().map(|p| p.stations()).max().unwrap_or(0);
        let maxm = program.pads.iter().map(|p| p.len()).max().unwrap_or(0);
        let total: usize = program.pads.iter().map(|p| p.stations() * p.len()).sum();
        eprintln!("pads: max {maxst} stations, max {maxm} modes, {total} station-modes in all");
    }
    if std::env::var("FS_PIANO_GPU_BENCH").is_ok() {
        bench(&mut *backend, &mut program);
    }
    let mut piano = Piano {
        backend,
        required_backend: matches!(choice.as_str(), "cuda" | "metal").then_some(choice.clone()),
        program,
        damper_key: Some((lifted, sustain.to_bits())),
        positions: Vec::new(),
        check,
        stats: Stats::default(),
        started: Instant::now(),
    };
    if let Err(e) = piano.backend.load_dampers(&piano.program) {
        return unavailable(&choice, &e);
    }
    Some(Box::new(piano))
}

/// FS_PIANO_GPU_BENCH: time stretches of several lengths from rest, with and without dampers.
fn bench(backend: &mut dyn Backend, program: &mut Program) {
    let (n, r, nc) = (program.n, program.r, program.contacts.len());
    let mut trace = vec![0.0; SEGMENT * r];
    let mut positions = vec![0.0; SEGMENT * nc.max(1)];
    let engaged = program.dampers.clone();
    // FS_PIANO_GPU_BENCH=a,b,...: the P_SKIP masks to time (kernel.h; 0 is the real kernel)
    let skips: Vec<i32> = std::env::var("FS_PIANO_GPU_BENCH")
        .ok()
        .map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect())
        .filter(|v: &Vec<i32>| !v.is_empty())
        .unwrap_or(vec![0, 1]);
    for (damp, skip) in skips.iter().flat_map(|&s| [(true, s), (false, s)]) {
        // SAFETY: the bench runs single-threaded before any rendering and exits afterwards
        unsafe {
            std::env::set_var("FS_PIANO_GPU_SKIP", skip.to_string());
        }
        program.dampers = engaged.clone();
        if !damp {
            program.dampers.any = false;
            program.dampers.active.fill(false);
        }
        if let Err(e) = backend.load_dampers(program) {
            eprintln!("bench: {e}");
            return;
        }
        for steps in [64usize, 2048] {
            let mut best = f64::INFINITY;
            let mut losses = (0.0, 0.0);
            for _ in 0..3 {
                let (mut q, mut v) = (vec![0.0; n + r], vec![1e-6; n + r]);
                let t = Instant::now();
                match backend.run(program, &mut q, &mut v, steps, &mut trace, &mut positions) {
                    Ok(l) => losses = l,
                    Err(e) => {
                        eprintln!("bench run failed: {e}");
                        return;
                    }
                }
                best = best.min(t.elapsed().as_secs_f64());
            }
            eprintln!(
                "bench {}: dampers {damp} skip {skip}: {steps} steps {:.3} ms ({:.1} us/step); losses {:.3e} {:.3e}",
                backend.name(),
                best * 1e3,
                best / steps as f64 * 1e6,
                losses.0,
                losses.1
            );
        }
    }
    std::process::exit(0);
}

fn select(choice: &str, program: &Program) -> Result<Box<dyn Backend>, String> {
    match choice {
        "cpu" => Ok(Box::new(Cpu(Emulator::new(program)))),
        "cores" => Ok(Box::new(crate::parallel::Cores(
            crate::parallel::Parallel::new(program, crate::parallel::threads()),
        ))),
        #[cfg(all(feature = "gpu", target_os = "macos"))]
        "auto" | "metal" => Ok(Box::new(crate::metal::Metal::new(program)?)),
        #[cfg(all(feature = "gpu", not(target_os = "macos")))]
        "auto" | "cuda" => Ok(Box::new(crate::cuda::Cuda::new(program)?)),
        #[cfg(not(feature = "gpu"))]
        "auto" => Err("no GPU backend in this build".into()),
        other => Err(format!(
            "FS_PIANO_GPU={other} is not available in this build (auto, metal on macOS, cuda elsewhere, cores, cpu, off)"
        )),
    }
}

impl Piano {
    /// One stretch on the backend. A backend failure changes nothing in the instrument, so the render
    /// may continue on the CPU phase emulator in automatic mode. Explicit GPU selection
    /// propagates the failure; `advance` restores the instrument's complete frame state.
    fn run(
        &mut self,
        q: &mut [f64],
        v: &mut [f64],
        steps: usize,
        trace: &mut [f64],
    ) -> Result<(f64, f64), String> {
        let (q0, v0) = (q.to_vec(), v.to_vec());
        match self
            .backend
            .run(&self.program, q, v, steps, trace, &mut self.positions)
        {
            Ok(losses) => Ok(losses),
            Err(e) => {
                if let Some(required) = &self.required_backend {
                    q.copy_from_slice(&q0);
                    v.copy_from_slice(&v0);
                    return Err(format!("required {required} backend failed: {e}"));
                }
                eprintln!(
                    "render-gpu-piano: {}: {e}; continuing on the CPU phase emulator",
                    self.backend.name()
                );
                let mut cpu = crate::parallel::Cores(crate::parallel::Parallel::new(
                    &self.program,
                    crate::parallel::threads(),
                ));
                cpu.load_dampers(&self.program)?;
                self.backend = Box::new(cpu);
                q.copy_from_slice(&q0);
                v.copy_from_slice(&v0);
                self.backend
                    .run(&self.program, q, v, steps, trace, &mut self.positions)
            }
        }
    }

    fn sync_dampers(&mut self, instrument: &mut Instrument) -> Result<(), String> {
        let (lifted, sustain) = instrument.damper_controls();
        let key = (lifted, sustain.to_bits());
        if self.damper_key.as_ref() == Some(&key) {
            return Ok(());
        }
        let t = Instant::now();
        let image = instrument.damper_image();
        self.program.set_dampers(image.as_ref(), &key.0, sustain)?;
        self.backend.load_dampers(&self.program)?;
        self.damper_key = Some(key);
        self.stats.rebuilds += 1;
        self.stats.rebuild_s += t.elapsed().as_secs_f64();
        Ok(())
    }

    /// Compare a backend stretch with the CPU's own steps from the same state; restores the bank.
    fn check(
        &mut self,
        instrument: &mut Instrument,
        initial: (&[f64], &[f64]),
        actual: (&[f64], &[f64]),
        steps: usize,
        trace: &[f64],
        losses: (f64, f64),
    ) -> Result<(), String> {
        let (q0, v0) = initial;
        let (q, v) = actual;
        let (n, r, nc) = (self.program.n, self.program.r, self.program.contacts.len());
        let (saved_q, saved_v) = (instrument.bank.q.clone(), instrument.bank.v.clone());
        instrument.bank.q.copy_from_slice(q0);
        instrument.bank.v.copy_from_slice(v0);
        let (mut modal, mut damper) = (0.0, 0.0);
        let scale = |x: &[f64]| x.iter().fold(0.0_f64, |m, y| m.max(y.abs())).max(1e-300);
        for s in 0..steps {
            let (m, d) = instrument
                .step_linear_reference()
                .map_err(|e| e.to_string())?;
            modal += m;
            damper += d;
            let reference = &instrument.bank.v[n..];
            let got = &trace[s * r..(s + 1) * r];
            let e = reference
                .iter()
                .zip(got)
                .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()))
                / scale(reference);
            self.stats.max_trace = self.stats.max_trace.max(e);
            for c in 0..nc {
                let x = instrument.bank.contact_position(c, &instrument.bank.q);
                self.stats.max_position_m = self
                    .stats
                    .max_position_m
                    .max((x - self.positions[s * nc + c]).abs());
            }
        }
        let rel = |a: &[f64], b: &[f64]| {
            a.iter()
                .zip(b)
                .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()))
                / scale(a)
        };
        let e = rel(&instrument.bank.q, q).max(rel(&instrument.bank.v, v));
        self.stats.max_state = self.stats.max_state.max(e);
        let le = ((modal - losses.0).abs() + (damper - losses.1).abs())
            / (modal.abs() + damper.abs()).max(1e-300);
        self.stats.max_loss = self.stats.max_loss.max(le);
        self.stats.checked += steps;
        if self.stats.checked == steps || e >= 1e-9 || e.is_nan() {
            println!(
                "render-gpu-piano check: {steps} steps, state rel {e:.2e}, board trace rel {:.2e}, contact positions {:.2e} m, losses rel {le:.2e}",
                self.stats.max_trace, self.stats.max_position_m
            );
        }
        instrument.bank.q.copy_from_slice(&saved_q);
        instrument.bank.v.copy_from_slice(&saved_v);
        Ok(())
    }
}

impl Accelerator for Piano {
    fn advance(
        &mut self,
        instrument: &mut Instrument,
        frames: usize,
        traces: &mut [f64],
    ) -> Result<usize, String> {
        let substeps = instrument.mechanics_substeps();
        let mut max = frames.min(SEGMENT / substeps) * substeps;
        if self.check > 0 {
            max = max.min(self.check);
        }
        max -= max % substeps;
        let planned = instrument.contact_free_steps(max, PLAN_CLEARANCE_M, PLAN_MARGIN);
        let steps = planned - planned % substeps;
        if steps < MIN_SEGMENT.min(max) || steps == 0 {
            return Ok(0);
        }
        let (r, nc) = (self.program.r, self.program.contacts.len());
        if traces.len() < steps * r {
            return Err("accelerator trace buffer cannot hold complete output frames".into());
        }
        self.sync_dampers(instrument)?;
        let (q0, v0) = (instrument.bank.q.clone(), instrument.bank.v.clone());
        let saved_hammers = instrument.checkpoint_external_hammers();
        let saved_modal = instrument.bank.last_modal_loss_j;
        let result = (|| {
            let (mut q, mut v) = (q0.clone(), v0.clone());
            let before = instrument.energy_j();
            self.positions.resize(steps * nc, 0.0);
            let t = Instant::now();
            let mut losses = self.run(&mut q, &mut v, steps, &mut traces[..steps * r])?;
            self.stats.backend_s += t.elapsed().as_secs_f64();
            if self.check > 0 {
                self.check(instrument, (&q0, &v0), (&q, &v), steps, traces, losses)?;
            }
            let t = Instant::now();
            let mut hammer = [0.0; 5];
            let mut taken = steps;
            for s in (0..steps).step_by(substeps) {
                match instrument
                    .advance_hammers_external_frame(
                        &self.positions[s * nc..(s + substeps) * nc],
                        EXACT_CLEARANCE_M,
                    )
                    .map_err(|e| e.to_string())?
                {
                    Some(e) => {
                        for (h, x) in hammer.iter_mut().zip(e) {
                            *h += x;
                        }
                    }
                    None => {
                        taken = s;
                        break;
                    }
                }
            }
            self.stats.hammer_s += t.elapsed().as_secs_f64();
            if taken < steps {
                self.stats.rollbacks += 1;
                if taken == 0 {
                    return Ok(0);
                }
                q.copy_from_slice(&q0);
                v.copy_from_slice(&v0);
                let t = Instant::now();
                losses = self.run(&mut q, &mut v, taken, &mut traces[..taken * r])?;
                self.stats.backend_s += t.elapsed().as_secs_f64();
            }
            instrument.bank.q.copy_from_slice(&q);
            instrument.bank.v.copy_from_slice(&v);
            instrument.absorb_external_losses(losses.0, losses.1);
            let after = instrument.energy_j();
            // jack work in; shank, felt (with its relaxation) and catch losses out
            let balance = after - before + losses.0 + losses.1 + hammer[1] + hammer[2] + hammer[4]
                - hammer[0];
            let tolerance = 1e-9 + 1e-10 * taken as f64 * before.abs().max(after.abs());
            if !balance.is_finite() || balance.abs() > tolerance {
                return Err(format!(
                    "accelerated stretch of {taken} steps does not close the energy ledger: defect {balance:.3e} J"
                ));
            }
            self.stats.stretches += 1;
            self.stats.steps += taken;
            Ok(taken / substeps)
        })();
        if !result.as_ref().is_ok_and(|done| *done > 0) {
            // No trace/frame is published on refusal, including a rerun or ledger error.
            instrument.bank.q.copy_from_slice(&q0);
            instrument.bank.v.copy_from_slice(&v0);
            instrument.bank.last_modal_loss_j = saved_modal;
            instrument.restore_external_hammers(&saved_hammers);
        }
        result
    }
}

impl Drop for Piano {
    fn drop(&mut self) {
        let s = &self.stats;
        let per = if s.steps > 0 {
            s.backend_s / s.steps as f64 * 1e6
        } else {
            0.0
        };
        println!(
            "render-gpu-piano: {} steps in {} stretches on {} ({per:.1} us/step, {:.1} s); {} rollbacks; {} reruns after a GPU stall; {} damper rebuilds ({:.2} s); hammer replay {:.2} s; {:.1} s since start.",
            s.steps,
            s.stretches,
            self.backend.name(),
            s.backend_s,
            s.rollbacks,
            self.backend.stalls(),
            s.rebuilds,
            s.rebuild_s,
            s.hammer_s,
            self.started.elapsed().as_secs_f64()
        );
        if s.checked > 0 {
            println!(
                "render-gpu-piano check: {} steps against the CPU step; worst state rel {:.2e}, board trace rel {:.2e}, contact positions {:.2e} m, losses rel {:.2e}.",
                s.checked, s.max_state, s.max_trace, s.max_position_m, s.max_loss
            );
        }
    }
}

#[cfg(test)]
#[path = "substeps_tests.rs"]
mod substeps_tests;
