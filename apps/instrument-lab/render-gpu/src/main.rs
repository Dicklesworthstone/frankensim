//! render-gpu: the fs-vkplate drum kit rendered on an NVIDIA GPU when one is present, on the CPU otherwise.
//!
//!   render-gpu PIECE HITS.csv OUT.wav SECONDS [--backend auto|cpu|cuda] [--airy f64|f32] [--device N|all]
//!   render-gpu --batch PIECE SECONDS OUT_PREFIX HITS1.csv [HITS2.csv ...]   (one simulation per list, in lockstep)
//!   render-gpu --check PIECE SECONDS [--airy f64|f32]                       (CPU vs GPU, same strokes)
//!   --exact-friction-order: frictions in the CPU's list order (slow; --check always uses it)
//!   render-gpu PIECE HITS.csv OUT.wav SECONDS --chunks N [--preroll S]       (N overlapping chunks at once)
//!
//! PIECE and HITS.csv are those of fs-vkplate's `perform` example. The CPU path IS `fs_vkplate::kit::render`;
//! the GPU path steps the same exported instrument (fs_vkplate::Instrument::image) in one persistent kernel.

#[cfg(all(feature = "cuda", not(target_os = "macos")))]
mod plan;
#[cfg(all(feature = "cuda", not(target_os = "macos")))]
mod cuda;

use fs_vkplate::kit::{self, Piece};
use std::time::Instant;

fn piece(name: &str, hits: &[(f64, f64)]) -> Result<Piece, String> {
    let cym = |r: f64, tip: f64| hits.iter().map(|&(t, v)| (t, v, r, tip)).collect::<Vec<_>>();
    match name {
        "kick" => kit::drum(&kit::kick(), hits),
        "snare" => kit::drum(&kit::snare(), hits),
        "snare-off" => kit::drum(&kit::snare_wires_off(), hits),
        "tom-12" => kit::drum(&kit::tom(12.0), hits),
        "tom-14" => kit::drum(&kit::tom(14.0), hits),
        "hihat-closed" => kit::hihat(0.0, 80.0, hits),
        "hihat-open" => kit::hihat(0.002, 0.0, hits),
        "crash" => kit::cymbal(&kit::crash(), &cym(0.95, 0.008)),
        "ride" => kit::cymbal(&kit::ride(), &cym(0.6, 0.004)),
        other => Err(format!("unknown piece {other}")),
    }
}

fn read_hits(path: &str) -> Vec<(f64, f64)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| die(&format!("{path}: {e}")));
    text.lines()
        .filter_map(|l| {
            let mut it = l.split(',').map(|s| s.trim().parse::<f64>());
            match (it.next(), it.next()) {
                (Some(Ok(t)), Some(Ok(v))) if t >= 0.0 && v > 0.0 => Some((t, v)),
                _ => None,
            }
        })
        .collect()
}

fn die(msg: &str) -> ! {
    eprintln!("render-gpu: {msg}");
    std::process::exit(1)
}

#[derive(Clone, Copy, PartialEq)]
enum Backend {
    Auto,
    Cpu,
    Cuda,
}

#[derive(Clone, Copy)]
struct Opts {
    backend: Backend,
    airy_f32: bool,
    device: usize,
    all_devices: bool,
    launch_s: f64,
    chunks: usize,
    preroll: f64,
    exact_friction: bool,
}

fn threads() -> usize {
    std::thread::available_parallelism().map_or(4, |n| n.get())
}

/// Whether to step on the GPU: `cuda` insists, `auto` falls back to the CPU (and says so).
fn use_gpu(o: &Opts) -> bool {
    match o.backend {
        Backend::Cpu => false,
        Backend::Cuda | Backend::Auto => {
            #[cfg(all(feature = "cuda", not(target_os = "macos")))]
            {
                match cuda::probe(o.device) {
                    Some(name) => {
                        eprintln!("GPU: {name} (CUDA device {})", o.device);
                        return true;
                    }
                    None if o.backend == Backend::Cuda => die("no CUDA driver or device"),
                    None => eprintln!("no CUDA device: rendering on the CPU"),
                }
            }
            #[cfg(not(all(feature = "cuda", not(target_os = "macos"))))]
            {
                if o.backend == Backend::Cuda {
                    die("no CUDA support in this build (built without the `cuda` feature, or for macOS)");
                }
                eprintln!("built without GPU support: rendering on the CPU");
            }
            false
        }
    }
}

/// Raw microphone samples of every simulation (internal rate, settle included) stepped on the GPU.
#[cfg(all(feature = "cuda", not(target_os = "macos")))]
fn gpu_raw(pieces: &mut [Piece], total: usize, o: &Opts) -> Vec<Vec<Vec<f64>>> {
    let t = Instant::now();
    let imgs: Vec<_> = pieces.iter_mut().map(|p| p.instrument.image()).collect();
    // nodes per work item (one block each): about one item per multiprocessor of a large GPU (128), at
    // most a full block
    let nodes: usize = imgs[0].bodies.iter().map(|b| b.n).sum();
    let item = (nodes * imgs.len()).div_ceil(128).next_multiple_of(32).clamp(64, 256);
    let pl = plan::Plan::new(&imgs[0], &pieces[0].instrument, threads(), item, o.exact_friction).unwrap_or_else(|e| die(&e));
    let sims: Vec<_> = imgs.iter().map(|i| pl.state(i).unwrap_or_else(|e| die(&e))).collect();
    let nvk = pl.vk_body.len();
    eprintln!(
        "plan: {} nodes, {} plate core nodes in {nvk} Airy inverse(s) ({:.0} MB as {}), {} contacts in {} batches, {} frictions in {} batches, {:.1} s",
        pl.nt,
        pl.nc,
        pl.airy.len() as f64 * if o.airy_f32 { 4.0 } else { 8.0 } / 1e6,
        if o.airy_f32 { "f32" } else { "f64" },
        pl.ct_gap.len(),
        pl.cb_ptr.len().saturating_sub(1),
        pl.fr_mu.len(),
        pl.fb_ptr.len().saturating_sub(1),
        t.elapsed().as_secs_f64()
    );
    // the batch split over the GPUs in use, one host thread per GPU
    let devices: Vec<usize> = if o.all_devices { (0..cuda::device_count()).collect() } else { vec![o.device] };
    let nm = pl.nm;
    let groups: Vec<std::ops::Range<usize>> = (0..devices.len())
        .map(|g| sims.len() * g / devices.len()..sims.len() * (g + 1) / devices.len())
        .filter(|r| !r.is_empty())
        .collect();
    let spl = ((o.launch_s * pl.rate) as usize).max(1);
    let first_dev = devices[0];
    let t = Instant::now();
    let results: Vec<(Vec<Vec<Vec<f64>>>, cuda::Run)> = std::thread::scope(|sc| {
        let handles: Vec<_> = groups
            .iter()
            .zip(&devices)
            .map(|(r, &dev)| {
                let (pl, sims) = (&pl, &sims[r.clone()]);
                sc.spawn(move || {
                    let b = sims.len();
                    let mut raw = vec![vec![Vec::with_capacity(total); nm]; b];
                    let mut last = Instant::now();
                    let run = cuda::run(pl, sims, total, o.airy_f32, dev, spl, |first, s| {
                        let steps = s.len() / (b * nm);
                        for i in 0..steps {
                            for sim in 0..b {
                                for m in 0..nm {
                                    raw[sim][m].push(s[(i * b + sim) * nm + m]);
                                }
                            }
                        }
                        if dev == first_dev && last.elapsed().as_secs_f64() > 2.0 {
                            last = Instant::now();
                            eprint!("\r{:5.1} %", 100.0 * (first + steps) as f64 / total as f64);
                        }
                    })
                    .unwrap_or_else(|e| die(&format!("device {dev}: {e}")));
                    (raw, run)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap_or_else(|_| die("a GPU thread panicked"))).collect()
    });
    let elapsed = t.elapsed().as_secs_f64();
    let audio = total as f64 / pl.rate;
    let mut raw = vec![];
    for (r, run) in results {
        eprintln!(
            "\r{} on {} blocks: {} x {:.2} s simulated in {:.2} s ({:.1} us/step, {:.2}x real time per simulation, {:.2}x for its batch; kernel {:.0} ms)",
            run.device_name,
            run.blocks,
            r.len(),
            audio,
            run.step_s,
            run.step_s / total as f64 * 1e6,
            audio / run.step_s,
            audio * r.len() as f64 / run.step_s,
            run.compile_ms
        );
        raw.extend(r);
    }
    eprintln!("{} simulation(s) on {} GPU(s) in {elapsed:.1} s ({:.2}x real time for the set)", raw.len(), groups.len(), audio * raw.len() as f64 / elapsed);
    raw
}

#[cfg(not(all(feature = "cuda", not(target_os = "macos"))))]
fn gpu_raw(_: &mut [Piece], _: usize, _: &Opts) -> Vec<Vec<Vec<f64>>> {
    unreachable!("use_gpu is false without the cuda feature")
}

fn cpu_raw(p: &mut Piece, total: usize) -> Vec<Vec<f64>> {
    let ins = &mut p.instrument;
    let mut raw = vec![Vec::with_capacity(total); ins.mics.len()];
    let mut buf = vec![0.0; ins.mics.len()];
    for _ in 0..total {
        ins.step(&mut buf);
        for (c, v) in raw.iter_mut().zip(&buf) {
            c.push(*v);
        }
    }
    raw
}

/// Raw microphone samples of every piece, stepped `total` times: one lockstep batch on the GPU, or one
/// thread per piece on the CPU (at most as many at once as there are cores).
fn raws(pieces: &mut [Piece], total: usize, o: &Opts) -> Vec<Vec<Vec<f64>>> {
    if use_gpu(o) {
        return gpu_raw(pieces, total, o);
    }
    let t = Instant::now();
    let mut out: Vec<Vec<Vec<f64>>> = vec![vec![]; pieces.len()];
    let lanes = threads().min(pieces.len()).max(1);
    for group in pieces.chunks_mut(lanes).zip(out.chunks_mut(lanes)) {
        std::thread::scope(|sc| {
            for (p, r) in group.0.iter_mut().zip(group.1.iter_mut()) {
                sc.spawn(move || *r = cpu_raw(p, total));
            }
        });
    }
    let audio = total as f64 / pieces[0].instrument.rate;
    let wall = t.elapsed().as_secs_f64();
    eprintln!(
        "CPU: {} x {audio:.2} s simulated in {wall:.1} s on {lanes} threads ({:.2}x real time for the set)",
        pieces.len(),
        audio * pieces.len() as f64 / wall
    );
    out
}

/// Half-width of the raised-cosine crossfade at each chunk boundary.
const XF: f64 = 0.002;

/// A long performance as `n` overlapping simulations stepped together (tools/strandbar's chunks.py,
/// validated against continuous renders within 0.06-0.3 dB per octave): each chunk starts `preroll`
/// seconds before its window with the real strokes of that pre-roll, so it enters the window already
/// ringing; cuts sit 10 ms before a stroke, and windows meet in a 4 ms raised-cosine crossfade. The joins
/// are not sample-exact: a cymbal pair is chaotic, so a restarted run decorrelates sample by sample while
/// keeping its spectrum and decay.
fn render_chunked(name: &str, hits: &[(f64, f64)], secs: f64, o: &Opts) -> Vec<Vec<f64>> {
    let n = o.chunks;
    let sr = kit::OUT_RATE;
    let mut cuts = vec![0.0];
    for c in 1..n {
        let (lo, hi) = (secs * c as f64 / n as f64, secs * (c + 1) as f64 / n as f64);
        let next = hits.iter().map(|h| h.0).filter(|&t| t >= lo).fold(f64::INFINITY, f64::min);
        let cut = if next < hi { next - 0.010 } else { lo };
        cuts.push(cut.max(cuts[c - 1] + 4.0 * XF));
    }
    cuts.push(secs);
    let mut pieces = vec![];
    let mut starts = vec![];
    let mut longest = 0.0f64;
    let mut lists = vec![];
    for c in 0..n {
        // whole output samples, so every stroke keeps its step relative to the output grid
        let s = ((cuts[c] - o.preroll).max(0.0) * sr).floor() / sr;
        let end = cuts[c + 1] + if c + 1 < n { XF } else { 0.0 };
        lists.push(hits.iter().filter(|&&(t, _)| t >= cuts[c] - o.preroll && t < cuts[c + 1] && t >= s).map(|&(t, v)| (t - s, v)).collect::<Vec<_>>());
        starts.push(s);
        longest = longest.max(end - s);
    }
    // every chunk gets the same number of strokes (a hi-hat adds a stick contact per stroke, and the batch
    // shares one contact list): the missing ones land after the end, where a not-yet-fired stroke is
    // exactly what a later real stroke is to the CPU model
    let most = lists.iter().map(|l| l.len()).max().unwrap_or(0);
    for l in lists.iter_mut() {
        l.resize(most, (longest + o.preroll + 10.0, 1.0));
    }
    // building an instrument factors its operators: the chunks' builds run side by side
    std::thread::scope(|sc| {
        let hs: Vec<_> = lists.iter().map(|l| sc.spawn(move || piece(name, l).unwrap_or_else(|e| die(&e)))).collect();
        pieces.extend(hs.into_iter().map(|h| h.join().unwrap_or_else(|_| die("a build panicked"))));
    });
    let rate = pieces[0].instrument.rate;
    let total = ((longest + pieces[0].settle) * rate).round() as usize;
    let skip = (pieces[0].settle * rate).round() as usize;
    eprintln!("{n} chunks of at most {longest:.2} s (pre-roll {:.2} s) for {secs} s", o.preroll);
    let raw = raws(&mut pieces, total, o);
    let nm = raw[0].len();
    let len = (secs * sr).round() as usize + 1;
    let mut y = vec![vec![0.0; len]; nm];
    let nx = (2.0 * XF * sr).round() as usize;
    let ramp: Vec<f64> = (0..nx).map(|i| 0.5 - 0.5 * (std::f64::consts::PI * i as f64 / nx as f64).cos()).collect();
    for c in 0..n {
        let a = cuts[c] - if c > 0 { XF } else { 0.0 };
        let b = cuts[c + 1] + if c + 1 < n { XF } else { 0.0 };
        let (i0, i1) = ((a * sr).round() as usize, ((b * sr).round() as usize).min(len));
        let off = ((a - starts[c]) * sr).round() as usize;
        for m in 0..nm {
            let x = fs_vkplate::decimate(&raw[c][m][skip.min(raw[c][m].len())..], pieces[c].oversample);
            for (j, yi) in y[m][i0..i1].iter_mut().enumerate() {
                let mut w = 1.0;
                if c > 0 && j < nx {
                    w = ramp[j];
                }
                if c + 1 < n && i1 - i0 - j <= nx {
                    w = ramp[i1 - i0 - j - 1];
                }
                *yi += x.get(off + j).copied().unwrap_or(0.0) * w;
            }
        }
    }
    y
}

fn finish(p: &Piece, raw: &[Vec<f64>], skip: usize, out: &str) {
    let ch: Vec<Vec<f64>> = raw.iter().map(|c| fs_vkplate::decimate(&c[skip.min(c.len())..], p.oversample)).collect();
    let peak = ch.iter().flatten().fold(0.0f64, |m, &x| m.max(x.abs()));
    if !peak.is_finite() {
        die("non-finite output");
    }
    std::fs::write(out, fs_vkplate::wav_f32(&ch, kit::OUT_RATE as u32)).unwrap_or_else(|e| die(&format!("{out}: {e}")));
    eprintln!("{out}: peak {peak:.3} Pa");
}

fn fft(re: &mut [f64], im: &mut [f64]) {
    let n = re.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let a = -2.0 * std::f64::consts::PI / len as f64;
        for s in (0..n).step_by(len) {
            for t in 0..len / 2 {
                let (c, sn) = ((a * t as f64).cos(), (a * t as f64).sin());
                let (u, v) = (s + t, s + t + len / 2);
                let (xr, xi) = (re[v] * c - im[v] * sn, re[v] * sn + im[v] * c);
                re[v] = re[u] - xr;
                im[v] = im[u] - xi;
                re[u] += xr;
                im[u] += xi;
            }
        }
        len <<= 1;
    }
}

/// Energy per octave band (63 Hz .. 16 kHz centres) of a signal at `rate`, in dB.
fn octaves(x: &[f64], rate: f64) -> Vec<(f64, f64)> {
    let n = x.len().next_power_of_two();
    let (mut re, mut im) = (vec![0.0; n], vec![0.0; n]);
    re[..x.len()].copy_from_slice(x);
    fft(&mut re, &mut im);
    let mut out = vec![];
    let mut fc = 62.5;
    while fc <= 16000.0 {
        let (lo, hi) = (fc / 2f64.sqrt(), fc * 2f64.sqrt());
        let (a, b) = ((lo / rate * n as f64) as usize, ((hi / rate * n as f64) as usize).min(n / 2));
        let e: f64 = (a..b).map(|i| re[i] * re[i] + im[i] * im[i]).sum();
        out.push((fc, 10.0 * e.max(1e-300).log10()));
        fc *= 2.0;
    }
    out
}

fn check(name: &str, secs: f64, o: &Opts) {
    // the CPU's own friction order, so any difference is the arithmetic's
    let o = &Opts { exact_friction: true, ..*o };
    let hits = [(0.002, 3.0), (secs * 0.5, 2.0)];
    let mut cpu = piece(name, &hits).unwrap_or_else(|e| die(&e));
    let mut gpu = vec![piece(name, &hits).unwrap_or_else(|e| die(&e))];
    let rate = cpu.instrument.rate;
    let total = ((secs + cpu.settle) * rate).round() as usize;
    if !use_gpu(o) {
        die("--check needs a GPU");
    }
    let tg = Instant::now();
    let g = gpu_raw(&mut gpu, total, o).remove(0);
    let wg = tg.elapsed().as_secs_f64();
    let tc = Instant::now();
    let c = cpu_raw(&mut cpu, total);
    let wc = tc.elapsed().as_secs_f64();
    println!("{name}: {total} steps at {rate:.0} Hz; CPU {wc:.2} s, GPU {wg:.2} s (incl. setup) -> {:.1}x", wc / wg);
    for m in 0..c.len() {
        let rms = (c[m].iter().map(|x| x * x).sum::<f64>() / c[m].len() as f64).sqrt();
        println!("mic {m}: rms {rms:.3e} Pa; relative difference |gpu - cpu| / rms over time:");
        let win = (0.002 * rate) as usize;
        let mut first_big = None;
        let mut line = String::new();
        for (wi, (a, b)) in c[m].chunks(win).zip(g[m].chunks(win)).enumerate() {
            let d = (a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum::<f64>() / a.len() as f64).sqrt() / rms.max(1e-300);
            if wi < 8 || wi % (c[m].len() / win / 12).max(1) == 0 {
                line += &format!(" {:.0}ms:{d:.1e}", wi as f64 * 2.0);
            }
            if first_big.is_none() && d > 1e-3 {
                first_big = Some(wi as f64 * 2.0);
            }
        }
        println!("  {line}");
        match first_big {
            Some(ms) => println!("  first 2-ms window above 1e-3: {ms:.0} ms"),
            None => println!("  never above 1e-3"),
        }
        let skip = (cpu.settle * rate) as usize;
        let (oc, og) = (octaves(&c[m][skip..], rate), octaves(&g[m][skip..], rate));
        let row: Vec<String> = oc.iter().zip(&og).map(|(a, b)| format!("{:.0}:{:+.2}", a.0, b.1 - a.1)).collect();
        println!("  octave bands, GPU - CPU dB: {}", row.join(" "));
    }
}

fn main() {
    let mut pos = vec![];
    let mut o = Opts { backend: Backend::Auto, airy_f32: false, device: 0, all_devices: false, launch_s: 0.25, chunks: 1, preroll: 1.0, exact_friction: false };
    let (mut batch, mut chk) = (false, false);
    let mut a = std::env::args().skip(1);
    while let Some(x) = a.next() {
        match x.as_str() {
            "--backend" => {
                o.backend = match a.next().as_deref() {
                    Some("auto") => Backend::Auto,
                    Some("cpu") => Backend::Cpu,
                    Some("cuda") => Backend::Cuda,
                    other => die(&format!("--backend {other:?}: auto, cpu or cuda")),
                }
            }
            "--airy" => {
                o.airy_f32 = match a.next().as_deref() {
                    Some("f64") => false,
                    Some("f32") => true,
                    other => die(&format!("--airy {other:?}: f64 or f32")),
                }
            }
            "--device" => match a.next().as_deref() {
                Some("all") => o.all_devices = true,
                v => o.device = v.and_then(|s| s.parse().ok()).unwrap_or_else(|| die("--device N|all")),
            },
            "--launch-seconds" => o.launch_s = a.next().and_then(|s| s.parse().ok()).unwrap_or_else(|| die("--launch-seconds S")),
            "--chunks" => o.chunks = a.next().and_then(|s| s.parse().ok()).filter(|&n| n >= 1).unwrap_or_else(|| die("--chunks N")),
            "--preroll" => o.preroll = a.next().and_then(|s| s.parse().ok()).unwrap_or_else(|| die("--preroll SECONDS")),
            #[cfg(all(feature = "cuda", not(target_os = "macos")))]
            "--precompile" => {
                let arch = a.next().unwrap_or_else(|| die("--precompile sm_XY"));
                cuda::precompile(&arch).unwrap_or_else(|e| die(&e));
                return;
            }
            "--exact-friction-order" => o.exact_friction = true,
            "--batch" => batch = true,
            "--check" => chk = true,
            _ => pos.push(x),
        }
    }
    let usage = "usage: render-gpu PIECE HITS.csv OUT.wav SECONDS | --batch PIECE SECONDS OUT_PREFIX HITS.csv... | --check PIECE SECONDS";
    if chk {
        if pos.len() != 2 {
            die(usage);
        }
        check(&pos[0], pos[1].parse().unwrap_or_else(|_| die("SECONDS")), &o);
        return;
    }
    if batch {
        if pos.len() < 4 {
            die(usage);
        }
        let secs: f64 = pos[1].parse().unwrap_or_else(|_| die("SECONDS"));
        let mut lists: Vec<Vec<(f64, f64)>> = pos[3..].iter().map(|h| read_hits(h)).collect();
        let most = lists.iter().map(|l| l.len()).max().unwrap_or(0);
        for l in lists.iter_mut() {
            l.resize(most, (secs + 10.0, 1.0)); // same stroke count for every member (see render_chunked)
        }
        let mut pieces: Vec<Piece> = lists.iter().map(|l| piece(&pos[0], l).unwrap_or_else(|e| die(&e))).collect();
        let rate = pieces[0].instrument.rate;
        let total = ((secs + pieces[0].settle) * rate).round() as usize;
        let skip = (pieces[0].settle * rate).round() as usize;
        let all = raws(&mut pieces, total, &o);
        for (i, raw) in all.iter().enumerate() {
            finish(&pieces[i], raw, skip, &format!("{}{i:03}.wav", pos[2]));
        }
        return;
    }
    if pos.len() != 4 {
        die(usage);
    }
    let hits = read_hits(&pos[1]);
    let secs: f64 = pos[3].parse().unwrap_or_else(|_| die("SECONDS"));
    let mut p = piece(&pos[0], &hits).unwrap_or_else(|e| die(&e));
    for n in &p.notes {
        eprintln!("{n}");
    }
    let (ratio, who) = p.instrument.stability();
    if ratio >= 0.95 {
        die(&format!("leapfrog would be unstable ({ratio:.3}, {who})"));
    }
    let rate = p.instrument.rate;
    let total = ((secs + p.settle) * rate).round() as usize;
    let skip = (p.settle * rate).round() as usize;
    eprintln!("{} strokes over {secs} s", hits.len());
    if o.chunks > 1 {
        let t = Instant::now();
        let y = render_chunked(&pos[0], &hits, secs, &o);
        let peak = y.iter().flatten().fold(0.0f64, |m, &x| m.max(x.abs()));
        if !peak.is_finite() {
            die("non-finite output");
        }
        std::fs::write(&pos[2], fs_vkplate::wav_f32(&y, kit::OUT_RATE as u32)).unwrap_or_else(|e| die(&e.to_string()));
        eprintln!("{}: peak {peak:.3} Pa; {secs} s in {:.1} s overall", pos[2], t.elapsed().as_secs_f64());
        return;
    }
    if use_gpu(&o) {
        let raw = gpu_raw(std::slice::from_mut(&mut p), total, &o).remove(0);
        finish(&p, &raw, skip, &pos[2]);
    } else {
        let t = Instant::now();
        let out = kit::render(&mut p, secs, |f| eprint!("\r{:5.1} %", 100.0 * f));
        let wall = t.elapsed().as_secs_f64();
        eprintln!("\r{secs} s rendered on the CPU in {wall:.1} s ({:.1}x real time)", secs / wall);
        std::fs::write(&pos[2], fs_vkplate::wav_f32(&out, kit::OUT_RATE as u32)).unwrap_or_else(|e| die(&e.to_string()));
    }
}
