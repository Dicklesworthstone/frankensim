//! render-gpu: the fs-vkplate drum kit rendered on an NVIDIA GPU when one is present, on the CPU otherwise.
//!
//!   render-gpu PIECE HITS.csv OUT.wav SECONDS [--backend auto|cpu|cuda] [--airy f64|f32] [--device N|all]
//!   render-gpu --batch PIECE SECONDS OUT_PREFIX HITS1.csv [HITS2.csv ...]   (one simulation per list, in lockstep)
//!   render-gpu --check PIECE SECONDS [--airy f64|f32]                       (CPU vs GPU, same strokes)
//!   --exact-friction-order: frictions in the CPU's list order (slow; --check always uses it)
//!   --preset-heads: render a drum's recording-fitted preset heads; without it, drum pieces refuse (.370)
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

/// The default refuses drums: their only card-free heads are the presets, whose T60 was fitted to
/// recordings (bead .370, DONE-WHEN 1). The sourced heads carry a head-radiation load this GPU path does
/// not have, so render them natively (fs-vkplate's examples). `--preset-heads` renders the fitted heads
/// on purpose, and says so in the notes.
pub const PRESET_HEADS_REFUSAL: &str = "drum heads: the default refuses the card-free preset heads, whose T60 was fitted to recordings (.370); render the sourced heads natively (fs-vkplate examples), or pass --preset-heads to render the fitted heads on purpose";
pub const PRESET_HEADS_RECEIPT: &str = "heads: PRESET, recording-fitted T60 (opt-in --preset-heads); not a sourced head loss";

fn piece(name: &str, hits: &[(f64, f64)], preset_heads: bool, grid_scale: f64) -> Result<Piece, String> {
    // Bound this diagnostic knob before allocating/factoring a mesh. In
    // particular, tiny finite scales can overflow the grid's integer extent.
    if !grid_scale.is_finite() || !(0.5..=8.0).contains(&grid_scale) {
        return Err("--grid-scale must be within 0.5..=8 (mesh resource bound)".into());
    }
    let cym = |r: f64, tip: f64| hits.iter().map(|&(t, v)| (t, v, r, tip)).collect::<Vec<_>>();
    let grid = |mut spec: kit::CymbalSpec| { spec.grid *= grid_scale; spec };
    let drum = |spec: kit::DrumSpec| {
        if grid_scale != 1.0 {
            return Err("--grid-scale currently supports hi-hat, crash and ride only".into());
        }
        if !preset_heads {
            return Err(PRESET_HEADS_REFUSAL.to_string());
        }
        let mut p = kit::drum(&spec, hits)?;
        p.notes.push(PRESET_HEADS_RECEIPT.into());
        Ok(p)
    };
    match name {
        "kick" => drum(kit::kick()),
        "snare" => drum(kit::snare()),
        "snare-off" => drum(kit::snare_wires_off()),
        "tom-12" => drum(kit::tom(12.0)),
        "tom-14" => drum(kit::tom(14.0)),
        "hihat-closed" => kit::hihat_scaled(0.0, 80.0, hits, grid_scale),
        "hihat-open" => kit::hihat_scaled(0.002, 0.0, hits, grid_scale),
        "crash" => kit::cymbal(&grid(kit::crash()), &cym(0.95, 0.008)),
        "ride" => kit::cymbal(&grid(kit::ride()), &cym(0.6, 0.004)),
        other => Err(format!("unknown piece {other}")),
    }
}

/// `time_s,speed_m_s` lines (a header and malformed lines are skipped), in time order: the GPU fires strokes
/// from a cursor, and a stroke list out of order would fire differently there and on the CPU.
fn read_hits(path: &str) -> Vec<(f64, f64)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| die(&format!("{path}: {e}")));
    let mut hits: Vec<(f64, f64)> = text
        .lines()
        .filter_map(|l| {
            let mut it = l.split(',').map(|s| s.trim().parse::<f64>());
            match (it.next(), it.next()) {
                (Some(Ok(t)), Some(Ok(v))) if t >= 0.0 && v > 0.0 => Some((t, v)),
                _ => None,
            }
        })
        .collect();
    hits.sort_by(|a, b| a.0.total_cmp(&b.0));
    hits
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
    preset_heads: bool,
    grid_scale: f64,
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
    let gpus = if o.all_devices { cuda::device_count().max(1) } else { 1 };
    let item = (nodes * imgs.len().div_ceil(gpus)).div_ceil(128).next_multiple_of(32).clamp(64, 256);
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
            "\r{} on {} blocks: {} x {:.2} s simulated in {:.2} s ({:.1} us/step, {:.2}x real time per simulation, {:.2}x for its batch; kernel compilation {:.0} ms)",
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
        let hs: Vec<_> = lists.iter().map(|l| sc.spawn(move || piece(name, l, o.preset_heads, o.grid_scale).unwrap_or_else(|e| die(&e)))).collect();
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
    // joins on whole output samples: chunk c owns [join_c - h, join_{c+1} + h), fading in over the first 2h
    // samples and out over the last 2h; the ramps are sampled at half-sample offsets so the two that meet
    // at a join sum to exactly one
    let h = (XF * sr).round() as usize;
    let join: Vec<usize> = cuts.iter().map(|&c| ((c * sr).round() as usize).min(len)).collect();
    let ramp: Vec<f64> = (0..2 * h).map(|i| 0.5 - 0.5 * (std::f64::consts::PI * (i as f64 + 0.5) / (2 * h) as f64).cos()).collect();
    for c in 0..n {
        let i0 = if c > 0 { join[c].saturating_sub(h) } else { 0 };
        let i1 = if c + 1 < n { (join[c + 1] + h).min(len) } else { len };
        let start = (starts[c] * sr).round() as usize;
        for m in 0..nm {
            let x = fs_vkplate::decimate(&raw[c][m][skip.min(raw[c][m].len())..], pieces[c].oversample);
            for (g, yi) in y[m].iter_mut().enumerate().take(i1).skip(i0) {
                let mut w = 1.0;
                if c > 0 && g < join[c] + h {
                    w *= ramp[g + h - join[c]];
                }
                if c + 1 < n && g + h >= join[c + 1] {
                    w *= ramp[2 * h - 1 - (g + h - join[c + 1])];
                }
                *yi += g.checked_sub(start).and_then(|i| x.get(i)).copied().unwrap_or(0.0) * w;
            }
        }
    }
    y
}

fn encode_audio(channels: &[Vec<f64>]) -> Result<(Vec<u8>, f64), String> {
    let mut peak = 0.0_f64;
    for (channel, samples) in channels.iter().enumerate() {
        for (frame, &sample) in samples.iter().enumerate() {
            // f64::max ignores NaN. Validate each sample, including its WAV
            // conversion, before encoding or opening the destination file.
            if !(sample as f32).is_finite() {
                return Err(format!("invalid float32 output at channel {channel}, frame {frame}"));
            }
            peak = peak.max(sample.abs());
        }
    }
    Ok((fs_vkplate::wav_f32(channels, kit::OUT_RATE as u32), peak))
}

fn write_audio(channels: &[Vec<f64>], out: &str) -> f64 {
    let (bytes, peak) = encode_audio(channels).unwrap_or_else(|e| die(&e));
    std::fs::write(out, bytes).unwrap_or_else(|e| die(&format!("{out}: {e}")));
    peak
}

fn finish(p: &Piece, raw: &[Vec<f64>], skip: usize, out: &str) {
    let ch: Vec<Vec<f64>> = raw.iter().map(|c| fs_vkplate::decimate(&c[skip.min(c.len())..], p.oversample)).collect();
    let peak = write_audio(&ch, out);
    eprintln!("{out}: peak {peak:.3} Pa");
}

#[cfg(test)]
mod audio_output_tests {
    use super::encode_audio;

    // G0: reject invalid output in any channel before WAV conversion can
    // hide a NaN in the peak reduction or overflow a finite f64 sample.
    #[test]
    fn refuses_nonfinite_and_float32_overflow_in_later_channels() {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY,
                        2.0 * f64::from(f32::MAX), -2.0 * f64::from(f32::MAX)] {
            let channels = [vec![0.0, 0.5, 0.0], vec![0.25, invalid, 0.0]];
            assert_eq!(encode_audio(&channels).unwrap_err(),
                "invalid float32 output at channel 1, frame 1");
        }
        assert!(encode_audio(&[vec![f64::NAN; 3]]).is_err());
    }

    #[test]
    fn preserves_valid_samples_and_interleaving_without_clipping_pressure() {
        let channels = [vec![-0.0, 2.5, -0.25], vec![0.125, -3.0, 0.5]];
        let (bytes, peak) = encode_audio(&channels).unwrap();
        assert_eq!(peak, 3.0);
        assert_eq!(&bytes[..4], b"RIFF");
        assert_eq!(&bytes[8..16], b"WAVEfmt ");
        assert_eq!(u16::from_le_bytes(bytes[20..22].try_into().unwrap()), 3);
        assert_eq!(u16::from_le_bytes(bytes[22..24].try_into().unwrap()), 2);
        assert_eq!(u32::from_le_bytes(bytes[24..28].try_into().unwrap()), 48000);
        assert_eq!(&bytes[36..40], b"data");
        assert_eq!(u32::from_le_bytes(bytes[40..44].try_into().unwrap()), 24);
        let expected = [-0.0_f32, 0.125, 2.5, -3.0, -0.25, 0.5];
        assert_eq!(bytes.len(), 44 + 4 * expected.len());
        for (encoded, sample) in bytes[44..].chunks_exact(4).zip(expected) {
            assert_eq!(u32::from_le_bytes(encoded.try_into().unwrap()), sample.to_bits());
        }
    }
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

/// Validate the complete trajectories before accepting the existing 2-ms,
/// 1e-3 relative-RMS comparison. Zip alone would hide a missing GPU tail,
/// and a comparison with NaN never satisfies `d > threshold`.
fn compare_raw(c: &[Vec<f64>], g: &[Vec<f64>], total: usize, rate: f64, skip: usize) -> Result<(), String> {
    if !rate.is_finite() || rate <= 0.0 || total == 0 || skip >= total {
        return Err("--check needs a finite positive rate and a nonempty post-settle trajectory".into());
    }
    if c.is_empty() || c.len() != g.len() {
        return Err("--check CPU/GPU microphone counts differ or are empty".into());
    }
    for (m, (cpu, gpu)) in c.iter().zip(g).enumerate() {
        if cpu.len() != total || gpu.len() != total {
            return Err(format!("--check mic {m}: expected {total} samples, got CPU {} and GPU {}", cpu.len(), gpu.len()));
        }
        if cpu.iter().chain(gpu).any(|x| !x.is_finite()) {
            return Err(format!("--check mic {m}: nonfinite trajectory"));
        }
    }
    let win = ((0.002 * rate) as usize).max(1);
    let mut failed = false;
    let mut audible = false;
    for m in 0..c.len() {
        let rms = (c[m].iter().map(|x| x * x).sum::<f64>() / c[m].len() as f64).sqrt();
        if !rms.is_finite() {
            return Err(format!("--check mic {m}: nonfinite reference energy"));
        }
        audible |= rms > 0.0;
        println!("mic {m}: rms {rms:.3e} Pa; relative difference |gpu - cpu| / rms over time:");
        let mut first_big = None;
        let mut line = String::new();
        for (wi, (a, b)) in c[m].chunks(win).zip(g[m].chunks(win)).enumerate() {
            let d = (a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum::<f64>() / a.len() as f64).sqrt() / rms.max(1e-300);
            if !d.is_finite() {
                return Err(format!("--check mic {m}: nonfinite difference in window {wi}"));
            }
            if wi < 8 || wi % (c[m].len() / win / 12).max(1) == 0 {
                line += &format!(" {:.0}ms:{d:.1e}", wi as f64 * 2.0);
            }
            if first_big.is_none() && d > 1e-3 {
                first_big = Some(wi as f64 * 2.0);
            }
        }
        println!("  {line}");
        match first_big {
            Some(ms) => { println!("  first 2-ms window above 1e-3: {ms:.0} ms"); failed = true; }
            None => println!("  never above 1e-3"),
        }
        let (oc, og) = (octaves(&c[m][skip..], rate), octaves(&g[m][skip..], rate));
        let row: Vec<String> = oc.iter().zip(&og).map(|(a, b)| format!("{:.0}:{:+.2}", a.0, b.1 - a.1)).collect();
        println!("  octave bands, GPU - CPU dB: {}", row.join(" "));
    }
    if failed { return Err("CPU/CUDA comparison exceeded the 1e-3 relative-RMS threshold".into()); }
    if !audible { return Err("--check reference is entirely silent; parity is unqualified".into()); }
    Ok(())
}

/// Independent trajectories from one constructed initial state. Constructing the
/// same mesh twice repeats factorization and the shell stability estimate;
/// cloning preserves their results without sharing mutable simulation state.
fn comparison_pieces(name: &str, hits: &[(f64, f64)], preset_heads: bool, grid_scale: f64) -> Result<(Piece, Piece), String> {
    let cpu = piece(name, hits, preset_heads, grid_scale)?;
    let gpu = Piece {
        instrument: cpu.instrument.clone(),
        oversample: cpu.oversample,
        settle: cpu.settle,
        notes: cpu.notes.clone(),
    };
    Ok((cpu, gpu))
}

fn check(name: &str, secs: f64, o: &Opts) {
    if !secs.is_finite() || secs <= 0.0 {
        die("--check SECONDS must be finite and positive");
    }
    // the CPU's own friction order, so any difference is the arithmetic's
    let o = &Opts { exact_friction: true, ..*o };
    let hits = [(0.002, 3.0), (secs * 0.5, 2.0)];
    let started = Instant::now();
    let (mut cpu, gpu) = comparison_pieces(name, &hits, o.preset_heads, o.grid_scale).unwrap_or_else(|e| die(&e));
    eprintln!("comparison initial state: {:.2} s (one construction, two independent trajectories)", started.elapsed().as_secs_f64());
    let mut gpu = vec![gpu];
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
    compare_raw(&c, &g, total, rate, (cpu.settle * rate) as usize).unwrap_or_else(|e| die(&e));
}

#[cfg(test)]
mod trajectory_check_tests {
    use super::{compare_raw, comparison_pieces, cpu_raw, piece};

    #[test]
    fn copied_initial_state_matches_independent_build_without_sharing_trajectory() {
        for name in ["hihat-open", "crash"] {
            let hits = [(0.002, 3.0), (0.006, 2.0)];
            let (mut original, mut copied) = comparison_pieces(name, &hits, false, 6.0).unwrap();
            let mut fresh = piece(name, &hits, false, 6.0).unwrap();
            assert_eq!(original.oversample, copied.oversample);
            assert_eq!(original.settle.to_bits(), copied.settle.to_bits());
            assert_eq!(original.notes, copied.notes);
            assert_eq!(original.instrument.rate.to_bits(), copied.instrument.rate.to_bits());
            let total = (0.012 * original.instrument.rate).round() as usize;
            let expected = cpu_raw(&mut original, total);
            assert_eq!(copied.instrument.step_index, 0);
            assert!(expected.iter().flatten().any(|&v| v != 0.0));
            for actual in [cpu_raw(&mut fresh, total), cpu_raw(&mut copied, total)] {
                assert_eq!(actual.len(), expected.len());
                for (a, e) in actual.iter().zip(&expected) {
                    assert_eq!(a.len(), e.len());
                    for (&a, &e) in a.iter().zip(e) {
                        assert!(a.is_finite() && e.is_finite());
                        assert_eq!(a.to_bits(), e.to_bits(), "{name}");
                    }
                }
            }
            assert_eq!(copied.instrument.step_index, original.instrument.step_index);
        }
    }

    #[test]
    fn full_matching_trajectory_passes_and_partial_last_window_is_checked() {
        let cpu = vec![vec![1.0; 193]]; // 96 samples/window, including a final singleton.
        assert!(compare_raw(&cpu, &cpu, 193, 48000.0, 0).is_ok());
        let mut gpu = cpu.clone();
        gpu[0][192] += 0.002;
        assert!(compare_raw(&cpu, &gpu, 193, 48000.0, 0).unwrap_err().contains("threshold"));
        gpu[0][192] = 1.0005;
        assert!(compare_raw(&cpu, &gpu, 193, 48000.0, 0).is_ok());
    }

    #[test]
    fn missing_samples_mics_nonfinite_and_silent_evidence_refuse() {
        let cpu = vec![vec![1.0; 193]];
        let short = vec![vec![1.0; 192]];
        assert!(compare_raw(&cpu, &short, 193, 48000.0, 0).is_err());
        assert!(compare_raw(&short, &cpu, 193, 48000.0, 0).is_err());
        assert!(compare_raw(&cpu, &[], 193, 48000.0, 0).is_err());
        assert!(compare_raw(&[], &[], 193, 48000.0, 0).is_err());
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut bad = cpu.clone(); bad[0][192] = invalid;
            assert!(compare_raw(&cpu, &bad, 193, 48000.0, 0).is_err());
            assert!(compare_raw(&bad, &cpu, 193, 48000.0, 0).is_err());
        }
        let silent = vec![vec![0.0; 193]];
        assert!(compare_raw(&silent, &silent, 193, 48000.0, 0).is_err());
        assert!(compare_raw(&cpu, &cpu, 193, f64::NAN, 0).is_err());
        assert!(compare_raw(&cpu, &cpu, 193, 48000.0, 193).is_err());
    }
}

fn main() {
    let mut pos = vec![];
    let mut o = Opts { backend: Backend::Auto, airy_f32: false, device: 0, all_devices: false, launch_s: 0.25, chunks: 1, preroll: 1.0, exact_friction: false, preset_heads: false, grid_scale: 1.0 };
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
            "--preset-heads" => o.preset_heads = true,
            "--grid-scale" => o.grid_scale = a.next().and_then(|s| s.parse::<f64>().ok())
                .filter(|s| s.is_finite() && (0.5..=8.0).contains(s))
                .unwrap_or_else(|| die("--grid-scale must be within 0.5..=8 (mesh resource bound)")),
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
        let mut pieces: Vec<Piece> = lists.iter().map(|l| piece(&pos[0], l, o.preset_heads, o.grid_scale).unwrap_or_else(|e| die(&e))).collect();
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
    let mut p = piece(&pos[0], &hits, o.preset_heads, o.grid_scale).unwrap_or_else(|e| die(&e));
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
        let peak = write_audio(&y, &pos[2]);
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
        write_audio(&out, &pos[2]);
    }
}
