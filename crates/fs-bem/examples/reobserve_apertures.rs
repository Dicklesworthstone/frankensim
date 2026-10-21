//! Re-observe solved aperture flows through a closed exterior BEM body.
//!
//! `bake BODY RECEIVERS OUT [--threads N]` meshes a lathed body with named
//! aperture patches (`fs_bem::lathed_apertures`), binds every aperture's unit
//! volume flow and every receiver to `BemObservation`, solves the frequency
//! samples in parallel, and runs the stock `bake_observation_transfers` over
//! them. It writes OUT.kernels (pressure kernels only) and OUT.json.
//!
//! `replay KERNELS FLOWS BIND OUTDIR` reads an aperture-flow file written by
//! `music_render wind --aperture-flows`, decimates its mechanical-rate flows to
//! the kernel clock with a declared Blackman-windowed sinc, binds each tapped
//! node to a body aperture through BIND, and writes one 32-bit float WAV of
//! pressure [Pa] per receiver plus OUTDIR/replay.json.
//!
//! No source is changed: the flows are replayed exactly as solved. The body,
//! receivers, band and DFT ladder are explicit inputs; nothing is defaulted.
use fs_bem::helmholtz::Medium;
use fs_bem::lathed_apertures::{LathedApertureSpec, WallAperture, lathed_aperture_body};
use fs_bem::observation_transfers::{
    BemObservation, ProducerUncertainty, ReceiverSpec, SourceId, TransferDeclaration,
    bake_observation_transfers,
};
use fs_math::c64::C64;
use fs_vfit::observation::{FirKernel, MAX_OBSERVATION_ENTRIES, ObservationFir};
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::time::Instant;

type Res<T> = Result<T, String>;
const KERNELS_SCHEMA: &str = "frankensim-observation-kernels-v1";
const FLOWS_SCHEMA: &str = "frankensim-aperture-flows-v1";
/// Decimation filter length at the mechanical clock; odd, so its delay is exact.
const DECIMATION_TAPS: usize = 289;

fn fields(line: &str) -> Vec<&str> {
    line.split('\t').map(str::trim).collect()
}
fn num(s: &str, what: &str) -> Res<f64> {
    s.parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("{what}: not a finite number: {s:?}"))
}
fn read_text(path: &str) -> Res<String> {
    std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))
}
fn create(path: &str) -> Res<std::fs::File> {
    std::fs::File::create_new(path).map_err(|e| format!("{path}: {e}"))
}
fn json_str(s: &str) -> String {
    format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
}

struct Body {
    spec: LathedApertureSpec,
    medium: Medium,
}

fn parse_body(text: &str) -> Res<Body> {
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    if lines.next() != Some("frankensim-lathed-aperture-body-v1") {
        return Err("body: schema line frankensim-lathed-aperture-body-v1 required".into());
    }
    let (mut density, mut speed, mut around, mut step, mut fine) = (None, None, None, None, None);
    let (mut profile, mut wall, mut end) = (Vec::new(), Vec::new(), None);
    for line in lines {
        let f = fields(line);
        match (f[0], f.len()) {
            ("medium", 3) => {
                density = Some(num(f[1], "density")?);
                speed = Some(num(f[2], "sound speed")?);
            }
            ("around", 2) => around = Some(f[1].parse::<usize>().map_err(|e| e.to_string())?),
            ("axial_step_m", 2) => step = Some(num(f[1], "axial step")?),
            ("aperture_step_m", 2) => fine = Some(num(f[1], "aperture step")?),
            ("profile", 3) => profile.push((num(f[1], "profile x")?, num(f[2], "profile r")?)),
            ("wall", 5) => wall.push(WallAperture {
                name: f[1].to_string(),
                axial_m: num(f[2], "wall x")?,
                azimuth_rad: num(f[3], "wall azimuth")?,
                radius_m: num(f[4], "wall radius")?,
            }),
            ("end", 3) => {
                if end
                    .replace((f[1].to_string(), num(f[2], "end radius")?))
                    .is_some()
                {
                    return Err("body: one end aperture at most".into());
                }
            }
            _ => return Err(format!("body: unknown or malformed row {line:?}")),
        }
    }
    let missing = |what: &str| format!("body: {what} row required");
    Ok(Body {
        spec: LathedApertureSpec {
            profile,
            around: around.ok_or_else(|| missing("around"))?,
            axial_step_m: step.ok_or_else(|| missing("axial_step_m"))?,
            aperture_step_m: fine.ok_or_else(|| missing("aperture_step_m"))?,
            wall,
            end,
        },
        medium: Medium {
            density: density.ok_or_else(|| missing("medium"))?,
            sound_speed: speed.ok_or_else(|| missing("medium"))?,
        },
    })
}

fn parse_receivers(text: &str) -> Res<Vec<ReceiverSpec>> {
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    if lines.next() != Some("frankensim-receivers-v1") {
        return Err("receivers: schema line frankensim-receivers-v1 required".into());
    }
    lines
        .map(|line| {
            let f = fields(line);
            if f.len() != 5 || f[0] != "receiver" || f[1].is_empty() {
                return Err(format!("receivers: malformed row {line:?}"));
            }
            Ok(ReceiverSpec {
                id: f[1].to_string(),
                position_m: [num(f[2], "x")?, num(f[3], "y")?, num(f[4], "z")?],
            })
        })
        .collect()
}

/// Flag value after `name`, required.
fn flag<'a>(args: &'a [String], name: &str) -> Res<&'a str> {
    let i = args
        .iter()
        .position(|a| a == name)
        .ok_or_else(|| format!("{name} is required"))?;
    args.get(i + 1)
        .map(String::as_str)
        .ok_or_else(|| format!("{name} needs a value"))
}

const SAMPLES_SCHEMA: &str = "frankensim-observation-samples-v1";

fn write_samples(
    path: &str,
    operator: &str,
    fs_hz: f64,
    finest: usize,
    cache: &BTreeMap<usize, Vec<C64>>,
) -> Res<()> {
    let channels = cache.values().next().map_or(0, Vec::len);
    let mut bytes = format!(
        "{SAMPLES_SCHEMA}\noperator {operator}\nfs_hz {fs_hz:e}\nfinest {finest}\nchannels {channels}\nbins {}\ndata\n",
        cache.len()
    )
    .into_bytes();
    for (key, row) in cache {
        if row.len() != channels {
            return Err("samples: ragged rows".into());
        }
        bytes.extend_from_slice(&(*key as u64).to_le_bytes());
        for z in row {
            bytes.extend_from_slice(&z.re.to_le_bytes());
            bytes.extend_from_slice(&z.im.to_le_bytes());
        }
    }
    create(path)?
        .write_all(&bytes)
        .map_err(|e| format!("{path}: {e}"))
}

fn read_samples(
    path: &str,
    operator: &str,
    fs_hz: f64,
    finest: usize,
) -> Res<BTreeMap<usize, Vec<C64>>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)
        .and_then(|mut f| f.read_to_end(&mut bytes))
        .map_err(|e| format!("{path}: {e}"))?;
    let (header, mut at) = split_header(&bytes, "samples")?;
    let mut lines = header.lines();
    if lines.next() != Some(SAMPLES_SCHEMA) {
        return Err("samples: wrong schema".into());
    }
    let (mut op, mut fs, mut fin, mut channels, mut bins) = (String::new(), 0.0, 0, 0, 0);
    for line in lines {
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "operator" => op = f[1].to_string(),
            "fs_hz" => fs = num(f[1], "fs_hz")?,
            "finest" => fin = f[1].parse().map_err(|_| "samples: finest")?,
            "channels" => channels = f[1].parse().map_err(|_| "samples: channels")?,
            "bins" => bins = f[1].parse().map_err(|_| "samples: bins")?,
            _ => {}
        }
    }
    if op != operator || fs.to_bits() != fs_hz.to_bits() || fin != finest {
        return Err(format!("samples {path}: bound to another operator or grid"));
    }
    let mut cache = BTreeMap::new();
    for _ in 0..bins {
        let key = u64_at(&bytes, &mut at)? as usize;
        let row = (0..channels)
            .map(|_| {
                let re = f64::from_bits(u64_at(&bytes, &mut at)?);
                let im = f64::from_bits(u64_at(&bytes, &mut at)?);
                Ok(C64::new(re, im))
            })
            .collect::<Res<Vec<C64>>>()?;
        cache.insert(key, row);
    }
    if at != bytes.len() {
        return Err("samples: trailing bytes".into());
    }
    Ok(cache)
}

fn bake(args: &[String]) -> Res<()> {
    let [body_path, receivers_path, out] = [&args[0], &args[1], &args[2]];
    let fs_hz = num(flag(args, "--fs")?, "--fs")?;
    let pass = num(flag(args, "--pass")?, "--pass")?;
    let stop = num(flag(args, "--stop")?, "--stop")?;
    let latency: usize = flag(args, "--latency")?.parse().map_err(|_| "--latency")?;
    let clearance = num(flag(args, "--clearance")?, "--clearance")?;
    let threads: usize = flag(args, "--threads")?.parse().map_err(|_| "--threads")?;
    let fft: Vec<usize> = flag(args, "--fft")?
        .split(',')
        .map(|v| v.parse().map_err(|_| format!("--fft: {v}")))
        .collect::<Res<_>>()?;
    let body_text = read_text(body_path)?;
    let body = parse_body(&body_text)?;
    let receivers = parse_receivers(&read_text(receivers_path)?)?;
    let t0 = Instant::now();
    let mesh = lathed_aperture_body(&body.spec).map_err(|e| format!("mesh: {e}"))?;
    let sources: Vec<SourceId> = mesh
        .names
        .iter()
        .map(|n| SourceId::ApertureVolumeFlow(n.clone()))
        .collect();
    let obs = BemObservation::new(
        &mesh.surface,
        body.medium,
        &sources,
        &mesh.traces,
        &receivers,
        clearance,
    )
    .map_err(|e| format!("observation: {e}"))?;
    let decl = TransferDeclaration {
        producer_hash: obs.operator_hash(),
        sources: sources.clone(),
        receivers: receivers.clone(),
        fs_hz,
        geometric_delays_s: obs.geometric_delays_s(),
        pass_band_upper_hz: pass,
        taper_stop_hz: stop,
        filter_latency_samples: latency,
        fft_lengths: fft.clone(),
        tail_energy_tolerance: 1e-10,
        dtft_relative_tolerance: 1e-4,
        producer_uncertainty: ProducerUncertainty::Unresolved {
            source_hash: obs.operator_hash(),
        },
        maximum_frequency_entries: MAX_OBSERVATION_ENTRIES,
    };
    decl.admit_work().map_err(|e| format!("declaration: {e}"))?;
    let finest = *fft.last().ok_or("--fft is empty")?;
    // Every finest-grid bin below the taper stop; coarser rungs reuse them.
    let keys: Vec<usize> = (1..=finest / 2)
        .filter(|&b| fs_hz * b as f64 / (finest as f64) < stop)
        .collect();
    eprintln!(
        "panels {} apertures {} receivers {} solves {} threads {threads}",
        mesh.surface.areas().len(),
        sources.len(),
        receivers.len(),
        keys.len()
    );
    // Optional on-disk cache of the solved samples, bound to this operator and
    // grid, so a refused bake (band, latency, ladder) never repeats the solves.
    let samples_path = args
        .iter()
        .position(|a| a == "--samples")
        .and_then(|i| args.get(i + 1))
        .cloned();
    let operator_hex = obs.operator_hash().to_hex();
    let mut cache = match &samples_path {
        Some(path) if std::path::Path::new(path).exists() => {
            read_samples(path, &operator_hex, fs_hz, finest)?
        }
        _ => BTreeMap::new(),
    };
    let loaded = cache.len();
    let todo: Vec<usize> = keys
        .iter()
        .copied()
        .filter(|k| !cache.contains_key(k))
        .collect();
    eprintln!("cached samples {loaded}, to solve {}", todo.len());
    let next = std::sync::atomic::AtomicUsize::new(0);
    let samples: Vec<(usize, Res<Vec<C64>>)> = std::thread::scope(|s| {
        let workers: Vec<_> = (0..threads.max(1))
            .map(|_| {
                s.spawn(|| {
                    let mut done = Vec::new();
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        let Some(&key) = todo.get(i) else { break };
                        let omega = core::f64::consts::TAU * fs_hz * key as f64 / finest as f64;
                        let row = obs
                            .frequency_sample(omega)
                            .map_err(|e| format!("solve: {e}"));
                        eprintln!("bin {key} done at {:.0} s", t0.elapsed().as_secs_f64());
                        done.push((key, row));
                    }
                    done
                })
            })
            .collect();
        workers
            .into_iter()
            .flat_map(|w| w.join().expect("solver thread"))
            .collect()
    });
    let solved = samples.len();
    for (key, row) in samples {
        cache.insert(key, row?);
    }
    if let Some(path) = &samples_path
        && solved > 0
    {
        let target = if loaded > 0 {
            format!("{path}.extended")
        } else {
            path.clone()
        };
        write_samples(&target, &operator_hex, fs_hz, finest, &cache)?;
    }
    let baked = bake_observation_transfers(
        |omega| {
            let key = (omega / core::f64::consts::TAU * finest as f64 / fs_hz).round() as usize;
            cache
                .get(&key)
                .cloned()
                .ok_or(fs_bem::helmholtz::HelmholtzError::BadParameter {
                    what: "frequency outside the precomputed grid",
                })
        },
        &decl,
        || false,
    )
    .map_err(|e| format!("bake: {e}"))?;
    // Pressure kernels only; the gradient outputs are not replayed.
    let kernels: Vec<&FirKernel> = baked.kernels.iter().filter(|k| k.output % 4 == 0).collect();
    let mut file = std::io::BufWriter::new(create(&format!("{out}.kernels"))?);
    let mut header = format!(
        "{KERNELS_SCHEMA}\nfs_hz {fs_hz:e}\nlatency_samples {latency}\nsources {}\n",
        sources.len()
    );
    for n in &mesh.names {
        header += &format!("source {n}\n");
    }
    header += &format!("receivers {}\n", receivers.len());
    for r in &receivers {
        header += &format!(
            "receiver {} {:e} {:e} {:e}\n",
            r.id, r.position_m[0], r.position_m[1], r.position_m[2]
        );
    }
    header += &format!("kernels {}\ndata\n", kernels.len());
    let mut bytes = header.into_bytes();
    for k in &kernels {
        bytes.extend_from_slice(&(k.input as u64).to_le_bytes());
        bytes.extend_from_slice(&((k.output / 4) as u64).to_le_bytes());
        bytes.extend_from_slice(&(k.integer_delay_samples as u64).to_le_bytes());
        bytes.extend_from_slice(&(k.taps.len() as u64).to_le_bytes());
        for t in &k.taps {
            bytes.extend_from_slice(&t.to_le_bytes());
        }
    }
    file.write_all(&bytes)
        .and_then(|()| file.flush())
        .map_err(|e| format!("kernels write: {e}"))?;
    let ladder: Vec<String> = baked
        .ladder
        .iter()
        .map(|r| {
            format!(
                "{{\"n_fft\":{},\"tail_energy_fraction\":{:e},\"wrap_energy_fraction\":{:e},\"dtft_relative_change\":{}}}",
                r.n_fft,
                r.tail_energy_fraction,
                r.wrap_energy_fraction,
                r.dtft_relative_change.map_or("null".into(), |v| format!("{v:e}"))
            )
        })
        .collect();
    let apertures: Vec<String> = mesh
        .names
        .iter()
        .zip(&mesh.patch_areas_m2)
        .map(|(n, a)| format!("{{\"name\":{},\"patch_area_m2\":{a:e}}}", json_str(n)))
        .collect();
    let report = format!(
        "{{\"schema\":\"frankensim-aperture-reobservation-bake-v1\",\"body_blake3\":\"{}\",\"operator_blake3\":\"{}\",\"declaration_blake3\":\"{}\",\"panels\":{},\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e},\"fs_hz\":{fs_hz:e},\"pass_band_upper_hz\":{pass:e},\"taper_stop_hz\":{stop:e},\"filter_latency_samples\":{latency},\"solves\":{},\"seconds\":{:.1},\"producer_uncertainty\":\"unresolved: compare panel-ladder rungs\",\"apertures\":[{}],\"ladder\":[{}]}}\n",
        fs_blake3::hash_domain("frankensim-lathed-aperture-body-v1", body_text.as_bytes()).to_hex(),
        obs.operator_hash().to_hex(),
        baked.declaration_hash.to_hex(),
        mesh.surface.areas().len(),
        body.medium.density,
        body.medium.sound_speed,
        keys.len(),
        t0.elapsed().as_secs_f64(),
        apertures.join(","),
        ladder.join(",")
    );
    create(&format!("{out}.json"))?
        .write_all(report.as_bytes())
        .map_err(|e| format!("report write: {e}"))?;
    println!("{}", report.trim_end());
    Ok(())
}

/// Split a "...\ndata\n" text header from its binary payload.
fn split_header(bytes: &[u8], what: &str) -> Res<(String, usize)> {
    let end = bytes
        .windows(6)
        .position(|w| w == b"\ndata\n")
        .ok_or_else(|| format!("{what}: no data line"))?
        + 6;
    let header =
        std::str::from_utf8(&bytes[..end]).map_err(|_| format!("{what}: header is not text"))?;
    Ok((header.to_string(), end))
}
fn u64_at(bytes: &[u8], at: &mut usize) -> Res<u64> {
    let v = bytes
        .get(*at..*at + 8)
        .ok_or("kernels: truncated")?
        .try_into()
        .map_err(|_| "kernels: truncated")?;
    *at += 8;
    Ok(u64::from_le_bytes(v))
}

fn wav_f32(path: &str, rate: u32, samples: &[f64]) -> Res<()> {
    let data = (samples.len() * 4) as u32;
    let mut b = Vec::with_capacity(44 + samples.len() * 4);
    b.extend_from_slice(b"RIFF");
    b.extend_from_slice(&(36 + data).to_le_bytes());
    b.extend_from_slice(b"WAVEfmt ");
    b.extend_from_slice(&16u32.to_le_bytes());
    b.extend_from_slice(&3u16.to_le_bytes());
    b.extend_from_slice(&1u16.to_le_bytes());
    b.extend_from_slice(&rate.to_le_bytes());
    b.extend_from_slice(&(rate * 4).to_le_bytes());
    b.extend_from_slice(&4u16.to_le_bytes());
    b.extend_from_slice(&32u16.to_le_bytes());
    b.extend_from_slice(b"data");
    b.extend_from_slice(&data.to_le_bytes());
    for &s in samples {
        b.extend_from_slice(&(s as f32).to_le_bytes());
    }
    create(path)?
        .write_all(&b)
        .map_err(|e| format!("{path}: {e}"))
}

fn replay(args: &[String]) -> Res<()> {
    let [kernels_path, flows_path, bind_path, out_dir] = [&args[0], &args[1], &args[2], &args[3]];
    let mut kbytes = Vec::new();
    std::fs::File::open(kernels_path)
        .and_then(|mut f| f.read_to_end(&mut kbytes))
        .map_err(|e| format!("{kernels_path}: {e}"))?;
    let (kheader, mut at) = split_header(&kbytes, "kernels")?;
    let mut lines = kheader.lines();
    if lines.next() != Some(KERNELS_SCHEMA) {
        return Err("kernels: wrong schema".into());
    }
    let (mut fs_hz, mut latency, mut sources, mut receivers, mut count) =
        (0.0, 0usize, Vec::new(), Vec::new(), 0usize);
    for line in lines {
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "fs_hz" => fs_hz = num(f[1], "fs_hz")?,
            "latency_samples" => latency = f[1].parse().map_err(|_| "latency")?,
            "source" => sources.push(f[1].to_string()),
            "receiver" => receivers.push(f[1].to_string()),
            "kernels" => count = f[1].parse().map_err(|_| "kernels")?,
            _ => {}
        }
    }
    let mut kernels = Vec::with_capacity(count);
    for _ in 0..count {
        let input = u64_at(&kbytes, &mut at)? as usize;
        let output = u64_at(&kbytes, &mut at)? as usize;
        let integer_delay_samples = u64_at(&kbytes, &mut at)? as usize;
        let n = u64_at(&kbytes, &mut at)? as usize;
        let taps = (0..n)
            .map(|_| u64_at(&kbytes, &mut at).map(f64::from_bits))
            .collect::<Res<Vec<f64>>>()?;
        kernels.push(FirKernel {
            input,
            output,
            integer_delay_samples,
            taps,
        });
    }
    if at != kbytes.len() {
        return Err("kernels: trailing bytes".into());
    }
    // Flows: header names each aperture by role and node, in data order.
    let mut fbytes = Vec::new();
    std::fs::File::open(flows_path)
        .and_then(|mut f| f.read_to_end(&mut fbytes))
        .map_err(|e| format!("{flows_path}: {e}"))?;
    let (fheader, fstart) = split_header(&fbytes, "flows")?;
    let mut flines = fheader.lines();
    if flines.next() != Some(FLOWS_SCHEMA) {
        return Err("flows: wrong schema".into());
    }
    let (mut mech_hz, mut nodes) = (0.0, Vec::new());
    for line in flines {
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "mechanics_rate_hz" => mech_hz = num(f[1], "mechanics rate")?,
            "aperture" => nodes.push(f[3].to_string()),
            _ => {}
        }
    }
    let ratio = mech_hz / fs_hz;
    if !(ratio >= 1.0 && ratio.fract() == 0.0) {
        return Err(format!(
            "flows at {mech_hz} Hz do not decimate to {fs_hz} Hz by an integer"
        ));
    }
    let ratio = ratio as usize;
    let data = &fbytes[fstart..];
    let width = nodes.len();
    if width == 0 || data.len() % (8 * width) != 0 {
        return Err("flows: payload is not whole frames".into());
    }
    let frames = data.len() / (8 * width);
    // BIND: "bind <node> <source>" maps a tapped node to a body aperture.
    let mut column_of_source = vec![None; sources.len()];
    for line in read_text(bind_path)?
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
    {
        let f = fields(line);
        if f.len() != 3 || f[0] != "bind" {
            return Err(format!("bind: malformed row {line:?}"));
        }
        let column = nodes
            .iter()
            .position(|n| n == f[1])
            .ok_or_else(|| format!("bind: node {} is not in the flow file", f[1]))?;
        let source = sources
            .iter()
            .position(|s| s == f[2])
            .ok_or_else(|| format!("bind: aperture {} is not in the body", f[2]))?;
        if column_of_source[source].replace(column).is_some() {
            return Err(format!("bind: aperture {} bound twice", f[2]));
        }
    }
    if column_of_source.iter().flatten().count() != width {
        return Err("bind: every tapped aperture must be bound exactly once".into());
    }
    // Blackman-windowed sinc, cutoff at 0.42 of the output rate's band.
    let cutoff = 0.42 * fs_hz / mech_hz;
    let half = (DECIMATION_TAPS / 2) as f64;
    let mut lowpass: Vec<f64> = (0..DECIMATION_TAPS)
        .map(|i| {
            let t = i as f64 - half;
            let sinc = if t == 0.0 {
                2.0 * cutoff
            } else {
                (core::f64::consts::TAU * cutoff * t).sin() / (core::f64::consts::PI * t)
            };
            let w = core::f64::consts::TAU * i as f64 / (DECIMATION_TAPS - 1) as f64;
            sinc * (0.42 - 0.5 * w.cos() + 0.08 * (2.0 * w).cos())
        })
        .collect();
    let gain: f64 = lowpass.iter().sum();
    lowpass.iter_mut().for_each(|v| *v /= gain);
    let flow = |frame: usize, column: usize| {
        let o = 8 * (frame * width + column);
        f64::from_le_bytes(data[o..o + 8].try_into().expect("8 bytes"))
    };
    let outputs = frames / ratio;
    let mut fir = ObservationFir::new(sources.len(), receivers.len(), kernels)
        .map_err(|e| format!("observer: {e}"))?;
    let mut pressure = vec![Vec::with_capacity(outputs); receivers.len()];
    let mut input = vec![0.0; sources.len()];
    for n in 0..outputs {
        let centre = n * ratio;
        for (s, column) in column_of_source.iter().enumerate() {
            input[s] = column.map_or(0.0, |c| {
                lowpass
                    .iter()
                    .enumerate()
                    .filter_map(|(i, h)| {
                        (centre + i)
                            .checked_sub(DECIMATION_TAPS / 2)
                            .filter(|&f| f < frames)
                            .map(|f| h * flow(f, c))
                    })
                    .sum()
            });
        }
        let y = fir.step(&input).map_err(|e| format!("replay: {e}"))?;
        for (r, v) in y.into_iter().enumerate() {
            pressure[r].push(v);
        }
    }
    std::fs::create_dir(out_dir).map_err(|e| format!("{out_dir}: {e}"))?;
    let mut rows = Vec::new();
    for (id, p) in receivers.iter().zip(&pressure) {
        wav_f32(&format!("{out_dir}/{id}.wav"), fs_hz as u32, p)?;
        let rms = (p.iter().map(|v| v * v).sum::<f64>() / p.len().max(1) as f64).sqrt();
        rows.push(format!(
            "{{\"receiver\":{},\"rms_pa\":{rms:e}}}",
            json_str(id)
        ));
    }
    let report = format!(
        "{{\"schema\":\"frankensim-aperture-reobservation-replay-v1\",\"kernels\":{},\"flows\":{},\"mechanics_rate_hz\":{mech_hz:e},\"output_rate_hz\":{fs_hz:e},\"frames\":{frames},\"output_samples\":{outputs},\"decimation\":\"Blackman-windowed sinc, {DECIMATION_TAPS} taps at the mechanics clock, cutoff 0.42 fs_out, unit DC gain, centred (zero phase); first output sample is mechanical frame 0\",\"observation_latency_samples\":{latency},\"receivers\":[{}]}}\n",
        json_str(kernels_path),
        json_str(flows_path),
        rows.join(",")
    );
    create(&format!("{out_dir}/replay.json"))?
        .write_all(report.as_bytes())
        .map_err(|e| format!("report: {e}"))?;
    println!("{}", report.trim_end());
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("bake") if args.len() >= 4 => bake(&args[1..]),
        Some("replay") if args.len() == 5 => replay(&args[1..]),
        _ => Err("usage: reobserve_apertures bake BODY RECEIVERS OUT --fs HZ --pass HZ --stop HZ --latency N --fft N1,N2,N3,N4 --clearance M --threads N [--samples CACHE]\n       reobserve_apertures replay KERNELS FLOWS BIND OUTDIR".into()),
    };
    if let Err(e) = result {
        eprintln!("reobserve_apertures: {e}");
        std::process::exit(2);
    }
}
