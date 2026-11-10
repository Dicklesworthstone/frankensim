//! File-driven lip-valve performances: `music_render lip INPUT.performance OUT.wav`.
//!
//! The mechanics clock of a bore with a large radiating mouth is bounded by the
//! flanged-radiation ceiling (`ka <= 20` at Nyquist), so the file declares its
//! own clock and this command writes the WAV at that clock. Nothing is
//! resampled; the sidecar states the clock and the declared acoustic band.
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use fs_blake3::DomainHasher;
use fs_couple::brass_loop::performance::{
    LipBoreInfo, LipDescription, LipInfo, MAX_LIP_PERFORMANCE_BYTES,
};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::stream::{MAX_PCM16_WAV_SAMPLES, Pcm16WavStream};
use fs_material::gas_species::input::ThermochemicalSources;

use super::{WAV_HASH_DOMAIN, create_outputs, json_string};

type Options<'a> = (&'a str, &'a str, usize, Option<&'a str>, Option<&'a str>);

fn options(args: &[String]) -> Result<Options<'_>, String> {
    let mut paths = Vec::new();
    let mut block = None;
    let mut mouthpiece = None;
    let mut impedance = None;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--block" => {
                if block.is_some() {
                    return Err("--block may only be supplied once".into());
                }
                block = Some(
                    iter.next()
                        .and_then(|v| v.parse::<usize>().ok())
                        .filter(|v| (1..=65_536).contains(v))
                        .ok_or_else(|| "--block needs an integer in 1..=65536".to_string())?,
                );
            }
            "--impedance" => {
                if impedance.is_some() {
                    return Err("--impedance may only be supplied once".into());
                }
                impedance = Some(
                    iter.next()
                        .map(String::as_str)
                        .filter(|v| !v.starts_with('-'))
                        .ok_or_else(|| "--impedance needs a path".to_string())?,
                );
            }
            "--mouthpiece-wav" => {
                if mouthpiece.is_some() {
                    return Err("--mouthpiece-wav may only be supplied once".into());
                }
                mouthpiece = Some(
                    iter.next()
                        .map(String::as_str)
                        .filter(|v| !v.starts_with('-'))
                        .ok_or_else(|| "--mouthpiece-wav needs a path".to_string())?,
                );
            }
            value if value.starts_with('-') => {
                return Err(format!(
                    "unsupported lip option {value:?}; geometry, lip, pressure, duration and full-scale belong to the source file"
                ));
            }
            value => paths.push(value),
        }
    }
    let [input, output] = paths.as_slice() else {
        return Err(
            "usage: music_render lip INPUT.performance OUT.wav [--block N] [--mouthpiece-wav PATH] [--impedance PATH.tsv]"
                .into(),
        );
    };
    Ok((*input, *output, block.unwrap_or(512), mouthpiece, impedance))
}

struct Rendered {
    clipped: u64,
    peak_pa: f64,
    rms_pa: f64,
    hash: fs_blake3::ContentHash,
}

// The shared command writer pins 48 kHz; this is the same staging, hashing and
// never-normalize contract on the source's declared clock.
fn render_at_rate<W: Read + Write + Seek>(
    renderer: &mut impl PressureRenderer,
    output: &mut W,
    rate: u32,
    samples: usize,
    block: usize,
    full_scale_pa: f64,
) -> Result<Rendered, String> {
    if samples == 0
        || samples as u128 > u128::from(MAX_PCM16_WAV_SAMPLES)
        || block == 0
        || block > renderer.max_block_len()
    {
        return Err(
            "render output needs a positive RIFF-sized history and an admitted block".into(),
        );
    }
    renderer
        .validate_sample_rate(rate)
        .map_err(|e| e.to_string())?;
    if renderer.samples_rendered() != 0 {
        return Err("the command output starts at sample zero".into());
    }
    renderer
        .validate_sample_count(samples as u64)
        .map_err(|e| e.to_string())?;
    let buffered = BufWriter::new(&mut *output);
    let mut stream =
        Pcm16WavStream::new(buffered, rate, full_scale_pa, block).map_err(|e| e.to_string())?;
    let mut pressure = vec![0.0_f64; block.min(samples)];
    let (mut peak, mut squares, mut completed) = (0.0_f64, 0.0_f64, 0);
    while completed < samples {
        let count = pressure.len().min(samples - completed);
        let chunk = &mut pressure[..count];
        renderer
            .block(chunk)
            .map_err(|e| format!("render refused: {e}"))?;
        for &value in chunk.iter() {
            peak = peak.max(value.abs());
            squares += value * value;
        }
        if !peak.is_finite() || !squares.is_finite() {
            return Err(
                "pressure statistics exceeded the finite range; output is incomplete".into(),
            );
        }
        stream.write_block(chunk).map_err(|e| e.to_string())?;
        completed += count;
    }
    let (buffered, summary) = stream.finish().map_err(|e| e.to_string())?;
    drop(buffered);
    output
        .seek(SeekFrom::Start(0))
        .map_err(|e| format!("WAV hash seek failed: {e}"))?;
    let mut remaining = 44 + 2 * summary.samples;
    let mut bytes = [0_u8; 8192];
    let mut hasher = DomainHasher::new(WAV_HASH_DOMAIN);
    while remaining > 0 {
        let count = remaining.min(bytes.len() as u64) as usize;
        output
            .read_exact(&mut bytes[..count])
            .map_err(|e| format!("WAV hash read failed: {e}"))?;
        hasher.update(&bytes[..count]);
        remaining -= count as u64;
    }
    Ok(Rendered {
        clipped: summary.clipped_samples,
        peak_pa: peak,
        rms_pa: (squares / samples as f64).sqrt(),
        hash: hasher.finalize(),
    })
}

fn number(x: f64) -> String {
    if x.is_finite() {
        format!("{x:e}")
    } else {
        "null".to_string()
    }
}

// The bore engine's receipts and the scope of what it claims.
fn bore_json(bore: &LipBoreInfo) -> String {
    match bore {
        LipBoreInfo::Linear(lines) => {
            let ladder = lines
                .mouth_flow_ladder
                .iter()
                .map(|(taps, tail, change)| {
                    format!(
                        "{{\"taps\":{taps},\"tail_energy\":{},\"tap_change\":{}}}",
                        number(*tail),
                        number(*change)
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            let line = &lines.line;
            // A player's tract admitted by the passivity gain names it (MB 20665); every other card keeps its provenance bytes.
            let tract_gain = lines.player_tract_passivity_gain.map_or_else(String::new, |gain| {
                format!(
                    "\"player_tract_passivity_gain\":{{\"passes\":{},\"top_band_db\":{:e}}},",
                    gain.passes,
                    gain.top_band_db
                )
            });
            let (radiation, scope) = if let Some(r) = &lines.modal_radiation {
                let changes = |values: &[Vec<f64>; 2]| values.iter().map(|row|
                    format!("[{}]", row.iter().map(|v| number(*v)).collect::<Vec<_>>().join(",")))
                    .collect::<Vec<_>>().join(",");
                (format!("{{\"modes\":{},\"radial_nodes\":{:?},\"angular_points\":{:?},\"radial_kernel_changes\":[{}],\"angular_kernel_changes\":[{}],\"basis_gram_defect\":{},\"kernel_terms\":{},\"tolerance\":{}}}",
                    r.modes, r.radial_nodes, r.angular_points, changes(&r.radial_changes),
                    changes(&r.angular_changes), number(r.basis_gram_defect), r.kernel_terms, number(r.tolerance)),
                 "full flanged modal reflection load and finite-distance Rayleigh observation of every original signed axisymmetric mouth mode, with the same gas and aperture; infinite rigid baffle; omitted modes and finite-body scattering unresolved")
            } else {
                ("null".into(), "uniform baffled-disk Rayleigh radiation of total mouth flow; higher mouth-mode directivity unresolved; source termination retains its declared scalar load")
            };
            format!(
                "{tract_gain}\"reflectance_line\":{{\"taps\":{},\"round_trip_samples\":{},\"peak_reflectance\":{},\"rungs_tried\":{},\"tail_energy\":{},\"peak_shift_cents\":{},\
                 \"image_peaks\":{},\"image_worst_cents\":{},\"image_worst_db\":{}}},\
                 \"mouth_flow_line\":{{\"taps\":{},\"rungs\":{},\"tail_tol\":{},\"ladder\":[{}],\"modes\":{}}},\
                 \"modal_radiation\":{},\
                 \"model_scope\":\"one held bore; 1-DOF outward-striking lip island; plane-source multimodal reflectance and mouth flow from the same retained modal bore; {}; no valve or slide motion, no bore wall vibration, no breath heating, no nonlinear propagation; spatial and temporal convergence do not establish instrument realism\"",
                line.n_fft,
                number(line.round_trip_samples),
                number(line.peak_reflectance),
                line.rungs_tried,
                number(line.tail_energy),
                number(line.peak_shift_cents),
                line.image_peaks,
                number(line.image_worst_cents),
                number(line.image_worst_db),
                lines.mouth_flow_taps,
                lines.mouth_flow_rungs,
                number(lines.mouth_flow_tail_tol),
                ladder,
                lines.mouth_flow_modes,
                radiation,
                scope,
            )
        }
        LipBoreInfo::Nonlinear(n) => format!(
            "\"nonlinear_bore\":{{\"maximum_cell_m\":{},\"cells\":{},\"mechanics_rate_hz\":{},\"gas_steps_per_mechanics\":{},\"gas_step_s\":{},\"radiation_taps\":{}{},\
             \"wall\":{{\"poles\":{},\"band_hz\":[{},{}],\"tolerance\":{},\"distinct_radii\":{},\"maximum_arms\":{},\"maximum_complex_error\":{},\"maximum_loss_error\":{}}}}},\
             \"model_scope\":\"one held bore; 1-DOF outward-striking lip island metering an isentropic compressible lip jet that mixes out at the cup face; conservative nonlinear Euler gas of the ambient NASA-9 mixture on the declared cells (steepening and shocks resolved by those cells); Zwikker-Kosten wall memory admitted on its declared band; rigid baffled-piston bell radiation; baffled-disk Rayleigh radiation of the bell volume flow at the declared point; no valve or slide motion, no bore wall vibration, no breath heating, no transverse duct modes\"",
            number(n.record.maximum_cell_m),
            n.cells,
            n.mechanics_rate_hz,
            n.record.gas_steps_per_mechanics,
            number(n.gas_step_s),
            n.radiation_taps,
            mouth_load_json(n.mouth_load.as_ref()),
            n.record.wall_poles,
            number(n.record.wall_band_hz[0]),
            number(n.record.wall_band_hz[1]),
            number(n.record.wall_tolerance),
            n.wall_distinct_radii,
            n.wall_maximum_arms,
            number(n.wall_maximum_complex_error),
            number(n.wall_maximum_loss_error),
        ),
    }
}

/// The `mouth_load modal` kernel's receipt object, with its leading comma; empty for the baffled piston, so every receipt
/// written without the row is unchanged.
fn mouth_load_json(info: Option<&fs_couple::brass_loop::performance::MouthLoadInfo>) -> String {
    match info {
        None => String::new(),
        Some(m) => format!(
            ",\"mouth_load\":{{\"modes\":{},\"extra_slices\":{},\"band_hz\":{},\"bin_hz\":{},\"tail_tolerance\":{},\"taps\":{},\"mouth_radius_m\":{},\"z_inf\":{},\
             \"dropped_tail_energy\":{},\"dropped_sum\":{},\"evaluated_hz\":{},\"minimum_resistance_ratio\":{},\
             \"scope\":\"the bore's radiating end is the plane-mode equivalent of the complete flanged multimodal mouth (lossless modal solve); the receiver is still the baffled-piston Rayleigh integral of the plane mouth flow\"}}",
            m.record.n_modes,
            m.record.extra_slices,
            number(m.record.band_hz),
            number(m.bin_hz),
            number(m.record.tail_tolerance),
            m.taps,
            number(m.mouth_radius_m),
            number(m.z_inf),
            number(m.dropped_tail_energy),
            number(m.dropped_sum),
            number(m.evaluated_hz),
            number(m.minimum_resistance_ratio),
        ),
    }
}

/// The lip receipt's optional second-mode object, empty for a performance without the `lip_response_mode2` record (so
/// every receipt written before the record existed is unchanged).
fn lip_mode2_json(i: &LipInfo) -> String {
    match &i.lip_mode2 {
        None => String::new(),
        Some(m) => format!(
            "\"lip_mode2\":{{\"mass_kg\":{},\"stiffness_n_m\":{},\"damping_n_s_m\":{},\"sign\":{},\"model\":\"this lip has TWO modes: the lip record above is the outward-striking first mode and this is the inward-striking second (aperture = rest gap + q1 + q2, q2 driven by sign x the same force); the model_scope's 1-DOF wording describes the lip without this record\"}},",
            number(m.mass_kg),
            number(m.stiffness_n_m),
            number(m.damping_n_s_m),
            number(m.sign),
        ),
    }
}

/// The lip receipt's optional flow-area-exponent object, empty for a performance without the `lip_area_exponent` record (so
/// every receipt written before the record existed is unchanged).
fn lip_area_exponent_json(i: &LipInfo) -> String {
    match i.lip_area_exponent {
        None => String::new(),
        Some(q) => format!(
            "\"lip_area_exponent\":{{\"q\":{},\"model\":\"flow aperture = width x gap x (gap / rest gap)^(q-1), a stated bracket (Copley and Strong 1996 call lip area against height almost linear and print no exponent); the normalisation at the rest gap is ours; the lip's motion is unchanged\"}},",
            number(q),
        ),
    }
}

/// The lip receipt's optional rim-stop object, empty for a performance without the `lip_rim_stop` record (so every receipt written
/// before the record existed is unchanged).
fn rim_stop_json(i: &LipInfo) -> String {
    match i.rim_stop {
        None => String::new(),
        Some(stop) => format!(
            "\"rim_stop\":{{\"kind\":\"{}\",\"max_gap_m\":{},\"rule\":\"the lip gap may not exceed twice the entry radius of the first bore cone (the mouthpiece rim opening); the first lip mode's displacement is held at the bound and its outward momentum is removed (inelastic) or reversed (elastic)\",\"restitution_measured\":false,\"note\":\"the two kinds are limit cases of a restitution no source measures; neither is a measured contact law\"}},",
            stop.kind.label(),
            number(stop.max_gap_m),
        ),
    }
}

fn input_json(i: &LipInfo) -> String {
    let lip = &i.lip;
    let mic = &i.microphone;
    let second = lip_mode2_json(i);
    let exponent = lip_area_exponent_json(i);
    let rim_stop = rim_stop_json(i);
    let observer_gas = i.observer_gas.as_ref().map_or_else(String::new, |report| {
        format!("\"observer_gas\":{},",
            super::duct_metadata::observer_gas_provenance(Some(report)))
    });
    format!(
        "\"lip_input\":{{\"schema\":\"{}\",\"blake3\":\"{}\",\"bore\":{{\"segments\":{},\"length_m\":{},\"mouth_radius_m\":{},\"gas_temperature_k\":{},\"sound_speed_m_s\":{},\"source\":{}}},\
         \"lip\":{{\"mass_kg\":{},\"stiffness_n_m\":{},\"damping_n_s_m\":{},\"width_m\":{},\"rest_gap_m\":{},\"face_area_m2\":{},\"source\":{}}},{}{}{}\
         \"cup_volume_m3\":{},\"control_rate_hz\":{},\"lip_tension_track\":{},\
         \"microphone\":{{\"position_m\":[{},{},{}],\"radial_rings\":{},\"angular_points\":{},\"resolved_band_hz\":{}}},{observer_gas}\
         {}}}",
        fs_couple::brass_loop::performance::LIP_PERFORMANCE_SCHEMA,
        i.input_hash.to_hex(),
        i.segments,
        number(i.bore_length_m),
        number(i.mouth_radius_m),
        number(i.gas_temperature_k),
        number(i.sound_speed_m_s),
        json_string(&i.bore_source),
        number(lip.mass_kg),
        number(lip.stiffness_n_m),
        number(lip.damping_n_s_m),
        number(lip.width_m),
        number(lip.rest_gap_m),
        number(lip.face_area_m2),
        json_string(&i.lip_source),
        second,
        exponent,
        rim_stop,
        number(i.cup_volume_m3),
        i.control_rate_hz,
        i.has_lip_tension_track,
        number(mic.position_m[0]),
        number(mic.position_m[1]),
        number(mic.position_m[2]),
        mic.radial_rings,
        mic.angular_points,
        number(mic.maximum_frequency_hz),
        bore_json(&i.bore),
    )
}

pub(super) fn run_with_gas_sources(
    args: &[String],
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<(), String> {
    let (input, output, block, mouthpiece_path, impedance_path) = options(args)?;
    let output = Path::new(output);
    let sidecar = output.with_extension("provenance.json");
    if output == sidecar.as_path() || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let file = std::fs::File::open(input).map_err(|e| format!("lip input open failed: {e}"))?;
    let mut bytes = Vec::new();
    file.take((MAX_LIP_PERFORMANCE_BYTES + 1) as u64)
        .read_to_end(&mut bytes)
        .map_err(|e| format!("lip input read failed: {e}"))?;
    let description = LipDescription::from_bytes_with_gas_sources(&bytes, sources)?;
    if let Some(path) = impedance_path {
        // Passive-bore scan, 20 Hz to 3 kHz in 1 Hz steps; no loop, lines or render are built and OUT.wav is not created.
        let frequencies: Vec<f64> = (20..=3000).map(f64::from).collect();
        let scan = description.bore_scan(&frequencies)?;
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .map_err(|e| format!("impedance table create refused: {e}"))?;
        let mut text = format!(
            "# bore scan of {input}: zc_in {:e} Pa s/m3, rho {:e} kg/m3, microphone distance {:e} m\n\
             f_hz\tzmm_re\tzmm_im\tzsc_re\tzsc_im\tflow_re\tflow_im\tradiated_per_input_pa_per_pa\n",
            scan.input_characteristic_impedance, scan.density_kg_m3, scan.microphone_distance_m
        );
        for r in &scan.rows {
            let z = r.impedance_scalar[0].hypot(r.impedance_scalar[1]);
            let flow = r.mouth_over_input_flow[0].hypot(r.mouth_over_input_flow[1]);
            // far field of a baffled piston on axis: p = rho omega Q / (2 pi r), Q = mouth flow = flow_ratio * p_in / Z_in
            let radiated = scan.density_kg_m3 * core::f64::consts::TAU * r.frequency_hz * flow
                / (core::f64::consts::TAU * scan.microphone_distance_m * z);
            text.push_str(&format!(
                "{}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n",
                r.frequency_hz,
                r.impedance_multimodal[0],
                r.impedance_multimodal[1],
                r.impedance_scalar[0],
                r.impedance_scalar[1],
                r.mouth_over_input_flow[0],
                r.mouth_over_input_flow[1],
                radiated
            ));
        }
        file.write_all(text.as_bytes())
            .and_then(|()| file.flush())
            .map_err(|e| format!("impedance table write failed: {e}"))?;
        println!(
            "{{\"suite\":\"music-render\",\"impedance_table\":{},\"rows\":{}}}",
            json_string(path),
            scan.rows.len()
        );
        return Ok(());
    }
    let mut performance = description.build()?;
    let info = performance.info().clone();
    let samples =
        usize::try_from(info.samples).map_err(|_| "sample count exceeds this host".to_string())?;
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    if mouthpiece_path.is_some() {
        performance.record_mouthpiece();
    }
    let rendered = render_at_rate(
        &mut performance,
        &mut audio,
        info.sample_rate_hz,
        samples,
        block,
        info.full_scale_pa,
    )?;
    let wav_hash = rendered.hash.to_hex();
    let (mechanics_rate_hz, ratio) = match &info.bore {
        LipBoreInfo::Linear(_) => (info.sample_rate_hz, 1),
        LipBoreInfo::Nonlinear(n) => (n.mechanics_rate_hz, n.record.mechanics_per_output),
    };
    let windows_json = performance
        .windows()
        .iter()
        .map(|w| {
            format!(
                "{{\"start_s\":{},\"bore_mean_pa\":{},\"bore_ac_rms_pa\":{},\"bore_ac_peak_pa\":{},\"gap_min_m\":{},\"gap_mean_m\":{},\"gap_max_m\":{},\"closed_chunk_fraction\":{},\"lip_tension_work_j\":[{},{}]}}",
                number(w.start_s),
                number(w.bore_mean_pa),
                number(w.bore_ac_rms_pa),
                number(w.bore_ac_peak_pa),
                number(w.gap_min_m),
                number(w.gap_mean_m),
                number(w.gap_max_m),
                number(w.closed_chunk_fraction),
                number(w.lip_tension_work_j[0]),
                number(w.lip_tension_work_j[1]),
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let gas_sources =
        super::thermochemical::provenance(core::slice::from_ref(description.gas_record()));
    // Only a performance with the `bore_gas` row writes the two gases it separates (the bore's and the room's complete
    // states); every receipt written before the row existed is unchanged.
    let bore_gas = description.bore_gas_record().map_or_else(String::new, |record| {
        format!(
            ",\"bore_gas\":{{\"bore\":{},\"room\":{}}}",
            super::duct_metadata::gas_state(record.state()),
            super::duct_metadata::gas_state(description.ambient_gas()),
        )
    });
    // Nonlinear bore only: gas steps accepted after a first-order retry of the failed second-order update.
    let fallbacks = performance.first_order_fallbacks().map_or_else(String::new, |n| {
        format!(",\"nonlinear_first_order_fallback_steps\":{n}")
    });
    // Only a performance with the `lip_rim_stop` row writes the samples at which its stop acted.
    let rim_hits = performance.rim_stop_hits().map_or_else(String::new, |n| format!(",\"rim_stop_hits\":{n}"));
    // Only a performance with the `jet_noise` row writes noise-band coverage.
    let jet_band = performance.jet_noise_coverage().map_or_else(String::new, |c| {
        format!(
            "\"jet_noise_band\":{{\"open_gate_steps\":{},\"whole\":{},\"truncated\":{},\"silent\":{}}},",
            c.open_gate(),
            c.whole,
            c.truncated,
            c.silent
        )
    });
    let tension_work = performance.lip_tension_work_j();
    let tension_work_json = format!(
        "\"lip_tension_work_j\":[{},{}]",
        number(tension_work[0]),
        number(tension_work[1]),
    );
    let (observation_filter, observation_delay) = if info.observer_gas.is_some() {
        (
            "shared differentiating Rayleigh receiver with material gas propagation on every surface path",
            format!(",\"filter_delay_output_samples\":{}", number(performance.filter_delay_samples())),
        )
    } else {
        (
            "none: the baffled receiver differentiates and delays the mouth flow of every mechanical sample",
            String::new(),
        )
    };
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"lip-input\",\
         \"sample_rate_hz\":{},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",{},\
         \"lip_windows_0p1s\":[{}],{tension_work_json},{jet_band}\
         \"observation\":{{\"mechanics_sample_rate_hz\":{mechanics_rate_hz},\"output_sample_rate_hz\":{},\"ratio\":{ratio},\
         \"filter\":\"{observation_filter}\",\
         \"resampled\":false,\"declared_band_hz\":{}{observation_delay}}}{fallbacks}{rim_hits}{gas_sources}{bore_gas}}}",
        info.sample_rate_hz,
        info.full_scale_pa,
        rendered.clipped,
        rendered.peak_pa,
        rendered.rms_pa,
        input_json(&info),
        windows_json,
        info.sample_rate_hz,
        number(info.microphone.maximum_frequency_hz),
    );
    writeln!(metadata, "{provenance}")
        .and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
    if let Some(path) = mouthpiece_path {
        // Diagnostic output: the INTERNAL pressure at the mouthpiece, mean (the mouth's DC operating point) removed. The
        // PCM scale is the next power of ten above its peak, printed here because the render's own scale is far too small.
        let mut signal = performance.mouthpiece_pressure().to_vec();
        let mean = signal.iter().sum::<f64>() / signal.len().max(1) as f64;
        for p in &mut signal {
            *p -= mean;
        }
        let peak = signal.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
        let full_scale = 10f64.powf((peak * 1.01).max(1.0).log10().ceil());
        let (wav, clipped) =
            fs_couple::pcm_wav::encode_pcm16_wav(&signal, info.sample_rate_hz, full_scale)
                .map_err(|e| format!("mouthpiece wav encode failed: {e:?}"))?;
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .map_err(|e| format!("mouthpiece wav create refused: {e}"))?;
        file.write_all(&wav)
            .and_then(|()| file.flush())
            .map_err(|e| format!("mouthpiece wav write failed: {e}"))?;
        println!(
            "{{\"suite\":\"music-render\",\"mouthpiece_wav\":{},\"full_scale_pa\":{full_scale:e},\"dc_removed_pa\":{mean:e},\"peak_ac_pa\":{peak:e},\"clipped\":{clipped}}}",
            json_string(path)
        );
    }
    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"lip-input\",\
         \"wav\":{},\"sample_rate_hz\":{},\"samples\":{samples},\"clipped\":{},\"peak_pa\":{:e},\
         \"rms_pa\":{:e},\"wav_blake3\":\"{wav_hash}\",\"input_blake3\":\"{}\"}}",
        json_string(&output.display().to_string()),
        info.sample_rate_hz,
        rendered.clipped,
        rendered.peak_pa,
        rendered.rms_pa,
        info.input_hash.to_hex(),
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn options_refuse_physical_overrides_and_missing_paths() {
        let args = |v: &[&str]| v.iter().map(|s| (*s).to_string()).collect::<Vec<_>>();
        assert!(options(&args(&["in.performance", "out.wav"])).is_ok());
        assert_eq!(
            options(&args(&["a", "b", "--mouthpiece-wav", "m.wav"]))
                .unwrap()
                .3,
            Some("m.wav")
        );
        assert!(options(&args(&["a", "b", "--mouthpiece-wav"])).is_err());
        assert_eq!(
            options(&args(&["a", "b", "--impedance", "z.tsv"]))
                .unwrap()
                .4,
            Some("z.tsv")
        );
        assert!(options(&args(&["a", "b", "--impedance"])).is_err());
        assert!(
            options(&args(&[
                "a",
                "b",
                "--mouthpiece-wav",
                "m.wav",
                "--mouthpiece-wav",
                "n.wav"
            ]))
            .is_err()
        );
        assert_eq!(options(&args(&["a", "b", "--block", "64"])).unwrap().2, 64);
        assert!(options(&args(&["a", "b", "--seconds", "3"])).is_err());
        assert!(options(&args(&["a"])).is_err());
        assert!(options(&args(&["a", "b", "--block", "0"])).is_err());
        assert!(options(&args(&["a", "b", "--block", "8", "--block", "8"])).is_err());
    }
}
