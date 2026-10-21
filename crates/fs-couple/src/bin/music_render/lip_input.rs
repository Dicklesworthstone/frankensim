//! File-driven lip-valve performances: `music_render lip INPUT.performance OUT.wav`.
//!
//! The mechanics clock of a bore with a large radiating mouth is bounded by the
//! flanged-radiation ceiling (`ka <= 20` at Nyquist), so the file declares its
//! own clock and this command writes the WAV at that clock. Nothing is
//! resampled; the sidecar states the clock and the declared acoustic band.
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use fs_blake3::DomainHasher;
use fs_couple::brass_loop::performance::{LipDescription, LipInfo, MAX_LIP_PERFORMANCE_BYTES};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::stream::{MAX_PCM16_WAV_SAMPLES, Pcm16WavStream};

use super::{WAV_HASH_DOMAIN, create_outputs, json_string};

fn options(args: &[String]) -> Result<(&str, &str, usize), String> {
    let mut paths = Vec::new();
    let mut block = None;
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
            value if value.starts_with('-') => {
                return Err(format!(
                    "unsupported lip option {value:?}; geometry, lip, pressure, duration and full-scale belong to the source file"
                ));
            }
            value => paths.push(value),
        }
    }
    let [input, output] = paths.as_slice() else {
        return Err("usage: music_render lip INPUT.performance OUT.wav [--block N]".into());
    };
    Ok((*input, *output, block.unwrap_or(512)))
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

fn input_json(i: &LipInfo) -> String {
    let ladder = i
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
    let lip = &i.lip;
    let mic = &i.microphone;
    let line = &i.line;
    format!(
        "\"lip_input\":{{\"schema\":\"{}\",\"blake3\":\"{}\",\"bore\":{{\"segments\":{},\"length_m\":{},\"mouth_radius_m\":{},\"gas_temperature_k\":{},\"sound_speed_m_s\":{},\"source\":{}}},\
         \"lip\":{{\"mass_kg\":{},\"stiffness_n_m\":{},\"damping_n_s_m\":{},\"width_m\":{},\"rest_gap_m\":{},\"face_area_m2\":{},\"source\":{}}},\
         \"cup_volume_m3\":{},\"control_rate_hz\":{},\"lip_tension_track\":{},\
         \"reflectance_line\":{{\"taps\":{},\"round_trip_samples\":{},\"peak_reflectance\":{},\"rungs_tried\":{},\"tail_energy\":{},\"peak_shift_cents\":{}}},\
         \"mouth_flow_line\":{{\"taps\":{},\"rungs\":{},\"tail_tol\":{},\"ladder\":[{}]}},\
         \"microphone\":{{\"position_m\":[{},{},{}],\"radial_rings\":{},\"angular_points\":{},\"resolved_band_hz\":{}}},\
         \"model_scope\":\"one held bore; 1-DOF outward-striking lip island; plane-mode multimodal reflectance line; plane-mode mouth-flow transfer; baffled-disk Rayleigh radiation at the declared point; no valve or slide motion, no bore wall vibration, no breath heating, no nonlinear propagation; the declared band is the resolved band of the radiation quadrature and the clock's flanged-radiation ceiling\"}}",
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
        number(i.cup_volume_m3),
        i.control_rate_hz,
        i.has_lip_tension_track,
        line.n_fft,
        number(line.round_trip_samples),
        number(line.peak_reflectance),
        line.rungs_tried,
        number(line.tail_energy),
        number(line.peak_shift_cents),
        i.mouth_flow_taps,
        i.mouth_flow_rungs,
        number(i.mouth_flow_tail_tol),
        ladder,
        number(mic.position_m[0]),
        number(mic.position_m[1]),
        number(mic.position_m[2]),
        mic.radial_rings,
        mic.angular_points,
        number(mic.maximum_frequency_hz),
    )
}

pub(super) fn run(args: &[String]) -> Result<(), String> {
    let (input, output, block) = options(args)?;
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
    let description = LipDescription::from_bytes(&bytes)?;
    let mut performance = description.build()?;
    let info = performance.info().clone();
    let samples =
        usize::try_from(info.samples).map_err(|_| "sample count exceeds this host".to_string())?;
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = render_at_rate(
        &mut performance,
        &mut audio,
        info.sample_rate_hz,
        samples,
        block,
        info.full_scale_pa,
    )?;
    let wav_hash = rendered.hash.to_hex();
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"lip-input\",\
         \"sample_rate_hz\":{},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",{},\
         \"observation\":{{\"mechanics_sample_rate_hz\":{},\"output_sample_rate_hz\":{},\"ratio\":1,\
         \"filter\":\"none: the baffled receiver differentiates and delays the mouth flow of every mechanical sample\",\
         \"resampled\":false,\"declared_band_hz\":{}}}}}",
        info.sample_rate_hz,
        info.full_scale_pa,
        rendered.clipped,
        rendered.peak_pa,
        rendered.rms_pa,
        input_json(&info),
        info.sample_rate_hz,
        info.sample_rate_hz,
        number(info.microphone.maximum_frequency_hz),
    );
    writeln!(metadata, "{provenance}")
        .and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
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
        assert_eq!(options(&args(&["a", "b", "--block", "64"])).unwrap().2, 64);
        assert!(options(&args(&["a", "b", "--seconds", "3"])).is_err());
        assert!(options(&args(&["a"])).is_err());
        assert!(options(&args(&["a", "b", "--block", "0"])).is_err());
        assert!(options(&args(&["a", "b", "--block", "8", "--block", "8"])).is_err());
    }
}
