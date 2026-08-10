//! Offline, sequential denoising of Euler-disc FinalDiagnostic EXR frames.
//!
//! This command is deliberately a bounded display-derivative producer. It
//! admits only the 30 float planes of the `FinalDiagnostic` AOV profile at no
//! more than 3840x2160, uses its material palette indices as exact denoising
//! edge labels, keeps only the immediately preceding biased denoise frame, and
//! never relabels its PNGs as raw render estimates.

use std::{
    collections::BTreeMap,
    ffi::OsString,
    fmt::Write as FmtWrite,
    fs::{self, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
};

use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_euler_disc_e2e::cinematic_fixture::critique_color_config;
use fs_img::{
    Channel, CinematicColorConfig, CinematicColorLimits, DecodedExr, ExrAttribute,
    ExrInspectLimits, MAX_TEMPORAL_DENOISE_SPATIAL_ITERATIONS, PixelType, PngColor,
    TEMPORAL_DENOISE_PIPELINE_VERSION, TemporalDenoiseConfig, TemporalDenoiseInput,
    TemporalDenoiseLimits, TemporalDenoisedFrame, TemporalFrameBoundary, inspect_exr, read_exr,
    temporal_denoise_rgb, transform_cinematic_preview, write_png16,
};
use fs_render::{
    aov::{
        CINEMATIC_AOV_CHANNEL_SEMANTICS, CINEMATIC_AOV_INVALID_SEMANTICS,
        CINEMATIC_AOV_PALETTE_ZERO_SEMANTICS, CINEMATIC_AOV_SEMANTICS_VERSION, CinematicAovProfile,
    },
    tracer::{
        ADAPTIVE_SAMPLING_SEMANTICS_VERSION, INDEPENDENT_PILOT_ALLOCATION_SEMANTICS_VERSION,
        MATERIAL_CONTENT_IDENTITY_DOMAIN,
    },
};

const FINAL_DIAGNOSTIC_CHANNELS: &[(&str, PixelType)] =
    CinematicAovProfile::FinalDiagnostic.exr_channel_layout();
const MAX_4K_PIXELS: u64 = 3_840 * 2_160;
const MAX_FINAL_DIAGNOSTIC_DECODED_BYTES: u64 =
    MAX_4K_PIXELS * FINAL_DIAGNOSTIC_CHANNELS.len() as u64 * 4;
const MAX_FINAL_DIAGNOSTIC_ENCODED_BYTES: u64 = 1_024 * 1_024 * 1_024;
const MAX_EXR_HEADER_BYTES: u64 = 1024 * 1024;
const MAX_EXR_METADATA_BYTES: u64 = 1024 * 1024;
const MAX_EXACT_F32_INTEGER: u64 = 1 << 24;
const FINAL_DIAGNOSTIC_PROFILE: &str = CinematicAovProfile::FinalDiagnostic.code();

#[derive(Debug, Clone, PartialEq)]
struct Cli {
    input: PathBuf,
    output: PathBuf,
    frame_start: u64,
    frame_count: u64,
    initial_cut: bool,
    allow_uniform_spp_transition_at: Option<u64>,
    denoise_spatial_passes: u8,
    denoise_spatial_sigma: f32,
}

#[derive(Debug, Clone, PartialEq)]
enum Command {
    Denoise(Cli),
    InspectSamples(PathBuf),
}

#[derive(Debug)]
struct FinalDiagnosticFrame {
    width: u32,
    height: u32,
    red: Vec<f32>,
    green: Vec<f32>,
    blue: Vec<f32>,
    motion_prev_x: Vec<f32>,
    motion_prev_y: Vec<f32>,
    axial_depth_m: Vec<f32>,
    normal_x: Vec<f32>,
    normal_y: Vec<f32>,
    normal_z: Vec<f32>,
    primary_coverage: Vec<f32>,
    variance_luminance: Vec<f32>,
    sample_counts: Vec<u32>,
    object_ids: Vec<u64>,
    material_ids: Vec<u64>,
    sequence_identity: SequenceIdentity,
    timing: FrameTiming,
}

/// Header values which must name one coherent raw-render sequence before the
/// history-dependent denoiser is allowed to join frames.
#[derive(Debug, Clone, PartialEq, Eq)]
struct SequenceIdentity {
    source_trajectory: String,
    scene_hash: String,
    composition: String,
    aov_profile: String,
    sample_mode: String,
    sample_ceiling: u32,
    adaptive_policy: Option<String>,
    independent_pilot_policy: Option<String>,
    object_palette: String,
    material_palette: String,
    object_palette_entries: u64,
    material_palette_entries: u64,
    shot_id: u64,
    cut_side: String,
    shutter: String,
    sampler: String,
    strategy: String,
    max_depth: u64,
    render_versions: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct UniformSppTransition {
    frame: u64,
    from_spp: u32,
    to_spp: u32,
    from_composition: String,
    to_composition: String,
    from_shutter: String,
    to_shutter: String,
}

fn admit_uniform_spp_transition(
    frame: u64,
    from: &SequenceIdentity,
    to: &SequenceIdentity,
) -> Result<UniformSppTransition, String> {
    if from.sample_mode != "uniform"
        || to.sample_mode != "uniform"
        || from.adaptive_policy.is_some()
        || to.adaptive_policy.is_some()
        || from.independent_pilot_policy.is_some()
        || to.independent_pilot_policy.is_some()
    {
        return Err("an authorized SPP transition must remain uniform on both sides".to_owned());
    }
    if from.sample_ceiling == to.sample_ceiling {
        return Err(format!(
            "authorized SPP transition at absolute frame {frame} did not change uniform SPP"
        ));
    }
    let from_suffix = format!(";strata={}", from.sample_ceiling);
    let to_suffix = format!(";strata={}", to.sample_ceiling);
    let from_shutter_family = from.shutter.strip_suffix(&from_suffix).ok_or_else(|| {
        format!(
            "authorized SPP transition at absolute frame {frame} has a source shutter that is not bound to its SPP strata"
        )
    })?;
    let to_shutter_family = to.shutter.strip_suffix(&to_suffix).ok_or_else(|| {
        format!(
            "authorized SPP transition at absolute frame {frame} has a destination shutter that is not bound to its SPP strata"
        )
    })?;
    if from_shutter_family != to_shutter_family {
        return Err(format!(
            "authorized SPP transition at absolute frame {frame} changed shutter semantics beyond the SPP stratum count"
        ));
    }

    let mut normalized_from = from.clone();
    normalized_from.sample_ceiling = 0;
    normalized_from.composition.clear();
    normalized_from.shutter = from_shutter_family.to_owned();
    let mut normalized_to = to.clone();
    normalized_to.sample_ceiling = 0;
    normalized_to.composition.clear();
    normalized_to.shutter = to_shutter_family.to_owned();
    if normalized_from != normalized_to {
        return Err(format!(
            "authorized SPP transition at absolute frame {frame} changed provenance beyond uniform SPP, its composition binding, and the matching shutter stratum count"
        ));
    }

    Ok(UniformSppTransition {
        frame,
        from_spp: from.sample_ceiling,
        to_spp: to.sample_ceiling,
        from_composition: from.composition.clone(),
        to_composition: to.composition.clone(),
        from_shutter: from.shutter.clone(),
        to_shutter: to.shutter.clone(),
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct FrameTiming {
    frame_time_bits: u64,
    previous_time_bits: u64,
    next_time_bits: u64,
    shutter_open_bits: u64,
    shutter_close_bits: u64,
}

impl FinalDiagnosticFrame {
    fn temporal_input(&self, frame_index: u64) -> TemporalDenoiseInput<'_> {
        TemporalDenoiseInput {
            frame_index,
            samples_per_pixel: self.sequence_identity.sample_ceiling,
            sample_counts_per_pixel: Some(&self.sample_counts),
            width: self.width as usize,
            height: self.height as usize,
            red: &self.red,
            green: &self.green,
            blue: &self.blue,
            motion_prev_x: &self.motion_prev_x,
            motion_prev_y: &self.motion_prev_y,
            axial_depth_m: &self.axial_depth_m,
            normal_x: &self.normal_x,
            normal_y: &self.normal_y,
            normal_z: &self.normal_z,
            primary_coverage: &self.primary_coverage,
            variance_luminance: &self.variance_luminance,
            object_ids: Some(&self.object_ids),
            material_ids: Some(&self.material_ids),
        }
    }
}

fn main() {
    if let Err(error) = run() {
        eprintln!("status=error message={error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    match parse_cli(std::env::args().skip(1))? {
        Command::Denoise(cli) => run_cli(&cli, |message| eprintln!("status=progress {message}")),
        Command::InspectSamples(path) => {
            print!("{}", inspect_sample_counts(&path)?);
            Ok(())
        }
    }
}

fn run_cli(cli: &Cli, mut progress: impl FnMut(&str)) -> Result<(), String> {
    if !cli.input.is_dir() {
        return Err(format!("input is not a directory: {}", cli.input.display()));
    }
    if cli.frame_start != 0 && !cli.initial_cut {
        return Err(format!(
            "nonzero frame-start {} requires --initial-cut so missing temporal history is explicit",
            cli.frame_start
        ));
    }
    let mut denoise_config = TemporalDenoiseConfig::default();
    denoise_config.spatial_iterations = cli.denoise_spatial_passes;
    denoise_config.spatial_sigma_rgb = cli.denoise_spatial_sigma;
    denoise_config
        .identity()
        .map_err(|error| format!("invalid denoiser configuration: {error}"))?;
    let frame_end = cli
        .frame_start
        .checked_add(cli.frame_count)
        .ok_or_else(|| "frame-start plus frame-count overflows u64".to_owned())?;
    if cli.frame_count == 0 {
        return Err("frame-count must be positive".to_owned());
    }
    require_absent(&cli.output, "final output")?;
    let staging_output = staging_output_path(&cli.output)?;
    require_absent(&staging_output, "staging output")?;
    let output_parent = cli
        .output
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_parent).map_err(|error| {
        format!(
            "create output parent directory {}: {error}",
            output_parent.display()
        )
    })?;
    fs::create_dir(&staging_output).map_err(|error| {
        format!(
            "create staging output directory {}: {error}",
            staging_output.display()
        )
    })?;

    progress(&format!(
        "stage=denoise begin input={} output={} staging_output={} frame_start={} frame_count={} initial_boundary=cut spatial_passes={} spatial_sigma={}",
        cli.input.display(),
        cli.output.display(),
        staging_output.display(),
        cli.frame_start,
        cli.frame_count,
        denoise_config.spatial_iterations,
        denoise_config.spatial_sigma_rgb,
    ));
    let color = critique_color_config();
    let mut history: Option<TemporalDenoisedFrame> = None;
    let mut expected_dimensions: Option<(u32, u32)> = None;
    let mut first_sequence: Option<SequenceIdentity> = None;
    let mut active_sequence: Option<SequenceIdentity> = None;
    let mut admitted_spp_transition: Option<UniformSppTransition> = None;
    let mut previous_timing: Option<FrameTiming> = None;
    let mut raw_sequence =
        DomainHasher::new("org.frankensim.euler-critique.offline-denoise.raw-input-sequence.v1");
    let mut preview_sequence = DomainHasher::new(
        "org.frankensim.euler-critique.offline-denoise.preview-output-sequence.v1",
    );

    let denoise_result = (|| -> Result<(), String> {
        for (ordinal, frame) in (cli.frame_start..frame_end).enumerate() {
            let input_path = raw_path(&cli.input, frame);
            progress(&format!(
                "stage=denoise frame={}/{} absolute_frame={} action=read path={}",
                ordinal + 1,
                cli.frame_count,
                frame,
                input_path.display(),
            ));
            let (raw, encoded_identity) = read_final_diagnostic(&input_path, frame)?;
            raw_sequence.update(encoded_identity.as_bytes());
            let dimensions = (raw.width, raw.height);
            if let Some(expected) = expected_dimensions {
                if dimensions != expected {
                    return Err(format!(
                        "frame continuity violation at absolute frame {frame}: dimensions {}x{} differ from {}x{}",
                        dimensions.0, dimensions.1, expected.0, expected.1
                    ));
                }
            } else {
                expected_dimensions = Some(dimensions);
            }
            if let Some(expected) = &active_sequence {
                if cli.allow_uniform_spp_transition_at == Some(frame) {
                    let transition =
                        admit_uniform_spp_transition(frame, expected, &raw.sequence_identity)?;
                    active_sequence = Some(raw.sequence_identity.clone());
                    admitted_spp_transition = Some(transition);
                } else if raw.sequence_identity != *expected {
                    return Err(format!(
                        "frame continuity violation at absolute frame {frame}: EXR provenance identity differs from the active requested segment"
                    ));
                }
            } else {
                first_sequence = Some(raw.sequence_identity.clone());
                active_sequence = Some(raw.sequence_identity.clone());
            }
            let boundary = denoise_boundary(ordinal);
            validate_timing_continuity(previous_timing, raw.timing, frame, boundary)?;
            previous_timing = Some(raw.timing);
            let denoised = temporal_denoise_rgb(
                raw.temporal_input(frame),
                history.as_ref(),
                boundary,
                denoise_config,
                TemporalDenoiseLimits::reference_4k(),
            )
            .map_err(|error| format!("temporal denoise frame {frame}: {error}"))?;
            // The denoised history owns the only guides needed by the next frame;
            // release this frame's raw AOV planes before allocating display data.
            drop(raw);
            let [red, green, blue] = denoised.linear_rgb();
            let preview = transform_cinematic_preview(
                dimensions.0,
                dimensions.1,
                [red, green, blue],
                color,
                CinematicColorLimits::reference_4k(),
            )
            .map_err(|error| format!("display transform frame {frame}: {error}"))?;
            let samples = preview.samples().as_u16().ok_or_else(|| {
                "reference 16-bit colour configuration returned non-16-bit preview".to_owned()
            })?;
            let png = write_png16(dimensions.0, dimensions.1, PngColor::Rgb, samples)
                .map_err(|error| format!("PNG16 encode frame {frame}: {error}"))?;
            preview_sequence.update(hash_domain("frame", &png).as_bytes());
            let output_path = preview_path(&staging_output, frame);
            write_new(&output_path, &png)?;
            progress(&format!(
                "stage=denoise frame={}/{} absolute_frame={} boundary={} action=wrote path={} bytes={}",
                ordinal + 1,
                cli.frame_count,
                frame,
                match boundary {
                    TemporalFrameBoundary::Cut => "cut",
                    TemporalFrameBoundary::Continuous => "continuous",
                },
                output_path.display(),
                png.len(),
            ));
            history = Some(denoised);
        }
        Ok(())
    })();
    if let Err(error) = denoise_result {
        return Err(format!(
            "{error}; incomplete staged output was preserved at {}",
            staging_output.display()
        ));
    }
    let sequence_identity = first_sequence.as_ref().ok_or_else(|| {
        format!(
            "denoiser produced no sequence identity; incomplete staged output was preserved at {}",
            staging_output.display()
        )
    })?;
    if cli.allow_uniform_spp_transition_at.is_some() && admitted_spp_transition.is_none() {
        return Err(format!(
            "the authorized uniform-SPP transition was not observed; incomplete staged output was preserved at {}",
            staging_output.display()
        ));
    }
    let dimensions = expected_dimensions.ok_or_else(|| {
        format!(
            "denoiser produced no raster dimensions; incomplete staged output was preserved at {}",
            staging_output.display()
        )
    })?;
    let manifest = offline_preview_manifest(
        cli,
        sequence_identity,
        dimensions,
        raw_sequence.finalize(),
        preview_sequence.finalize(),
        admitted_spp_transition.as_ref(),
        color,
        denoise_config,
    )
    .map_err(|error| {
        format!(
            "{error}; incomplete staged output was preserved at {}",
            staging_output.display()
        )
    })?;
    write_new(
        &staging_output.join("denoise-manifest.json"),
        manifest.as_bytes(),
    )
    .map_err(|error| {
        format!(
            "{error}; incomplete staged output was preserved at {}",
            staging_output.display()
        )
    })?;
    // This CLI is a single-writer producer. Under that declared operating
    // contract, the complete sequence becomes visible with one same-filesystem
    // directory rename; failures retain only the explicitly incomplete staging
    // path. The second check diagnoses accidental ordinary reuse, but is not a
    // cross-process no-replace primitive.
    require_absent(&cli.output, "final output")?;
    fs::rename(&staging_output, &cli.output).map_err(|error| {
        format!(
            "publish complete preview sequence {} as {}: {error}; complete staged output remains at {}",
            staging_output.display(),
            cli.output.display(),
            staging_output.display()
        )
    })?;
    progress(&format!(
        "stage=denoise complete output={} frame_start={} frame_count={}",
        cli.output.display(),
        cli.frame_start,
        cli.frame_count
    ));
    Ok(())
}

const fn denoise_boundary(ordinal: usize) -> TemporalFrameBoundary {
    if ordinal == 0 {
        TemporalFrameBoundary::Cut
    } else {
        TemporalFrameBoundary::Continuous
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SampleCountSummary {
    histogram: BTreeMap<u32, u64>,
    pixels: u64,
    total_samples: u64,
    minimum: u32,
    maximum: u32,
    at_ceiling_spp_pixels: u64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct EstimatorVarianceSummary {
    sample_variance_total: f64,
    estimator_variance_total: f64,
    estimator_variance_mean: f64,
    estimator_standard_error_rms: f64,
    maximum_estimator_variance: f64,
    maximum_pixel_index: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PixelCrop {
    x_min: u32,
    x_max: u32,
    y_min: u32,
    y_max: u32,
}

fn inspect_sample_counts(path: &Path) -> Result<String, String> {
    let (decoded, encoded_identity) = read_final_diagnostic_exr(path)?;
    let frame_index = exr_attribute_u64(&decoded, "frankensim.frame.index", path)?;
    let frame = decode_final_diagnostic(decoded, path, frame_index)?;
    sample_inspection_report(&frame, frame_index, encoded_identity)
}

/// Emit stable, whitespace-delimited records. Integer numerator/denominator
/// pairs are the exact mean and fraction evidence; decimal fields are only
/// convenient renderings of those ratios. Quantiles use the one-based nearest-
/// rank definition over exact integer sample counts. A count equal to the SPP
/// ceiling does not imply nonconvergence: the samples plane cannot recover an
/// adaptive stopping decision or the discarded pilot observations which fixed
/// an independent production allocation.
fn sample_inspection_report(
    frame: &FinalDiagnosticFrame,
    frame_index: u64,
    encoded_identity: ContentHash,
) -> Result<String, String> {
    let expected_pixels = usize::try_from(u64::from(frame.width) * u64::from(frame.height))
        .map_err(|_| "sample inspection raster size does not fit usize".to_owned())?;
    if frame.sample_counts.len() != expected_pixels {
        return Err(format!(
            "sample inspection plane has {} pixels; expected {expected_pixels}",
            frame.sample_counts.len()
        ));
    }
    let ceiling = frame.sequence_identity.sample_ceiling;
    let full = summarize_sample_counts(frame.sample_counts.iter().copied(), ceiling)?;
    let variance = summarize_estimator_variance(&frame.variance_luminance, &frame.sample_counts)?;
    let adaptive_policy = frame
        .sequence_identity
        .adaptive_policy
        .as_deref()
        .unwrap_or("none");
    let independent_pilot_policy = frame
        .sequence_identity
        .independent_pilot_policy
        .as_deref()
        .unwrap_or("none");
    let mut report = String::new();
    writeln!(
        report,
        "record=metadata schema=frankensim-euler-sample-inspection-v2 source_identity={} frame={} width={} height={} sample_mode={} sample_ceiling={} adaptive_policy={} independent_pilot_policy={} quantile_method=nearest-rank allocation_decision=unavailable-from-samples-plane",
        encoded_identity.to_hex(),
        frame_index,
        frame.width,
        frame.height,
        frame.sequence_identity.sample_mode,
        ceiling,
        adaptive_policy,
        independent_pilot_policy,
    )
    .expect("writing to a String cannot fail");
    append_sample_summary(&mut report, "record=summary scope=full", &full);
    writeln!(
        report,
        "record=uncertainty scope=full channel=variance.Y meaning=unbiased-raw-CIE-Y-sample-variance sample_variance_total={:.17e} estimator_variance_total={:.17e} estimator_variance_mean={:.17e} estimator_standard_error_rms={:.17e} maximum_estimator_variance={:.17e} maximum_pixel_index={}",
        variance.sample_variance_total,
        variance.estimator_variance_total,
        variance.estimator_variance_mean,
        variance.estimator_standard_error_rms,
        variance.maximum_estimator_variance,
        variance.maximum_pixel_index,
    )
    .expect("writing to a String cannot fail");
    for (&samples, &pixels) in &full.histogram {
        writeln!(
            report,
            "record=histogram scope=full spp={samples} pixels={pixels} fraction_numerator={pixels} fraction_denominator={} fraction={:.9}",
            full.pixels,
            pixels as f64 / full.pixels as f64,
        )
        .expect("writing to a String cannot fail");
    }
    // These normalized crops are the canonical 320x180 calibration regions
    // used by the Euler-disc look-development matrix. Scaling their half-open
    // bounds makes the same named evidence available at 960p and 4K.
    for (name, crop) in [
        (
            "disc",
            reference_crop(frame.width, frame.height, 108, 212, 54, 118),
        ),
        (
            "front_glass",
            reference_crop(frame.width, frame.height, 10, 230, 105, 173),
        ),
        (
            "right_glass",
            reference_crop(frame.width, frame.height, 190, 290, 65, 165),
        ),
        (
            "background",
            reference_crop(frame.width, frame.height, 0, 100, 0, 40),
        ),
    ] {
        let summary = summarize_crop(&frame.sample_counts, frame.width, crop, ceiling)?;
        let prefix = format!(
            "record=crop name={name} x_min={} x_max={} y_min={} y_max={}",
            crop.x_min, crop.x_max, crop.y_min, crop.y_max
        );
        append_sample_summary(&mut report, &prefix, &summary);
    }
    Ok(report)
}

fn summarize_estimator_variance(
    sample_variance: &[f32],
    sample_counts: &[u32],
) -> Result<EstimatorVarianceSummary, String> {
    if sample_variance.is_empty() || sample_variance.len() != sample_counts.len() {
        return Err(
            "variance and sample-count planes must have one common nonempty raster".to_owned(),
        );
    }
    let mut sample_variance_total = 0.0_f64;
    let mut estimator_variance_total = 0.0_f64;
    let mut maximum_estimator_variance = -1.0_f64;
    let mut maximum_pixel_index = 0_usize;
    for (index, (&variance, &samples)) in sample_variance.iter().zip(sample_counts).enumerate() {
        if !variance.is_finite() || variance < 0.0 || samples == 0 {
            return Err(format!(
                "variance summary encountered invalid pixel {index}: variance={variance}, samples={samples}"
            ));
        }
        let variance = f64::from(variance);
        let estimator_variance = variance / f64::from(samples);
        sample_variance_total += variance;
        estimator_variance_total += estimator_variance;
        if estimator_variance > maximum_estimator_variance {
            maximum_estimator_variance = estimator_variance;
            maximum_pixel_index = index;
        }
    }
    if !sample_variance_total.is_finite() || !estimator_variance_total.is_finite() {
        return Err("variance summary arithmetic overflowed".to_owned());
    }
    let estimator_variance_mean = estimator_variance_total / sample_variance.len() as f64;
    Ok(EstimatorVarianceSummary {
        sample_variance_total,
        estimator_variance_total,
        estimator_variance_mean,
        estimator_standard_error_rms: estimator_variance_mean.sqrt(),
        maximum_estimator_variance,
        maximum_pixel_index,
    })
}

fn summarize_sample_counts(
    counts: impl IntoIterator<Item = u32>,
    ceiling: u32,
) -> Result<SampleCountSummary, String> {
    let mut histogram = BTreeMap::<u32, u64>::new();
    let mut pixels = 0_u64;
    let mut total_samples = 0_u64;
    let mut minimum = u32::MAX;
    let mut maximum = 0_u32;
    let mut at_ceiling_spp_pixels = 0_u64;
    for samples in counts {
        if samples == 0 || samples > ceiling {
            return Err(format!(
                "sample inspection encountered count {samples} outside 1..={ceiling}"
            ));
        }
        pixels = pixels
            .checked_add(1)
            .ok_or_else(|| "sample inspection pixel count overflowed u64".to_owned())?;
        total_samples = total_samples
            .checked_add(u64::from(samples))
            .ok_or_else(|| "sample inspection total sample count overflowed u64".to_owned())?;
        let bin = histogram.entry(samples).or_insert(0);
        *bin = bin
            .checked_add(1)
            .ok_or_else(|| "sample inspection histogram count overflowed u64".to_owned())?;
        minimum = minimum.min(samples);
        maximum = maximum.max(samples);
        at_ceiling_spp_pixels += u64::from(samples == ceiling);
    }
    if pixels == 0 {
        return Err("sample inspection region is empty".to_owned());
    }
    Ok(SampleCountSummary {
        histogram,
        pixels,
        total_samples,
        minimum,
        maximum,
        at_ceiling_spp_pixels,
    })
}

fn summarize_crop(
    counts: &[u32],
    raster_width: u32,
    crop: PixelCrop,
    ceiling: u32,
) -> Result<SampleCountSummary, String> {
    let width = usize::try_from(raster_width)
        .map_err(|_| "sample inspection raster width does not fit usize".to_owned())?;
    let x_min = crop.x_min as usize;
    let x_max = crop.x_max as usize;
    let y_min = crop.y_min as usize;
    let y_max = crop.y_max as usize;
    let crop_counts = (y_min..y_max).flat_map(|y| {
        let row = y * width;
        counts[row + x_min..row + x_max].iter().copied()
    });
    summarize_sample_counts(crop_counts, ceiling)
}

fn reference_crop(
    width: u32,
    height: u32,
    reference_x_min: u32,
    reference_x_max: u32,
    reference_y_min: u32,
    reference_y_max: u32,
) -> PixelCrop {
    const REFERENCE_WIDTH: u64 = 320;
    const REFERENCE_HEIGHT: u64 = 180;
    debug_assert!(reference_x_min < reference_x_max && reference_x_max <= 320);
    debug_assert!(reference_y_min < reference_y_max && reference_y_max <= 180);
    let x_min = (u64::from(width) * u64::from(reference_x_min) / REFERENCE_WIDTH) as u32;
    let x_max = (u64::from(width) * u64::from(reference_x_max)).div_ceil(REFERENCE_WIDTH) as u32;
    let y_min = (u64::from(height) * u64::from(reference_y_min) / REFERENCE_HEIGHT) as u32;
    let y_max = (u64::from(height) * u64::from(reference_y_max)).div_ceil(REFERENCE_HEIGHT) as u32;
    PixelCrop {
        x_min,
        x_max: x_max.max(x_min + 1).min(width),
        y_min,
        y_max: y_max.max(y_min + 1).min(height),
    }
}

fn append_sample_summary(report: &mut String, prefix: &str, summary: &SampleCountSummary) {
    writeln!(
        report,
        "{prefix} pixels={} total_samples={} mean_numerator={} mean_denominator={} mean_spp={:.9} min_spp={} p10_spp={} p25_spp={} p50_spp={} p75_spp={} p90_spp={} p95_spp={} p99_spp={} max_spp={} at_ceiling_spp_pixels={} at_ceiling_spp_fraction_numerator={} at_ceiling_spp_fraction_denominator={} at_ceiling_spp_fraction={:.9}",
        summary.pixels,
        summary.total_samples,
        summary.total_samples,
        summary.pixels,
        summary.total_samples as f64 / summary.pixels as f64,
        summary.minimum,
        nearest_rank_quantile(summary, 10, 100),
        nearest_rank_quantile(summary, 25, 100),
        nearest_rank_quantile(summary, 50, 100),
        nearest_rank_quantile(summary, 75, 100),
        nearest_rank_quantile(summary, 90, 100),
        nearest_rank_quantile(summary, 95, 100),
        nearest_rank_quantile(summary, 99, 100),
        summary.maximum,
        summary.at_ceiling_spp_pixels,
        summary.at_ceiling_spp_pixels,
        summary.pixels,
        summary.at_ceiling_spp_pixels as f64 / summary.pixels as f64,
    )
    .expect("writing to a String cannot fail");
}

fn nearest_rank_quantile(summary: &SampleCountSummary, numerator: u64, denominator: u64) -> u32 {
    debug_assert!(numerator > 0 && numerator <= denominator);
    let rank = (summary.pixels * numerator).div_ceil(denominator);
    let mut cumulative = 0_u64;
    for (&samples, &pixels) in &summary.histogram {
        cumulative += pixels;
        if cumulative >= rank {
            return samples;
        }
    }
    summary.maximum
}

fn read_final_diagnostic(
    path: &Path,
    expected_frame: u64,
) -> Result<(FinalDiagnosticFrame, ContentHash), String> {
    let (decoded, encoded_identity) = read_final_diagnostic_exr(path)?;
    decode_final_diagnostic(decoded, path, expected_frame).map(|frame| (frame, encoded_identity))
}

fn read_final_diagnostic_exr(path: &Path) -> Result<(DecodedExr, ContentHash), String> {
    let length = fs::metadata(path)
        .map_err(|error| format!("inspect input {}: {error}", path.display()))?
        .len();
    if length > MAX_FINAL_DIAGNOSTIC_ENCODED_BYTES {
        return Err(format!(
            "FinalDiagnostic EXR {} is {length} bytes; maximum encoded input is {MAX_FINAL_DIAGNOSTIC_ENCODED_BYTES} bytes",
            path.display()
        ));
    }
    let bytes =
        fs::read(path).map_err(|error| format!("read input {}: {error}", path.display()))?;
    let encoded_identity = hash_domain(
        "org.frankensim.euler-critique.final-diagnostic-frame.v1",
        &bytes,
    );
    let inspection = inspect_exr(
        &bytes,
        ExrInspectLimits {
            max_input_bytes: MAX_FINAL_DIAGNOSTIC_ENCODED_BYTES,
            max_header_bytes: MAX_EXR_HEADER_BYTES,
            max_decoded_bytes: MAX_FINAL_DIAGNOSTIC_DECODED_BYTES,
            max_metadata_bytes: MAX_EXR_METADATA_BYTES,
        },
    )
    .map_err(|error| format!("inspect FinalDiagnostic EXR {}: {error}", path.display()))?;
    validate_dimensions(inspection.width, inspection.height, path)?;
    let decoded = read_exr(&bytes)
        .map_err(|error| format!("decode FinalDiagnostic EXR {}: {error}", path.display()))?;
    drop(bytes);
    Ok((decoded, encoded_identity))
}

fn decode_final_diagnostic(
    exr: DecodedExr,
    path: &Path,
    expected_frame: u64,
) -> Result<FinalDiagnosticFrame, String> {
    validate_dimensions(exr.width, exr.height, path)?;
    let sequence_identity = validate_final_diagnostic_attributes(&exr, path, expected_frame)?;
    let independent_pilot_plan = (sequence_identity.sample_mode == "independent-pilot-fixed-v1")
        .then(|| {
            validate_independent_pilot_plan_attribute(&exr, sequence_identity.sample_ceiling, path)
        })
        .transpose()?;
    let timing = validate_frame_timing(&exr, path)?;
    let pixels = usize::try_from(u64::from(exr.width) * u64::from(exr.height)).map_err(|_| {
        format!(
            "FinalDiagnostic EXR {} pixel count does not fit usize",
            path.display()
        )
    })?;
    let mut planes = BTreeMap::new();
    for Channel { name, ty, data } in exr.channels {
        if ty != PixelType::Float {
            return Err(format!(
                "FinalDiagnostic EXR {} channel {name} is not FLOAT",
                path.display()
            ));
        }
        if data.len() != pixels {
            return Err(format!(
                "FinalDiagnostic EXR {} channel {name} has {} samples; expected {pixels}",
                path.display(),
                data.len()
            ));
        }
        if planes.insert(name.clone(), data).is_some() {
            return Err(format!(
                "FinalDiagnostic EXR {} contains duplicate channel {name}",
                path.display()
            ));
        }
    }
    for &(name, _) in FINAL_DIAGNOSTIC_CHANNELS {
        if !planes.contains_key(name) {
            return Err(format!(
                "FinalDiagnostic EXR {} is missing required channel {name}",
                path.display()
            ));
        }
    }
    if let Some((unexpected, _)) = planes.iter().find(|(name, _)| {
        !FINAL_DIAGNOSTIC_CHANNELS
            .iter()
            .any(|(expected, _)| *expected == name.as_str())
    }) {
        return Err(format!(
            "FinalDiagnostic EXR {} contains unexpected channel {unexpected}",
            path.display()
        ));
    }
    let primary_coverage = take_plane(&mut planes, "primary.coverage", path)?;
    let sample_counts = exact_sample_counts(
        take_plane(&mut planes, "samples", path)?,
        &sequence_identity,
        independent_pilot_plan.as_ref(),
        path,
    )?;
    let object_ids = exact_palette_indices(
        "id.object",
        take_plane(&mut planes, "id.object", path)?,
        &primary_coverage,
        sequence_identity.object_palette_entries,
        path,
    )?;
    let material_ids = exact_palette_indices(
        "id.material",
        take_plane(&mut planes, "id.material", path)?,
        &primary_coverage,
        sequence_identity.material_palette_entries,
        path,
    )?;
    Ok(FinalDiagnosticFrame {
        width: exr.width,
        height: exr.height,
        blue: take_plane(&mut planes, "B", path)?,
        green: take_plane(&mut planes, "G", path)?,
        red: take_plane(&mut planes, "R", path)?,
        axial_depth_m: take_plane(&mut planes, "depth.Z", path)?,
        motion_prev_x: take_plane(&mut planes, "motion.prev.X", path)?,
        motion_prev_y: take_plane(&mut planes, "motion.prev.Y", path)?,
        normal_x: take_plane(&mut planes, "normal.X", path)?,
        normal_y: take_plane(&mut planes, "normal.Y", path)?,
        normal_z: take_plane(&mut planes, "normal.Z", path)?,
        primary_coverage,
        variance_luminance: take_plane(&mut planes, "variance.Y", path)?,
        sample_counts,
        object_ids,
        material_ids,
        sequence_identity,
        timing,
    })
}

fn exact_sample_counts(
    values: Vec<f32>,
    sequence: &SequenceIdentity,
    independent_pilot_plan: Option<&IndependentPilotPlanHeader>,
    path: &Path,
) -> Result<Vec<u32>, String> {
    let mut counts = Vec::new();
    counts
        .try_reserve_exact(values.len())
        .map_err(|_| format!("allocate sample-count plane for {}", path.display()))?;
    for (index, value) in values.into_iter().enumerate() {
        let samples = value as u32;
        if !value.is_finite()
            || value < 1.0
            || value > sequence.sample_ceiling as f32
            || samples as f32 != value
        {
            return Err(format!(
                "FinalDiagnostic EXR {} samples plane has invalid exact count at sample {index}: value={value}, ceiling={}",
                path.display(),
                sequence.sample_ceiling,
            ));
        }
        if sequence.sample_mode == "uniform" && samples != sequence.sample_ceiling {
            return Err(format!(
                "FinalDiagnostic EXR {} samples plane disagrees with uniform header SPP at sample {index}: value={value}, header={}",
                path.display(),
                sequence.sample_ceiling,
            ));
        }
        if let Some(plan) = independent_pilot_plan
            && samples < plan.minimum_samples
        {
            return Err(format!(
                "FinalDiagnostic EXR {} samples plane is below the independent-pilot minimum at sample {index}: value={value}, minimum={}",
                path.display(),
                plan.minimum_samples,
            ));
        }
        counts.push(samples);
    }
    if let Some(plan) = independent_pilot_plan {
        let total_samples = counts
            .iter()
            .try_fold(0_u64, |total, &count| total.checked_add(u64::from(count)));
        if total_samples != Some(plan.total_samples) {
            return Err(format!(
                "FinalDiagnostic EXR {} independent-pilot samples plane total {:?} disagrees with declared total {}",
                path.display(),
                total_samples,
                plan.total_samples,
            ));
        }
    }
    Ok(counts)
}

fn take_plane(
    planes: &mut BTreeMap<String, Vec<f32>>,
    name: &str,
    path: &Path,
) -> Result<Vec<f32>, String> {
    planes.remove(name).ok_or_else(|| {
        format!(
            "FinalDiagnostic EXR {} lost validated required channel {name} during reconstruction",
            path.display()
        )
    })
}

fn exact_palette_indices(
    channel: &'static str,
    values: Vec<f32>,
    primary_coverage: &[f32],
    maximum_palette_index: u64,
    path: &Path,
) -> Result<Vec<u64>, String> {
    if values.len() != primary_coverage.len() {
        return Err(format!(
            "FinalDiagnostic EXR {} {channel} and coverage plane lengths differ",
            path.display()
        ));
    }
    let mut indices = Vec::new();
    indices.try_reserve_exact(values.len()).map_err(|_| {
        format!(
            "allocate {channel} palette indices for FinalDiagnostic EXR {}",
            path.display()
        )
    })?;
    for (index, value) in values.into_iter().enumerate() {
        if !value.is_finite()
            || value < 0.0
            || value >= MAX_EXACT_F32_INTEGER as f32
            || value.fract() != 0.0
        {
            return Err(format!(
                "FinalDiagnostic EXR {} {channel} sample {index} is not an exact nonnegative f32 integer palette index: {value}",
                path.display()
            ));
        }
        let palette_index = value as u64;
        if palette_index > maximum_palette_index {
            return Err(format!(
                "FinalDiagnostic EXR {} {channel} sample {index} references palette index {palette_index} above declared maximum {maximum_palette_index}",
                path.display()
            ));
        }
        let covered = primary_coverage[index] > 0.0;
        if !covered && palette_index != 0 {
            return Err(format!(
                "FinalDiagnostic EXR {} {channel} sample {index} disagrees with primary coverage",
                path.display()
            ));
        }
        indices.push(palette_index);
    }
    Ok(indices)
}

fn validate_final_diagnostic_attributes(
    exr: &DecodedExr,
    path: &Path,
    expected_frame: u64,
) -> Result<SequenceIdentity, String> {
    require_string_attribute(exr, "frankensim.aov.authority", "raw-estimate", path)?;
    require_string_attribute(
        exr,
        "frankensim.aov.schemaVersion",
        &CINEMATIC_AOV_SEMANTICS_VERSION.to_string(),
        path,
    )?;
    require_string_attribute(
        exr,
        "frankensim.aov.channelSemantics",
        CINEMATIC_AOV_CHANNEL_SEMANTICS,
        path,
    )?;
    require_string_attribute(
        exr,
        "frankensim.aov.invalidSemantics",
        CINEMATIC_AOV_INVALID_SEMANTICS,
        path,
    )?;
    require_string_attribute(
        exr,
        "frankensim.aov.materialDomain",
        MATERIAL_CONTENT_IDENTITY_DOMAIN,
        path,
    )?;
    require_string_attribute(
        exr,
        "frankensim.aov.paletteZero",
        CINEMATIC_AOV_PALETTE_ZERO_SEMANTICS,
        path,
    )?;
    let frame_index = exr_attribute_u64(exr, "frankensim.frame.index", path)?;
    if frame_index != expected_frame {
        return Err(format!(
            "FinalDiagnostic EXR {} declares frame index {frame_index}; expected {expected_frame}",
            path.display()
        ));
    }
    let aov_profile = exr_attribute_string(exr, "frankensim.aov.profile", path)?;
    if aov_profile != FINAL_DIAGNOSTIC_PROFILE {
        return Err(format!(
            "FinalDiagnostic EXR {} declares AOV profile {aov_profile:?}; expected {FINAL_DIAGNOSTIC_PROFILE:?}",
            path.display()
        ));
    }
    // The per-frame AOV configuration hash is mandatory integrity metadata,
    // but it is not a sequence identity: the hash deliberately commits to the
    // absolute frame index and neighbouring frame times.  Requiring equality
    // across frames would reject every valid temporal sequence at frame 1.
    let _frame_config_hash = content_hash_attribute(exr, "frankensim.aov.configHash", path)?;
    let sample_mode = exr_attribute_string(exr, "frankensim.render.sampleMode", path)?;
    if sample_mode != "uniform"
        && sample_mode != "adaptive"
        && sample_mode != "independent-pilot-fixed-v1"
    {
        return Err(format!(
            "FinalDiagnostic EXR {} declares unsupported sample mode {sample_mode:?}; expected \"uniform\", \"adaptive\", or \"independent-pilot-fixed-v1\"",
            path.display()
        ));
    }
    let sample_ceiling = u32::try_from(exr_attribute_u64(
        exr,
        "frankensim.render.sppCeiling",
        path,
    )?)
    .map_err(|_| {
        format!(
            "FinalDiagnostic EXR {} render sample ceiling exceeds u32",
            path.display()
        )
    })?;
    if sample_ceiling == 0 {
        return Err(format!(
            "FinalDiagnostic EXR {} render sample ceiling must be positive",
            path.display()
        ));
    }
    if u64::from(sample_ceiling) > MAX_EXACT_F32_INTEGER {
        return Err(format!(
            "FinalDiagnostic EXR {} render sample ceiling {sample_ceiling} exceeds exact FLOAT integer ceiling {MAX_EXACT_F32_INTEGER}",
            path.display()
        ));
    }
    let rendered_spp = exr_attribute_string(exr, "frankensim.render.spp", path)?;
    let (adaptive_policy, independent_pilot_policy) = match sample_mode.as_str() {
        "uniform" => {
            let samples_per_pixel = rendered_spp.parse::<u32>().map_err(|_| {
                format!(
                    "FinalDiagnostic EXR {} uniform render SPP is not canonical u32",
                    path.display()
                )
            })?;
            if samples_per_pixel.to_string() != rendered_spp || samples_per_pixel != sample_ceiling
            {
                return Err(format!(
                    "FinalDiagnostic EXR {} uniform render SPP {rendered_spp:?} disagrees with SPP ceiling {sample_ceiling}",
                    path.display()
                ));
            }
            (None, None)
        }
        "adaptive" => {
            if rendered_spp != "per-pixel-channel" {
                return Err(format!(
                    "FinalDiagnostic EXR {} adaptive render SPP must be \"per-pixel-channel\"; got {rendered_spp:?}",
                    path.display()
                ));
            }
            (
                Some(validate_adaptive_policy_attribute(
                    exr,
                    sample_ceiling,
                    path,
                )?),
                None,
            )
        }
        "independent-pilot-fixed-v1" => {
            if rendered_spp != "per-pixel-channel" {
                return Err(format!(
                    "FinalDiagnostic EXR {} independent-pilot render SPP must be \"per-pixel-channel\"; got {rendered_spp:?}",
                    path.display()
                ));
            }
            let plan = validate_independent_pilot_plan_attribute(exr, sample_ceiling, path)?;
            (None, Some(plan.sequence_policy))
        }
        _ => unreachable!("sample mode was validated above"),
    };
    let (object_palette, object_palette_entries) =
        canonical_object_palette(exr, "frankensim.aov.objectPalette", path)?;
    let (material_palette, material_palette_entries) =
        canonical_material_palette(exr, "frankensim.aov.materialPalette", path)?;
    let cut_side = nonempty_string_attribute(exr, "frankensim.render.cutSide", path)?;
    if cut_side != "before" && cut_side != "after" {
        return Err(format!(
            "FinalDiagnostic EXR {} has unsupported cut side {cut_side:?}",
            path.display()
        ));
    }
    Ok(SequenceIdentity {
        source_trajectory: content_hash_attribute(exr, "frankensim.source.trajectory", path)?,
        scene_hash: content_hash_attribute(exr, "frankensim.source.sceneHash", path)?,
        composition: content_hash_attribute(exr, "frankensim.source.composition", path)?,
        aov_profile,
        sample_mode,
        sample_ceiling,
        adaptive_policy,
        independent_pilot_policy,
        object_palette,
        material_palette,
        object_palette_entries,
        material_palette_entries,
        shot_id: exr_attribute_u64(exr, "frankensim.render.shotId", path)?,
        cut_side,
        shutter: nonempty_string_attribute(exr, "frankensim.render.shutter", path)?,
        sampler: nonempty_string_attribute(exr, "frankensim.render.sampler", path)?,
        strategy: nonempty_string_attribute(exr, "frankensim.render.strategy", path)?,
        max_depth: exr_attribute_u64(exr, "frankensim.render.maxDepth", path)?,
        render_versions: nonempty_string_attribute(exr, "frankensim.render.versions", path)?,
    })
}

fn validate_adaptive_policy_attribute(
    exr: &DecodedExr,
    sample_ceiling: u32,
    path: &Path,
) -> Result<String, String> {
    let policy = exr_attribute_string(exr, "frankensim.render.adaptive", path)?;
    let mut fields = policy.split(';');
    let version = fields
        .next()
        .and_then(|field| field.strip_prefix("version="));
    let minimum = fields
        .next()
        .and_then(|field| field.strip_prefix("minimum="));
    let batch = fields.next().and_then(|field| field.strip_prefix("batch="));
    let absolute = fields
        .next()
        .and_then(|field| field.strip_prefix("absolute="));
    let relative = fields
        .next()
        .and_then(|field| field.strip_prefix("relative="));
    let dark_floor = fields
        .next()
        .and_then(|field| field.strip_prefix("darkFloor="));
    if fields.next().is_some()
        || canonical_u32(version) != Some(ADAPTIVE_SAMPLING_SEMANTICS_VERSION)
    {
        return Err(invalid_adaptive_policy(path));
    }
    let minimum = canonical_u32(minimum).ok_or_else(|| invalid_adaptive_policy(path))?;
    let batch = canonical_u32(batch).ok_or_else(|| invalid_adaptive_policy(path))?;
    if minimum < 2 || minimum > sample_ceiling || batch == 0 {
        return Err(invalid_adaptive_policy(path));
    }
    for value in [absolute, relative, dark_floor] {
        let value = value.ok_or_else(|| invalid_adaptive_policy(path))?;
        let bits = canonical_f64_string_bits(value).ok_or_else(|| invalid_adaptive_policy(path))?;
        if f64::from_bits(bits) < 0.0 {
            return Err(invalid_adaptive_policy(path));
        }
    }
    Ok(policy)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct IndependentPilotPlanHeader {
    /// Frame-invariant allocation policy. The per-frame retained path total is
    /// deliberately excluded so adjacent frames can share one sequence
    /// identity while still proving their own exact count-plane total.
    sequence_policy: String,
    minimum_samples: u32,
    total_samples: u64,
}

fn validate_independent_pilot_plan_attribute(
    exr: &DecodedExr,
    sample_ceiling: u32,
    path: &Path,
) -> Result<IndependentPilotPlanHeader, String> {
    let policy = exr_attribute_string(exr, "frankensim.render.pilotPlan", path)?;
    let mut fields = policy.split(';');
    let version = fields
        .next()
        .and_then(|field| field.strip_prefix("version="));
    let pilot_seed = fields
        .next()
        .and_then(|field| field.strip_prefix("pilotSeed="));
    let pilot_sampler = fields
        .next()
        .and_then(|field| field.strip_prefix("pilotSampler="));
    let minimum = fields
        .next()
        .and_then(|field| field.strip_prefix("minimum="));
    let maximum = fields
        .next()
        .and_then(|field| field.strip_prefix("maximum="));
    let absolute = fields
        .next()
        .and_then(|field| field.strip_prefix("absolute="));
    let relative = fields
        .next()
        .and_then(|field| field.strip_prefix("relative="));
    let dark_floor = fields
        .next()
        .and_then(|field| field.strip_prefix("darkFloor="));
    let safety = fields
        .next()
        .and_then(|field| field.strip_prefix("safety="));
    let total = fields.next().and_then(|field| field.strip_prefix("total="));
    if fields.next().is_some()
        || canonical_u32(version) != Some(INDEPENDENT_PILOT_ALLOCATION_SEMANTICS_VERSION)
    {
        return Err(invalid_independent_pilot_plan(path));
    }
    let pilot_seed =
        canonical_u64(pilot_seed).ok_or_else(|| invalid_independent_pilot_plan(path))?;
    let pilot_sampler = pilot_sampler
        .filter(|value| !value.is_empty())
        .ok_or_else(|| invalid_independent_pilot_plan(path))?;
    let minimum = canonical_u32(minimum).ok_or_else(|| invalid_independent_pilot_plan(path))?;
    let maximum = canonical_u32(maximum).ok_or_else(|| invalid_independent_pilot_plan(path))?;
    if minimum < 2 || minimum > maximum || maximum != sample_ceiling {
        return Err(invalid_independent_pilot_plan(path));
    }
    for value in [absolute, relative, dark_floor] {
        let value = value.ok_or_else(|| invalid_independent_pilot_plan(path))?;
        let bits =
            canonical_f64_string_bits(value).ok_or_else(|| invalid_independent_pilot_plan(path))?;
        if f64::from_bits(bits) < 0.0 {
            return Err(invalid_independent_pilot_plan(path));
        }
    }
    let safety = safety.ok_or_else(|| invalid_independent_pilot_plan(path))?;
    let safety_bits =
        canonical_f64_string_bits(safety).ok_or_else(|| invalid_independent_pilot_plan(path))?;
    if f64::from_bits(safety_bits) <= 0.0 {
        return Err(invalid_independent_pilot_plan(path));
    }
    let total_samples = canonical_u64(total).ok_or_else(|| invalid_independent_pilot_plan(path))?;
    if total_samples == 0 {
        return Err(invalid_independent_pilot_plan(path));
    }
    let sequence_policy = format!(
        "version={INDEPENDENT_PILOT_ALLOCATION_SEMANTICS_VERSION};pilotSeed={pilot_seed};pilotSampler={pilot_sampler};minimum={minimum};maximum={maximum};absolute={};relative={};darkFloor={};safety={safety}",
        absolute.expect("validated absolute field"),
        relative.expect("validated relative field"),
        dark_floor.expect("validated dark-floor field"),
    );
    Ok(IndependentPilotPlanHeader {
        sequence_policy,
        minimum_samples: minimum,
        total_samples,
    })
}

fn canonical_u32(value: Option<&str>) -> Option<u32> {
    let value = value?;
    let parsed = value.parse::<u32>().ok()?;
    (parsed.to_string() == value).then_some(parsed)
}

fn canonical_u64(value: Option<&str>) -> Option<u64> {
    let value = value?;
    let parsed = value.parse::<u64>().ok()?;
    (parsed.to_string() == value).then_some(parsed)
}

fn canonical_f64_string_bits(encoded: &str) -> Option<u64> {
    let (_, hex_bits) = encoded.split_once("@0x")?;
    if hex_bits.len() != 16
        || !hex_bits
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    {
        return None;
    }
    let bits = u64::from_str_radix(hex_bits, 16).ok()?;
    let value = f64::from_bits(bits);
    if !value.is_finite() {
        return None;
    }
    let canonical = if value == 0.0 { 0.0 } else { value };
    (encoded == format!("{canonical}@0x{:016x}", canonical.to_bits()))
        .then_some(canonical.to_bits())
}

fn invalid_adaptive_policy(path: &Path) -> String {
    format!(
        "FinalDiagnostic EXR {} attribute frankensim.render.adaptive is not the canonical adaptive policy grammar",
        path.display()
    )
}

fn invalid_independent_pilot_plan(path: &Path) -> String {
    format!(
        "FinalDiagnostic EXR {} attribute frankensim.render.pilotPlan is not the canonical independent-pilot plan grammar",
        path.display()
    )
}

fn nonempty_string_attribute(exr: &DecodedExr, name: &str, path: &Path) -> Result<String, String> {
    let value = exr_attribute_string(exr, name, path)?;
    if value.is_empty() {
        return Err(format!(
            "FinalDiagnostic EXR {} attribute {name} must not be empty",
            path.display()
        ));
    }
    Ok(value)
}

fn validate_frame_timing(exr: &DecodedExr, path: &Path) -> Result<FrameTiming, String> {
    let timing = FrameTiming {
        frame_time_bits: canonical_f64_bits_attribute(exr, "frankensim.frame.timeSeconds", path)?,
        previous_time_bits: canonical_f64_bits_attribute(
            exr,
            "frankensim.frame.previousTimeS",
            path,
        )?,
        next_time_bits: canonical_f64_bits_attribute(exr, "frankensim.frame.nextTimeS", path)?,
        shutter_open_bits: canonical_f64_bits_attribute(
            exr,
            "frankensim.render.shutterOpenS",
            path,
        )?,
        shutter_close_bits: canonical_f64_bits_attribute(
            exr,
            "frankensim.render.shutterCloseS",
            path,
        )?,
    };
    let frame = f64::from_bits(timing.frame_time_bits);
    let previous = f64::from_bits(timing.previous_time_bits);
    let next = f64::from_bits(timing.next_time_bits);
    let shutter_open = f64::from_bits(timing.shutter_open_bits);
    let shutter_close = f64::from_bits(timing.shutter_close_bits);
    if previous < 0.0
        || previous > frame
        || frame > next
        || previous > shutter_open
        || shutter_open > shutter_close
        || frame > shutter_close
        || shutter_close > next
    {
        return Err(format!(
            "FinalDiagnostic EXR {} frame/motion/shutter times are not coherently ordered",
            path.display()
        ));
    }
    Ok(timing)
}

fn validate_timing_continuity(
    previous: Option<FrameTiming>,
    current: FrameTiming,
    frame_index: u64,
    boundary: TemporalFrameBoundary,
) -> Result<(), String> {
    if boundary == TemporalFrameBoundary::Cut {
        if frame_index != 0 {
            return Ok(());
        }
        if current.frame_time_bits != 0 || current.previous_time_bits != 0 {
            return Err(
                "frame continuity violation at absolute frame 0: presentation and previous-reference clocks must begin at canonical +0"
                    .to_owned(),
            );
        }
        return Ok(());
    }
    if frame_index == 0 {
        return Err(
            "frame continuity violation at absolute frame 0: the first requested frame must be a temporal-history cut"
                .to_owned(),
        );
    }
    let previous = previous.ok_or_else(|| {
        format!(
            "frame continuity violation at absolute frame {frame_index}: prior timing is unavailable"
        )
    })?;
    if current.frame_time_bits != previous.next_time_bits
        || current.previous_time_bits != previous.frame_time_bits
    {
        return Err(format!(
            "frame continuity violation at absolute frame {frame_index}: presentation/motion-reference clocks do not link exactly to the preceding frame"
        ));
    }
    Ok(())
}

fn canonical_f64_bits_attribute(exr: &DecodedExr, name: &str, path: &Path) -> Result<u64, String> {
    let encoded = exr_attribute_string(exr, name, path)?;
    canonical_f64_string_bits(&encoded).ok_or_else(|| invalid_f64_attribute(path, name))
}

fn invalid_f64_attribute(path: &Path, name: &str) -> String {
    format!(
        "FinalDiagnostic EXR {} attribute {name} is not a finite canonical decimal@0xIEEE754 value",
        path.display()
    )
}

fn require_string_attribute(
    exr: &DecodedExr,
    name: &str,
    expected: &str,
    path: &Path,
) -> Result<(), String> {
    let value = exr_attribute_string(exr, name, path)?;
    if value != expected {
        return Err(format!(
            "FinalDiagnostic EXR {} attribute {name} is {value:?}; expected frozen value {expected:?}",
            path.display()
        ));
    }
    Ok(())
}

fn canonical_object_palette(
    exr: &DecodedExr,
    name: &str,
    path: &Path,
) -> Result<(String, u64), String> {
    let value = exr_attribute_string(exr, name, path)?;
    let mut rows = value.split(';');
    if rows.next() != Some("0=unavailable") {
        return Err(invalid_palette(path, name));
    }
    let mut count = 0_u64;
    let mut previous = None;
    for row in rows {
        count = count
            .checked_add(1)
            .ok_or_else(|| invalid_palette(path, name))?;
        if count >= MAX_EXACT_F32_INTEGER {
            return Err(invalid_palette(path, name));
        }
        let (index, object_id) = row
            .split_once('=')
            .ok_or_else(|| invalid_palette(path, name))?;
        if !is_canonical_decimal(index) || index.parse::<u64>().ok() != Some(count) {
            return Err(invalid_palette(path, name));
        }
        if !is_canonical_decimal(object_id) {
            return Err(invalid_palette(path, name));
        }
        let object_id = object_id
            .parse::<u64>()
            .map_err(|_| invalid_palette(path, name))?;
        if object_id == 0 || previous.is_some_and(|prior| object_id <= prior) {
            return Err(invalid_palette(path, name));
        }
        previous = Some(object_id);
    }
    Ok((value, count))
}

fn canonical_material_palette(
    exr: &DecodedExr,
    name: &str,
    path: &Path,
) -> Result<(String, u64), String> {
    let value = exr_attribute_string(exr, name, path)?;
    let mut rows = value.split(';');
    if rows.next() != Some("0=unavailable") {
        return Err(invalid_palette(path, name));
    }
    let mut count = 0_u64;
    let mut previous: Option<&str> = None;
    for row in rows {
        count = count
            .checked_add(1)
            .ok_or_else(|| invalid_palette(path, name))?;
        if count >= MAX_EXACT_F32_INTEGER {
            return Err(invalid_palette(path, name));
        }
        let (index, material_hash) = row
            .split_once('=')
            .ok_or_else(|| invalid_palette(path, name))?;
        if !is_canonical_decimal(index)
            || index.parse::<u64>().ok() != Some(count)
            || !is_canonical_content_hash(material_hash)
            || previous.is_some_and(|prior| material_hash <= prior)
        {
            return Err(invalid_palette(path, name));
        }
        previous = Some(material_hash);
    }
    Ok((value, count))
}

fn invalid_palette(path: &Path, name: &str) -> String {
    format!(
        "FinalDiagnostic EXR {} attribute {name} is not the producer's canonical sorted one-based palette grammar",
        path.display()
    )
}

fn is_canonical_decimal(value: &str) -> bool {
    !value.is_empty()
        && value.bytes().all(|byte| byte.is_ascii_digit())
        && (value == "0" || !value.starts_with('0'))
}

fn content_hash_attribute(exr: &DecodedExr, name: &str, path: &Path) -> Result<String, String> {
    let value = exr_attribute_string(exr, name, path)?;
    if !is_canonical_content_hash(&value) {
        return Err(format!(
            "FinalDiagnostic EXR {} attribute {name} is not a nonzero canonical 64-digit lowercase content hash",
            path.display()
        ));
    }
    Ok(value)
}

fn is_canonical_content_hash(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
        && value.bytes().any(|byte| byte != b'0')
}

fn exr_attribute_string(exr: &DecodedExr, name: &str, path: &Path) -> Result<String, String> {
    let attribute = required_exr_attribute(exr, name, path)?;
    if attribute.ty != "string" {
        return Err(format!(
            "FinalDiagnostic EXR {} attribute {name} has type {:?}; expected string",
            path.display(),
            attribute.ty
        ));
    }
    String::from_utf8(attribute.value.clone()).map_err(|_| {
        format!(
            "FinalDiagnostic EXR {} attribute {name} is not valid UTF-8 string data",
            path.display()
        )
    })
}

fn exr_attribute_u64(exr: &DecodedExr, name: &str, path: &Path) -> Result<u64, String> {
    let attribute = required_exr_attribute(exr, name, path)?;
    match attribute.ty.as_str() {
        "string" => std::str::from_utf8(&attribute.value)
            .map_err(|_| {
                format!(
                    "FinalDiagnostic EXR {} attribute {name} is not valid UTF-8 string data",
                    path.display()
                )
            })?
            .parse::<u64>()
            .map_err(|_| {
                format!(
                    "FinalDiagnostic EXR {} attribute {name} is not a decimal u64",
                    path.display()
                )
            }),
        "uint64" => {
            let bytes: [u8; 8] = attribute.value.as_slice().try_into().map_err(|_| {
                format!(
                    "FinalDiagnostic EXR {} attribute {name} uint64 payload is not exactly eight bytes",
                    path.display()
                )
            })?;
            Ok(u64::from_le_bytes(bytes))
        }
        ty => Err(format!(
            "FinalDiagnostic EXR {} attribute {name} has type {ty:?}; expected string or uint64",
            path.display()
        )),
    }
}

fn required_exr_attribute<'a>(
    exr: &'a DecodedExr,
    name: &str,
    path: &Path,
) -> Result<&'a ExrAttribute, String> {
    let mut found = exr
        .attributes
        .iter()
        .filter(|attribute| attribute.name == name);
    let attribute = found.next().ok_or_else(|| {
        format!(
            "FinalDiagnostic EXR {} is missing required attribute {name}",
            path.display()
        )
    })?;
    if found.next().is_some() {
        return Err(format!(
            "FinalDiagnostic EXR {} has duplicate required attribute {name}",
            path.display()
        ));
    }
    Ok(attribute)
}

fn validate_dimensions(width: u32, height: u32, path: &Path) -> Result<(), String> {
    let pixels = u64::from(width)
        .checked_mul(u64::from(height))
        .ok_or_else(|| {
            format!(
                "FinalDiagnostic EXR {} dimension product overflows",
                path.display()
            )
        })?;
    if width == 0 || height == 0 || width > 3_840 || height > 2_160 || pixels > MAX_4K_PIXELS {
        return Err(format!(
            "FinalDiagnostic EXR {} dimensions {}x{} exceed the bounded 4K envelope",
            path.display(),
            width,
            height,
        ));
    }
    Ok(())
}

fn offline_preview_manifest(
    cli: &Cli,
    source: &SequenceIdentity,
    dimensions: (u32, u32),
    raw_sequence_identity: ContentHash,
    preview_sequence_identity: ContentHash,
    spp_transition: Option<&UniformSppTransition>,
    color: CinematicColorConfig,
    denoise: TemporalDenoiseConfig,
) -> Result<String, String> {
    let denoise_config = denoise
        .identity()
        .map_err(|error| format!("identify temporal denoiser configuration: {error}"))?;
    let denoise_config_identity = hash_domain(
        "org.frankensim.euler-critique.temporal-denoise-config.v1",
        denoise_config.as_bytes(),
    );
    let color_bytes = color
        .canonical_bytes()
        .map_err(|error| format!("identify cinematic color configuration: {error}"))?;
    let color_identity = hash_domain(
        "org.frankensim.euler-critique.display-color-config.v1",
        &color_bytes,
    );
    let object_palette_identity = hash_domain(
        "org.frankensim.euler-critique.object-palette.v1",
        source.object_palette.as_bytes(),
    );
    let material_palette_identity = hash_domain(
        "org.frankensim.euler-critique.material-palette.v1",
        source.material_palette.as_bytes(),
    );
    let mut render_contract =
        DomainHasher::new("org.frankensim.euler-critique.source-render-contract.v2");
    render_contract.update(&source.shot_id.to_le_bytes());
    render_contract.update(source.cut_side.as_bytes());
    render_contract.update(source.shutter.as_bytes());
    render_contract.update(source.sampler.as_bytes());
    render_contract.update(source.strategy.as_bytes());
    render_contract.update(&source.max_depth.to_le_bytes());
    render_contract.update(source.render_versions.as_bytes());
    render_contract.update(source.sample_mode.as_bytes());
    render_contract.update(&source.sample_ceiling.to_le_bytes());
    if let Some(policy) = &source.adaptive_policy {
        render_contract.update(policy.as_bytes());
    }
    if let Some(policy) = &source.independent_pilot_policy {
        render_contract.update(policy.as_bytes());
    }
    if let Some(transition) = spp_transition {
        render_contract.update(b"authorized-uniform-spp-transition-v1");
        render_contract.update(&transition.frame.to_le_bytes());
        render_contract.update(&transition.from_spp.to_le_bytes());
        render_contract.update(&transition.to_spp.to_le_bytes());
        render_contract.update(transition.from_composition.as_bytes());
        render_contract.update(transition.to_composition.as_bytes());
        render_contract.update(transition.from_shutter.as_bytes());
        render_contract.update(transition.to_shutter.as_bytes());
    }
    let render_contract_identity = render_contract.finalize();
    let mut derivative =
        DomainHasher::new("org.frankensim.euler-critique.offline-denoised-preview-sequence.v2");
    derivative.update(raw_sequence_identity.as_bytes());
    derivative.update(preview_sequence_identity.as_bytes());
    derivative.update(denoise_config_identity.as_bytes());
    derivative.update(color_identity.as_bytes());
    derivative.update(render_contract_identity.as_bytes());
    derivative.update(&cli.frame_start.to_le_bytes());
    derivative.update(&cli.frame_count.to_le_bytes());
    derivative.update(&dimensions.0.to_le_bytes());
    derivative.update(&dimensions.1.to_le_bytes());
    let derivative_identity = derivative.finalize();
    let adaptive_policy_json = source
        .adaptive_policy
        .as_ref()
        .map_or_else(|| "null".to_owned(), |policy| format!("\"{policy}\""));
    let independent_pilot_policy_json = source
        .independent_pilot_policy
        .as_ref()
        .map_or_else(|| "null".to_owned(), |policy| format!("\"{policy}\""));
    let spp_transition_json = spp_transition.map_or_else(
        || "null".to_owned(),
        |transition| {
            format!(
                concat!(
                    "{{\"frame\":{},\"authority\":\"explicit-cli-opt-in\",",
                    "\"from_uniform_spp\":{},\"to_uniform_spp\":{},",
                    "\"from_composition_identity\":\"{}\",",
                    "\"to_composition_identity\":\"{}\",",
                    "\"from_shutter\":\"{}\",\"to_shutter\":\"{}\"}}"
                ),
                transition.frame,
                transition.from_spp,
                transition.to_spp,
                transition.from_composition,
                transition.to_composition,
                transition.from_shutter,
                transition.to_shutter,
            )
        },
    );

    Ok(format!(
        concat!(
            "{{\n",
            "  \"schema\": \"frankensim-euler-offline-denoise-v4\",\n",
            "  \"authority\": \"biased-display-derivative-not-raw-estimate\",\n",
            "  \"publication\": \"single-writer-staged-directory-rename\",\n",
            "  \"initial_boundary\": \"cut\",\n",
            "  \"nonzero_initial_cut_authorized\": {},\n",
            "  \"frame_start\": {},\n",
            "  \"frame_count\": {},\n",
            "  \"width\": {},\n",
            "  \"height\": {},\n",
            "  \"source_trajectory_identity\": \"{}\",\n",
            "  \"source_scene_identity\": \"{}\",\n",
            "  \"source_composition_identity\": \"{}\",\n",
            "  \"source_aov_profile\": \"{}\",\n",
            "  \"source_aov_semantics_version\": {},\n",
            "  \"source_sample_mode\": \"{}\",\n",
            "  \"source_sample_ceiling\": {},\n",
            "  \"source_adaptive_policy\": {},\n",
            "  \"source_independent_pilot_policy\": {},\n",
            "  \"authorized_uniform_spp_transition\": {},\n",
            "  \"source_object_palette_entries\": {},\n",
            "  \"source_material_palette_entries\": {},\n",
            "  \"source_object_palette_identity\": \"{}\",\n",
            "  \"source_material_palette_identity\": \"{}\",\n",
            "  \"source_raw_sequence_identity\": \"{}\",\n",
            "  \"source_render_contract_identity\": \"{}\",\n",
            "  \"temporal_pipeline_version\": \"{}\",\n",
            "  \"spatial_iterations\": {},\n",
            "  \"spatial_sigma_rgb\": {},\n",
            "  \"temporal_config_identity\": \"{}\",\n",
            "  \"display_color_config_identity\": \"{}\",\n",
            "  \"preview_sequence_identity\": \"{}\",\n",
            "  \"derivative_identity\": \"{}\"\n",
            "}}\n"
        ),
        cli.initial_cut && cli.frame_start != 0,
        cli.frame_start,
        cli.frame_count,
        dimensions.0,
        dimensions.1,
        source.source_trajectory,
        source.scene_hash,
        source.composition,
        source.aov_profile,
        CINEMATIC_AOV_SEMANTICS_VERSION,
        source.sample_mode,
        source.sample_ceiling,
        adaptive_policy_json,
        independent_pilot_policy_json,
        spp_transition_json,
        source.object_palette_entries,
        source.material_palette_entries,
        object_palette_identity.to_hex(),
        material_palette_identity.to_hex(),
        raw_sequence_identity.to_hex(),
        render_contract_identity.to_hex(),
        TEMPORAL_DENOISE_PIPELINE_VERSION,
        denoise.spatial_iterations,
        denoise.spatial_sigma_rgb,
        denoise_config_identity.to_hex(),
        color_identity.to_hex(),
        preview_sequence_identity.to_hex(),
        derivative_identity.to_hex(),
    ))
}

fn raw_path(input: &Path, frame: u64) -> PathBuf {
    input.join(format!("frame-{frame:06}.exr"))
}

fn preview_path(output: &Path, frame: u64) -> PathBuf {
    output.join(format!("frame-{frame:06}.png"))
}

fn staging_output_path(output: &Path) -> Result<PathBuf, String> {
    let file_name = output.file_name().ok_or_else(|| {
        format!(
            "output must name a directory rather than a filesystem root: {}",
            output.display()
        )
    })?;
    let mut staging_name = OsString::from(".");
    staging_name.push(file_name);
    staging_name.push(format!(".incomplete-{}", std::process::id()));
    Ok(output.with_file_name(staging_name))
}

fn require_absent(path: &Path, purpose: &str) -> Result<(), String> {
    match fs::symlink_metadata(path) {
        Ok(_) => Err(format!(
            "refusing to overwrite existing {purpose}: {}",
            path.display()
        )),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(format!("inspect {purpose} {}: {error}", path.display())),
    }
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|error| format!("create preview {}: {error}", path.display()))?;
    file.write_all(bytes)
        .map_err(|error| format!("write preview {}: {error}", path.display()))?;
    file.flush()
        .map_err(|error| format!("flush preview {}: {error}", path.display()))?;
    Ok(())
}

fn parse_cli(args: impl Iterator<Item = String>) -> Result<Command, String> {
    let mut input = None;
    let mut output = None;
    let mut frame_start = None;
    let mut frame_count = None;
    let mut initial_cut = false;
    let mut allow_uniform_spp_transition_at = None;
    let mut denoise_spatial_passes = None;
    let mut denoise_spatial_sigma = None;
    let mut args = args;
    while let Some(argument) = args.next() {
        match argument.as_str() {
            "--inspect-samples" => {
                if input.is_some()
                    || output.is_some()
                    || frame_start.is_some()
                    || frame_count.is_some()
                    || initial_cut
                    || allow_uniform_spp_transition_at.is_some()
                    || denoise_spatial_passes.is_some()
                    || denoise_spatial_sigma.is_some()
                {
                    return Err(format!(
                        "--inspect-samples cannot be combined with denoise arguments\n{}",
                        usage()
                    ));
                }
                let path = PathBuf::from(next_value(&mut args, "--inspect-samples")?);
                if let Some(extra) = args.next() {
                    return Err(format!(
                        "unexpected argument after --inspect-samples EXR: {extra}\n{}",
                        usage()
                    ));
                }
                return Ok(Command::InspectSamples(path));
            }
            "--input" => input = Some(PathBuf::from(next_value(&mut args, "--input")?)),
            "--output" => output = Some(PathBuf::from(next_value(&mut args, "--output")?)),
            "--frame-start" => {
                frame_start = Some(parse_u64(
                    &next_value(&mut args, "--frame-start")?,
                    "frame-start",
                )?)
            }
            "--frame-count" => {
                frame_count = Some(parse_u64(
                    &next_value(&mut args, "--frame-count")?,
                    "frame-count",
                )?)
            }
            "--initial-cut" => initial_cut = true,
            "--allow-uniform-spp-transition-at" => {
                allow_uniform_spp_transition_at = Some(parse_u64(
                    &next_value(&mut args, "--allow-uniform-spp-transition-at")?,
                    "allow-uniform-spp-transition-at",
                )?)
            }
            "--denoise-spatial-passes" => {
                denoise_spatial_passes = Some(parse_u8(
                    &next_value(&mut args, "--denoise-spatial-passes")?,
                    "denoise-spatial-passes",
                )?)
            }
            "--denoise-spatial-sigma" => {
                denoise_spatial_sigma = Some(parse_positive_f32(
                    &next_value(&mut args, "--denoise-spatial-sigma")?,
                    "denoise-spatial-sigma",
                )?)
            }
            "--help" | "-h" => return Err(usage().to_owned()),
            _ => return Err(format!("unknown argument: {argument}\n{}", usage())),
        }
    }
    let defaults = TemporalDenoiseConfig::default();
    let cli = Cli {
        input: input.ok_or_else(|| format!("missing --input\n{}", usage()))?,
        output: output.ok_or_else(|| format!("missing --output\n{}", usage()))?,
        frame_start: frame_start.ok_or_else(|| format!("missing --frame-start\n{}", usage()))?,
        frame_count: frame_count.ok_or_else(|| format!("missing --frame-count\n{}", usage()))?,
        initial_cut,
        allow_uniform_spp_transition_at,
        denoise_spatial_passes: denoise_spatial_passes.unwrap_or(defaults.spatial_iterations),
        denoise_spatial_sigma: denoise_spatial_sigma.unwrap_or(defaults.spatial_sigma_rgb),
    };
    if cli.frame_count == 0 {
        return Err("frame-count must be positive".to_owned());
    }
    if cli.frame_start != 0 && !cli.initial_cut {
        return Err(format!(
            "nonzero frame-start {} requires --initial-cut so missing temporal history is explicit",
            cli.frame_start
        ));
    }
    let frame_end = cli
        .frame_start
        .checked_add(cli.frame_count)
        .ok_or_else(|| "frame range overflows u64".to_owned())?;
    if let Some(transition_frame) = cli.allow_uniform_spp_transition_at
        && (transition_frame <= cli.frame_start || transition_frame >= frame_end)
    {
        return Err(format!(
            "allow-uniform-spp-transition-at must identify a frame strictly inside the requested range {}..{}; found {transition_frame}",
            cli.frame_start, frame_end
        ));
    }
    if cli.denoise_spatial_passes == 0
        || cli.denoise_spatial_passes > MAX_TEMPORAL_DENOISE_SPATIAL_ITERATIONS
    {
        return Err(format!(
            "denoise-spatial-passes must be in 1..={MAX_TEMPORAL_DENOISE_SPATIAL_ITERATIONS}; found {}",
            cli.denoise_spatial_passes
        ));
    }
    Ok(Command::Denoise(cli))
}

fn next_value(args: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    args.next()
        .ok_or_else(|| format!("missing value for {flag}"))
}

fn parse_u64(value: &str, name: &str) -> Result<u64, String> {
    value
        .parse()
        .map_err(|_| format!("invalid {name}: {value}"))
}

fn parse_u8(value: &str, name: &str) -> Result<u8, String> {
    value
        .parse()
        .map_err(|_| format!("invalid {name}: {value}"))
}

fn parse_positive_f32(value: &str, name: &str) -> Result<f32, String> {
    let parsed: f32 = value
        .parse()
        .map_err(|_| format!("invalid {name}: {value}"))?;
    if !parsed.is_finite() || parsed <= 0.0 {
        return Err(format!("{name} must be finite and positive; found {value}"));
    }
    Ok(parsed)
}

const fn usage() -> &'static str {
    "Usage:\n  euler_cinematic_denoise --input DIR --output DIR --frame-start N --frame-count N [--initial-cut] [--allow-uniform-spp-transition-at FRAME] [--denoise-spatial-passes 1..8] [--denoise-spatial-sigma X]\n  euler_cinematic_denoise --inspect-samples EXR"
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set_attribute(exr: &mut DecodedExr, name: &str, value: &str) {
        exr.attributes
            .iter_mut()
            .find(|attribute| attribute.name == name)
            .unwrap_or_else(|| panic!("fixture attribute {name}"))
            .value = value.as_bytes().to_vec();
    }

    fn final_diagnostic() -> DecodedExr {
        let schema_version = CINEMATIC_AOV_SEMANTICS_VERSION.to_string();
        let channels = FINAL_DIAGNOSTIC_CHANNELS
            .iter()
            .map(|(name, ty)| Channel {
                name: (*name).to_owned(),
                ty: *ty,
                data: match *name {
                    "samples" => vec![16.0, 16.0],
                    "primary.coverage" => vec![1.0, 1.0],
                    "motion.prev.X" | "motion.prev.Y" | "normal.X" | "normal.Y" => {
                        vec![0.0, 0.0]
                    }
                    "normal.Z" => vec![1.0, 1.0],
                    _ => vec![1.0, 2.0],
                },
            })
            .collect();
        DecodedExr {
            width: 2,
            height: 1,
            channels,
            attributes: [
                ("frankensim.aov.authority", "raw-estimate"),
                ("frankensim.aov.schemaVersion", schema_version.as_str()),
                (
                    "frankensim.aov.channelSemantics",
                    CINEMATIC_AOV_CHANNEL_SEMANTICS,
                ),
                (
                    "frankensim.aov.invalidSemantics",
                    CINEMATIC_AOV_INVALID_SEMANTICS,
                ),
                (
                    "frankensim.aov.materialDomain",
                    MATERIAL_CONTENT_IDENTITY_DOMAIN,
                ),
                (
                    "frankensim.aov.paletteZero",
                    CINEMATIC_AOV_PALETTE_ZERO_SEMANTICS,
                ),
                (
                    "frankensim.aov.objectPalette",
                    "0=unavailable;1=101;2=202",
                ),
                (
                    "frankensim.aov.materialPalette",
                    "0=unavailable;1=5555555555555555555555555555555555555555555555555555555555555555;2=6666666666666666666666666666666666666666666666666666666666666666",
                ),
                ("frankensim.frame.index", "0"),
                ("frankensim.frame.timeSeconds", "0@0x0000000000000000"),
                ("frankensim.frame.previousTimeS", "0@0x0000000000000000"),
                ("frankensim.frame.nextTimeS", "1@0x3ff0000000000000"),
                (
                    "frankensim.source.trajectory",
                    "1111111111111111111111111111111111111111111111111111111111111111",
                ),
                (
                    "frankensim.source.sceneHash",
                    "2222222222222222222222222222222222222222222222222222222222222222",
                ),
                (
                    "frankensim.source.composition",
                    "3333333333333333333333333333333333333333333333333333333333333333",
                ),
                ("frankensim.aov.profile", FINAL_DIAGNOSTIC_PROFILE),
                (
                    "frankensim.aov.configHash",
                    "4444444444444444444444444444444444444444444444444444444444444444",
                ),
                ("frankensim.render.sampleMode", "uniform"),
                ("frankensim.render.spp", "16"),
                ("frankensim.render.sppCeiling", "16"),
                ("frankensim.render.shotId", "7"),
                ("frankensim.render.cutSide", "after"),
                (
                    "frankensim.render.shutter",
                    "convention=back-loaded;distribution=stratified-counter-v1;strata=16",
                ),
                ("frankensim.render.shutterOpenS", "0@0x0000000000000000"),
                ("frankensim.render.shutterCloseS", "1@0x3ff0000000000000"),
                ("frankensim.render.sampler", "sobol-owen-v1"),
                ("frankensim.render.strategy", "next-event-mis"),
                ("frankensim.render.maxDepth", "8"),
                ("frankensim.render.versions", "fixture-v1"),
            ]
            .into_iter()
            .map(|(name, value)| ExrAttribute {
                name: name.to_owned(),
                ty: "string".to_owned(),
                value: value.as_bytes().to_vec(),
            })
            .collect(),
        }
    }

    #[test]
    fn final_diagnostic_schema_reconstructs_exact_denoise_planes_and_ids() {
        let frame = decode_final_diagnostic(final_diagnostic(), Path::new("fixture.exr"), 0)
            .expect("FinalDiagnostic");
        assert_eq!(frame.temporal_input(42).frame_index, 42);
        assert_eq!(frame.temporal_input(42).samples_per_pixel, 16);
        assert_eq!(
            frame.temporal_input(42).sample_counts_per_pixel,
            Some(&[16, 16][..])
        );
        assert_eq!(frame.temporal_input(42).red, &[1.0, 2.0]);
        assert_eq!(frame.temporal_input(42).variance_luminance, &[1.0, 2.0]);
        assert_eq!(frame.temporal_input(42).object_ids, Some(&[1, 2][..]));
        assert_eq!(frame.temporal_input(42).material_ids, Some(&[1, 2][..]));
    }

    #[test]
    fn missing_final_diagnostic_plane_refuses_before_denoising() {
        let mut exr = final_diagnostic();
        exr.channels
            .retain(|channel| channel.name != "motion.prev.Y");
        let error =
            decode_final_diagnostic(exr, Path::new("fixture.exr"), 0).expect_err("missing plane");
        assert!(error.contains("missing required channel motion.prev.Y"));
    }

    #[test]
    fn header_frame_index_must_match_the_absolute_requested_frame() {
        let error = decode_final_diagnostic(final_diagnostic(), Path::new("fixture.exr"), 1)
            .expect_err("mismatched metadata must refuse temporal history");
        assert!(error.contains("declares frame index 0; expected 1"));
    }

    #[test]
    fn uint64_frame_index_attribute_is_accepted_when_it_matches() {
        let mut exr = final_diagnostic();
        let attribute = exr
            .attributes
            .iter_mut()
            .find(|attribute| attribute.name == "frankensim.frame.index")
            .expect("fixture frame index");
        attribute.ty = "uint64".to_owned();
        attribute.value = 7_u64.to_le_bytes().to_vec();
        decode_final_diagnostic(exr, Path::new("fixture.exr"), 7)
            .expect("matching uint64 metadata must be accepted");
    }

    #[test]
    fn consecutive_frames_allow_distinct_frame_bound_aov_config_hashes() {
        let first = final_diagnostic();
        let first_identity =
            validate_final_diagnostic_attributes(&first, Path::new("frame-0.exr"), 0)
                .expect("first frame provenance");

        let mut second = final_diagnostic();
        let frame_index = second
            .attributes
            .iter_mut()
            .find(|attribute| attribute.name == "frankensim.frame.index")
            .expect("fixture frame index");
        frame_index.value = b"1".to_vec();
        let config_hash = second
            .attributes
            .iter_mut()
            .find(|attribute| attribute.name == "frankensim.aov.configHash")
            .expect("fixture AOV config hash");
        config_hash.value =
            b"5555555555555555555555555555555555555555555555555555555555555555".to_vec();

        let second_identity =
            validate_final_diagnostic_attributes(&second, Path::new("frame-1.exr"), 1)
                .expect("second frame provenance");
        assert_eq!(first_identity, second_identity);
    }

    #[test]
    fn frame_and_motion_reference_clocks_link_exactly_across_the_sequence() {
        let first = decode_final_diagnostic(final_diagnostic(), Path::new("frame-0.exr"), 0)
            .expect("first frame");
        validate_timing_continuity(None, first.timing, 0, TemporalFrameBoundary::Cut)
            .expect("canonical clock origin");

        let mut second_exr = final_diagnostic();
        set_attribute(&mut second_exr, "frankensim.frame.index", "1");
        set_attribute(
            &mut second_exr,
            "frankensim.frame.timeSeconds",
            "1@0x3ff0000000000000",
        );
        set_attribute(
            &mut second_exr,
            "frankensim.frame.previousTimeS",
            "0@0x0000000000000000",
        );
        set_attribute(
            &mut second_exr,
            "frankensim.frame.nextTimeS",
            "2@0x4000000000000000",
        );
        set_attribute(
            &mut second_exr,
            "frankensim.render.shutterOpenS",
            "1@0x3ff0000000000000",
        );
        set_attribute(
            &mut second_exr,
            "frankensim.render.shutterCloseS",
            "2@0x4000000000000000",
        );
        let second = decode_final_diagnostic(second_exr.clone(), Path::new("frame-1.exr"), 1)
            .expect("second frame");
        validate_timing_continuity(
            Some(first.timing),
            second.timing,
            1,
            TemporalFrameBoundary::Continuous,
        )
        .expect("exact current/previous/next linkage");
        validate_timing_continuity(None, second.timing, 1, TemporalFrameBoundary::Cut)
            .expect("an explicitly cut nonzero range validates its first clock locally");
        let error =
            validate_timing_continuity(None, second.timing, 1, TemporalFrameBoundary::Continuous)
                .expect_err("a continuous nonzero frame still requires its exact predecessor");
        assert!(error.contains("prior timing is unavailable"));

        set_attribute(
            &mut second_exr,
            "frankensim.frame.previousTimeS",
            "1@0x3ff0000000000000",
        );
        let second = decode_final_diagnostic(second_exr, Path::new("frame-1.exr"), 1)
            .expect("locally ordered but wrongly linked frame");
        let error = validate_timing_continuity(
            Some(first.timing),
            second.timing,
            1,
            TemporalFrameBoundary::Continuous,
        )
        .expect_err("wrong motion-reference clock must break sequence admission");
        assert!(error.contains("do not link exactly"));
    }

    #[test]
    fn noncanonical_or_incoherent_f64_clock_metadata_is_refused() {
        let mut noncanonical = final_diagnostic();
        set_attribute(
            &mut noncanonical,
            "frankensim.frame.nextTimeS",
            "1.0@0x3ff0000000000000",
        );
        let error = decode_final_diagnostic(noncanonical, Path::new("frame-0.exr"), 0)
            .expect_err("noncanonical clock string");
        assert!(error.contains("finite canonical decimal@0xIEEE754"));

        let mut incoherent = final_diagnostic();
        set_attribute(
            &mut incoherent,
            "frankensim.render.shutterCloseS",
            "2@0x4000000000000000",
        );
        let error = decode_final_diagnostic(incoherent, Path::new("frame-0.exr"), 0)
            .expect_err("shutter beyond next motion reference");
        assert!(error.contains("times are not coherently ordered"));
    }

    #[test]
    fn every_frame_still_requires_its_aov_config_hash() {
        let mut exr = final_diagnostic();
        exr.attributes
            .retain(|attribute| attribute.name != "frankensim.aov.configHash");
        let error = validate_final_diagnostic_attributes(&exr, Path::new("frame-0.exr"), 0)
            .expect_err("missing per-frame integrity metadata must be refused");
        assert!(error.contains("missing required attribute frankensim.aov.configHash"));
    }

    #[test]
    fn content_hashes_must_be_canonical_nonzero_lowercase_sha256_values() {
        for invalid in [
            "0".repeat(64),
            "a".repeat(63),
            "A".repeat(64),
            format!("{}g", "a".repeat(63)),
        ] {
            let mut exr = final_diagnostic();
            let hash = exr
                .attributes
                .iter_mut()
                .find(|attribute| attribute.name == "frankensim.source.trajectory")
                .expect("fixture trajectory hash");
            hash.value = invalid.into_bytes();
            let error = validate_final_diagnostic_attributes(&exr, Path::new("frame-0.exr"), 0)
                .expect_err("noncanonical content hash must be refused");
            assert!(error.contains("nonzero canonical 64-digit lowercase content hash"));
        }
    }

    #[test]
    fn frozen_aov_semantics_are_required_before_variance_is_interpreted() {
        for (attribute_name, replacement) in [
            ("frankensim.aov.authority", "filtered-estimate"),
            ("frankensim.aov.schemaVersion", "1"),
            ("frankensim.aov.channelSemantics", "variance.Y=unknown"),
            ("frankensim.aov.invalidSemantics", "nan"),
            ("frankensim.aov.materialDomain", "unknown-material-domain"),
            ("frankensim.aov.paletteZero", "0=object"),
        ] {
            let mut exr = final_diagnostic();
            exr.attributes
                .iter_mut()
                .find(|attribute| attribute.name == attribute_name)
                .expect("fixture semantic attribute")
                .value = replacement.as_bytes().to_vec();
            let error = validate_final_diagnostic_attributes(&exr, Path::new("frame-0.exr"), 0)
                .expect_err("changed AOV semantics must be refused");
            assert!(error.contains("expected frozen value"));
        }
    }

    #[test]
    fn palettes_are_canonical_sequence_identity_and_bound_id_ranges() {
        let first =
            validate_final_diagnostic_attributes(&final_diagnostic(), Path::new("frame-0.exr"), 0)
                .expect("first identity");
        let mut remapped = final_diagnostic();
        remapped
            .attributes
            .iter_mut()
            .find(|attribute| attribute.name == "frankensim.aov.objectPalette")
            .expect("fixture object palette")
            .value = b"0=unavailable;1=102;2=202".to_vec();
        let remapped_identity =
            validate_final_diagnostic_attributes(&remapped, Path::new("frame-0.exr"), 0)
                .expect("canonical but remapped identity");
        assert_ne!(first, remapped_identity);

        let mut malformed = final_diagnostic();
        malformed
            .attributes
            .iter_mut()
            .find(|attribute| attribute.name == "frankensim.aov.objectPalette")
            .expect("fixture object palette")
            .value = b"0=unavailable;2=101;1=202".to_vec();
        let error = validate_final_diagnostic_attributes(&malformed, Path::new("frame-0.exr"), 0)
            .expect_err("noncanonical palette grammar must be refused");
        assert!(error.contains("canonical sorted one-based palette grammar"));

        let mut out_of_range = final_diagnostic();
        out_of_range
            .channels
            .iter_mut()
            .find(|channel| channel.name == "id.object")
            .expect("fixture object IDs")
            .data[1] = 3.0;
        let error = decode_final_diagnostic(out_of_range, Path::new("frame-0.exr"), 0)
            .expect_err("ID beyond the declared palette must be refused");
        assert!(error.contains("above declared maximum 2"));
    }

    #[test]
    fn uniform_header_spp_must_match_every_sample_plane_value() {
        let mut exr = final_diagnostic();
        exr.channels
            .iter_mut()
            .find(|channel| channel.name == "samples")
            .expect("fixture samples plane")
            .data[1] = 15.0;
        let error = decode_final_diagnostic(exr, Path::new("fixture.exr"), 0)
            .expect_err("sample-count mismatch must be refused");
        assert!(error.contains("samples plane disagrees with uniform header SPP"));
    }

    #[test]
    fn adaptive_sample_metadata_and_exact_count_plane_are_admitted() {
        let mut adaptive = final_diagnostic();
        adaptive
            .attributes
            .iter_mut()
            .find(|attribute| attribute.name == "frankensim.render.sampleMode")
            .expect("fixture sample mode")
            .value = b"adaptive".to_vec();
        adaptive
            .attributes
            .iter_mut()
            .find(|attribute| attribute.name == "frankensim.render.spp")
            .expect("fixture SPP")
            .value = b"per-pixel-channel".to_vec();
        adaptive.attributes.push(ExrAttribute {
            name: "frankensim.render.adaptive".to_owned(),
            ty: "string".to_owned(),
            value: b"version=1;minimum=2;batch=2;absolute=0@0x0000000000000000;relative=0.1@0x3fb999999999999a;darkFloor=0.01@0x3f847ae147ae147b".to_vec(),
        });
        adaptive
            .channels
            .iter_mut()
            .find(|channel| channel.name == "samples")
            .expect("fixture samples plane")
            .data[1] = 8.0;
        let frame = decode_final_diagnostic(adaptive, Path::new("adaptive.exr"), 0)
            .expect("canonical adaptive metadata and count plane");
        assert_eq!(frame.sequence_identity.sample_mode, "adaptive");
        assert_eq!(frame.sequence_identity.sample_ceiling, 16);
        assert_eq!(frame.sample_counts, [16, 8]);
        assert_eq!(
            frame.temporal_input(0).sample_counts_per_pixel,
            Some(&[16, 8][..])
        );
    }

    #[test]
    fn independent_pilot_metadata_binds_policy_and_exact_count_total() {
        let mut pilot = final_diagnostic();
        set_attribute(
            &mut pilot,
            "frankensim.render.sampleMode",
            "independent-pilot-fixed-v1",
        );
        set_attribute(&mut pilot, "frankensim.render.spp", "per-pixel-channel");
        pilot.attributes.push(ExrAttribute {
            name: "frankensim.render.pilotPlan".to_owned(),
            ty: "string".to_owned(),
            value: b"version=1;pilotSeed=41;pilotSampler=owen-sobol-full-path-v1;minimum=2;maximum=16;absolute=0@0x0000000000000000;relative=0.1@0x3fb999999999999a;darkFloor=0.01@0x3f847ae147ae147b;safety=1.5@0x3ff8000000000000;total=24".to_vec(),
        });
        pilot
            .channels
            .iter_mut()
            .find(|channel| channel.name == "samples")
            .expect("fixture samples plane")
            .data[1] = 8.0;
        let frame = decode_final_diagnostic(pilot.clone(), Path::new("pilot.exr"), 0)
            .expect("canonical independent-pilot metadata and count plane");
        assert_eq!(
            frame.sequence_identity.sample_mode,
            "independent-pilot-fixed-v1"
        );
        assert_eq!(frame.sample_counts, [16, 8]);
        assert_eq!(
            frame.sequence_identity.independent_pilot_policy.as_deref(),
            Some(
                "version=1;pilotSeed=41;pilotSampler=owen-sobol-full-path-v1;minimum=2;maximum=16;absolute=0@0x0000000000000000;relative=0.1@0x3fb999999999999a;darkFloor=0.01@0x3f847ae147ae147b;safety=1.5@0x3ff8000000000000"
            )
        );

        set_attribute(
            &mut pilot,
            "frankensim.render.pilotPlan",
            "version=1;pilotSeed=41;pilotSampler=owen-sobol-full-path-v1;minimum=2;maximum=16;absolute=0@0x0000000000000000;relative=0.1@0x3fb999999999999a;darkFloor=0.01@0x3f847ae147ae147b;safety=1.5@0x3ff8000000000000;total=25",
        );
        let error = decode_final_diagnostic(pilot, Path::new("pilot.exr"), 0)
            .expect_err("declared pilot total must match the count plane");
        assert!(error.contains("samples plane total Some(24) disagrees with declared total 25"));
    }

    #[test]
    fn explicit_uniform_spp_transition_changes_only_sampling_bound_identity() {
        let from =
            validate_final_diagnostic_attributes(&final_diagnostic(), Path::new("frame-0.exr"), 0)
                .expect("source identity");
        let mut higher_spp = final_diagnostic();
        set_attribute(&mut higher_spp, "frankensim.render.spp", "32");
        set_attribute(&mut higher_spp, "frankensim.render.sppCeiling", "32");
        set_attribute(
            &mut higher_spp,
            "frankensim.render.shutter",
            "convention=back-loaded;distribution=stratified-counter-v1;strata=32",
        );
        set_attribute(
            &mut higher_spp,
            "frankensim.source.composition",
            "7777777777777777777777777777777777777777777777777777777777777777",
        );
        let to = validate_final_diagnostic_attributes(&higher_spp, Path::new("frame-1.exr"), 0)
            .expect("destination identity");
        let transition =
            admit_uniform_spp_transition(1, &from, &to).expect("authorized SPP-only rung");
        assert_eq!(transition.frame, 1);
        assert_eq!((transition.from_spp, transition.to_spp), (16, 32));

        let mut changed_scene = to.clone();
        changed_scene.scene_hash =
            "8888888888888888888888888888888888888888888888888888888888888888".to_owned();
        let error = admit_uniform_spp_transition(1, &from, &changed_scene)
            .expect_err("scene changes must not hide behind an SPP transition");
        assert!(error.contains("changed provenance beyond uniform SPP"));

        let mut changed_shutter = to.clone();
        changed_shutter.shutter =
            "convention=centered;distribution=stratified-counter-v1;strata=32".to_owned();
        let error = admit_uniform_spp_transition(1, &from, &changed_shutter)
            .expect_err("shutter semantics must remain fixed");
        assert!(error.contains("changed shutter semantics beyond"));

        let mut adaptive = to.clone();
        adaptive.sample_mode = "adaptive".to_owned();
        adaptive.adaptive_policy = Some("policy".to_owned());
        let error = admit_uniform_spp_transition(1, &from, &adaptive)
            .expect_err("adaptive policy changes need a different contract");
        assert!(error.contains("must remain uniform"));
    }

    #[test]
    fn sample_inspection_statistics_are_exact_and_use_nearest_rank_quantiles() {
        let summary = summarize_sample_counts([2, 2, 2, 2, 8, 8, 16, 16].into_iter(), 16)
            .expect("valid exact sample counts");
        assert_eq!(summary.pixels, 8);
        assert_eq!(summary.total_samples, 56);
        assert_eq!(summary.minimum, 2);
        assert_eq!(summary.maximum, 16);
        assert_eq!(summary.at_ceiling_spp_pixels, 2);
        assert_eq!(summary.histogram, BTreeMap::from([(2, 4), (8, 2), (16, 2)]));
        assert_eq!(nearest_rank_quantile(&summary, 10, 100), 2);
        assert_eq!(nearest_rank_quantile(&summary, 50, 100), 2);
        assert_eq!(nearest_rank_quantile(&summary, 75, 100), 8);
        assert_eq!(nearest_rank_quantile(&summary, 99, 100), 16);

        let variance = summarize_estimator_variance(&[1.0, 2.0], &[16, 8])
            .expect("finite nonnegative variance and positive counts");
        assert_eq!(variance.sample_variance_total.to_bits(), 3.0_f64.to_bits());
        assert_eq!(
            variance.estimator_variance_total.to_bits(),
            0.3125_f64.to_bits()
        );
        assert_eq!(variance.maximum_estimator_variance, 0.25);
        assert_eq!(variance.maximum_pixel_index, 1);
    }

    #[test]
    fn sample_inspection_report_exposes_exact_ratios_histogram_and_fixed_crops() {
        let frame = decode_final_diagnostic(final_diagnostic(), Path::new("fixture.exr"), 0)
            .expect("FinalDiagnostic fixture");
        let report =
            sample_inspection_report(&frame, 0, hash_domain("sample-inspection-test", b"fixture"))
                .expect("sample inspection report");
        assert!(
            report.starts_with("record=metadata schema=frankensim-euler-sample-inspection-v2 ")
        );
        assert!(report.contains("allocation_decision=unavailable-from-samples-plane"));
        assert!(report.contains(
            "record=uncertainty scope=full channel=variance.Y meaning=unbiased-raw-CIE-Y-sample-variance"
        ));
        assert!(report.contains(
            "record=summary scope=full pixels=2 total_samples=32 mean_numerator=32 mean_denominator=2 mean_spp=16.000000000"
        ));
        assert!(report.contains(
            "at_ceiling_spp_pixels=2 at_ceiling_spp_fraction_numerator=2 at_ceiling_spp_fraction_denominator=2 at_ceiling_spp_fraction=1.000000000"
        ));
        assert!(report.contains(
            "record=histogram scope=full spp=16 pixels=2 fraction_numerator=2 fraction_denominator=2 fraction=1.000000000"
        ));
        assert!(report.contains("record=crop name=disc x_min=0 x_max=2 y_min=0 y_max=1 pixels=2"));
        assert!(
            report
                .contains("record=crop name=front_glass x_min=0 x_max=2 y_min=0 y_max=1 pixels=2")
        );
        assert!(
            report
                .contains("record=crop name=right_glass x_min=1 x_max=2 y_min=0 y_max=1 pixels=1")
        );
        assert!(
            report.contains("record=crop name=background x_min=0 x_max=1 y_min=0 y_max=1 pixels=1")
        );
    }

    #[test]
    fn calibration_crops_scale_exactly_from_320p_to_960p() {
        assert_eq!(
            reference_crop(960, 540, 108, 212, 54, 118),
            PixelCrop {
                x_min: 324,
                x_max: 636,
                y_min: 162,
                y_max: 354,
            }
        );
        assert_eq!(
            reference_crop(960, 540, 10, 230, 105, 173),
            PixelCrop {
                x_min: 30,
                x_max: 690,
                y_min: 315,
                y_max: 519,
            }
        );
    }

    #[test]
    fn zero_or_inexact_sample_ceiling_is_refused() {
        let mut zero = final_diagnostic();
        for name in ["frankensim.render.spp", "frankensim.render.sppCeiling"] {
            zero.attributes
                .iter_mut()
                .find(|attribute| attribute.name == name)
                .expect("fixture SPP")
                .value = b"0".to_vec();
        }
        let error = validate_final_diagnostic_attributes(&zero, Path::new("zero.exr"), 0)
            .expect_err("zero SPP must be refused");
        assert!(error.contains("sample ceiling must be positive"));

        let mut inexact = final_diagnostic();
        for name in ["frankensim.render.spp", "frankensim.render.sppCeiling"] {
            inexact
                .attributes
                .iter_mut()
                .find(|attribute| attribute.name == name)
                .expect("fixture SPP")
                .value = b"16777217".to_vec();
        }
        let error = validate_final_diagnostic_attributes(&inexact, Path::new("inexact.exr"), 0)
            .expect_err("SPP beyond exact FLOAT integer range must be refused");
        assert!(error.contains("exceeds exact FLOAT integer ceiling"));
    }

    #[test]
    fn dimensions_enforce_each_native_4k_axis_not_only_pixel_product() {
        assert!(validate_dimensions(3_840, 2_160, Path::new("4k.exr")).is_ok());
        assert!(validate_dimensions(7_680, 1_080, Path::new("too-wide.exr")).is_err());
        assert!(validate_dimensions(1_920, 4_320, Path::new("too-tall.exr")).is_err());
    }

    #[test]
    fn categorical_identity_planes_must_agree_with_primary_coverage() {
        let mut exr = final_diagnostic();
        exr.channels
            .iter_mut()
            .find(|channel| channel.name == "id.object")
            .expect("fixture object IDs")
            .data[0] = 0.0;
        decode_final_diagnostic(exr, Path::new("fixture.exr"), 0)
            .expect("zero object ID is an admitted unavailable identity on covered raw meshes");

        let mut exr = final_diagnostic();
        exr.channels
            .iter_mut()
            .find(|channel| channel.name == "primary.coverage")
            .expect("fixture coverage")
            .data[0] = 0.0;
        let error = decode_final_diagnostic(exr, Path::new("fixture.exr"), 0)
            .expect_err("background pixels cannot name an object identity");
        assert!(error.contains("id.object sample 0 disagrees with primary coverage"));

        let mut exr = final_diagnostic();
        exr.channels
            .iter_mut()
            .find(|channel| channel.name == "id.material")
            .expect("fixture material IDs")
            .data[1] = 1.5;
        let error = decode_final_diagnostic(exr, Path::new("fixture.exr"), 0)
            .expect_err("fractional categorical identity must be refused");
        assert!(error.contains("is not an exact nonnegative f32 integer palette index"));
    }

    #[test]
    fn cli_requires_a_complete_positive_contiguous_range_request() {
        let error = parse_cli(
            [
                "--input",
                "raw",
                "--output",
                "preview",
                "--frame-start",
                "3",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .expect_err("missing count");
        assert!(error.contains("missing --frame-count"));
        let error = parse_cli(
            [
                "--input",
                "raw",
                "--output",
                "preview",
                "--frame-start",
                "3",
                "--frame-count",
                "0",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .expect_err("zero count");
        assert_eq!(error, "frame-count must be positive");

        let error = parse_cli(
            [
                "--input",
                "raw",
                "--output",
                "preview",
                "--frame-start",
                "3",
                "--frame-count",
                "4",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .expect_err("a mid-sequence cut must be explicit rather than silently temporal");
        assert!(error.contains("requires --initial-cut"));

        assert_eq!(
            parse_cli(
                [
                    "--input",
                    "raw",
                    "--output",
                    "preview",
                    "--frame-start",
                    "3",
                    "--frame-count",
                    "4",
                    "--initial-cut",
                    "--denoise-spatial-passes",
                    "4",
                    "--denoise-spatial-sigma",
                    "0.08",
                ]
                .into_iter()
                .map(str::to_owned),
            )
            .expect("an explicitly cut tuning range"),
            Command::Denoise(Cli {
                input: PathBuf::from("raw"),
                output: PathBuf::from("preview"),
                frame_start: 3,
                frame_count: 4,
                initial_cut: true,
                allow_uniform_spp_transition_at: None,
                denoise_spatial_passes: 4,
                denoise_spatial_sigma: 0.08,
            })
        );

        let transition = parse_cli(
            [
                "--input",
                "raw",
                "--output",
                "preview",
                "--frame-start",
                "0",
                "--frame-count",
                "3",
                "--allow-uniform-spp-transition-at",
                "2",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .expect("one explicit SPP transition inside the requested range");
        assert_eq!(
            transition,
            Command::Denoise(Cli {
                input: PathBuf::from("raw"),
                output: PathBuf::from("preview"),
                frame_start: 0,
                frame_count: 3,
                initial_cut: false,
                allow_uniform_spp_transition_at: Some(2),
                denoise_spatial_passes: TemporalDenoiseConfig::default().spatial_iterations,
                denoise_spatial_sigma: TemporalDenoiseConfig::default().spatial_sigma_rgb,
            })
        );
        for invalid_frame in ["0", "3"] {
            let error = parse_cli(
                [
                    "--input",
                    "raw",
                    "--output",
                    "preview",
                    "--frame-start",
                    "0",
                    "--frame-count",
                    "3",
                    "--allow-uniform-spp-transition-at",
                    invalid_frame,
                ]
                .into_iter()
                .map(str::to_owned),
            )
            .expect_err("transition must be strictly inside the requested range");
            assert!(error.contains("must identify a frame strictly inside"));
        }
    }

    #[test]
    fn cli_admits_one_read_only_sample_inspection_or_the_legacy_denoise_mode() {
        assert_eq!(
            parse_cli(
                ["--inspect-samples", "frame-000096.exr"]
                    .into_iter()
                    .map(str::to_owned)
            )
            .expect("sample inspection command"),
            Command::InspectSamples(PathBuf::from("frame-000096.exr"))
        );
        assert_eq!(
            parse_cli(
                [
                    "--input",
                    "raw",
                    "--output",
                    "preview",
                    "--frame-start",
                    "0",
                    "--frame-count",
                    "2",
                ]
                .into_iter()
                .map(str::to_owned)
            )
            .expect("legacy denoise command"),
            Command::Denoise(Cli {
                input: PathBuf::from("raw"),
                output: PathBuf::from("preview"),
                frame_start: 0,
                frame_count: 2,
                initial_cut: false,
                allow_uniform_spp_transition_at: None,
                denoise_spatial_passes: TemporalDenoiseConfig::default().spatial_iterations,
                denoise_spatial_sigma: TemporalDenoiseConfig::default().spatial_sigma_rgb,
            })
        );
        let error = parse_cli(
            ["--inspect-samples", "frame.exr", "--frame-count", "1"]
                .into_iter()
                .map(str::to_owned),
        )
        .expect_err("inspection must remain one-file read-only mode");
        assert!(error.contains("unexpected argument after --inspect-samples EXR"));
    }

    #[test]
    fn cli_rejects_invalid_spatial_strength_without_relaxing_temporal_guides() {
        for (flag, value, expected) in [
            (
                "--denoise-spatial-passes",
                "0",
                "denoise-spatial-passes must be in",
            ),
            (
                "--denoise-spatial-passes",
                "9",
                "denoise-spatial-passes must be in",
            ),
            (
                "--denoise-spatial-sigma",
                "0",
                "denoise-spatial-sigma must be finite and positive",
            ),
            (
                "--denoise-spatial-sigma",
                "NaN",
                "denoise-spatial-sigma must be finite and positive",
            ),
        ] {
            let error = parse_cli(
                [
                    "--input",
                    "raw",
                    "--output",
                    "preview",
                    "--frame-start",
                    "0",
                    "--frame-count",
                    "1",
                    flag,
                    value,
                ]
                .into_iter()
                .map(str::to_owned),
            )
            .expect_err("invalid spatial strength must refuse");
            assert!(error.contains(expected), "{error}");
        }
    }

    #[test]
    fn every_requested_range_begins_with_one_explicit_history_cut() {
        assert_eq!(denoise_boundary(0), TemporalFrameBoundary::Cut);
        assert_eq!(denoise_boundary(1), TemporalFrameBoundary::Continuous);
        assert_eq!(
            denoise_boundary(usize::MAX),
            TemporalFrameBoundary::Continuous
        );
    }

    #[test]
    fn staging_directory_is_a_sibling_explicitly_marked_incomplete() {
        let output = Path::new("renders/final-preview");
        let staging = staging_output_path(output).expect("staging path");
        assert_eq!(staging.parent(), output.parent());
        let name = staging
            .file_name()
            .and_then(|name| name.to_str())
            .expect("UTF-8 fixture name");
        assert!(name.starts_with(".final-preview.incomplete-"));
        assert_ne!(staging, output);
    }
}
