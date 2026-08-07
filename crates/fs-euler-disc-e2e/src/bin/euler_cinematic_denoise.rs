//! Offline, sequential denoising of Euler-disc DailyCore EXR frames.
//!
//! This command is deliberately a bounded display-derivative producer. It
//! admits only the 14 float planes of the `DailyCore` AOV profile at no more
//! than 3840x2160, keeps only the immediately preceding biased denoise frame,
//! and never relabels its PNGs as raw render estimates.

use std::{
    collections::BTreeMap,
    fs::{self, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
};

use fs_img::{
    Channel, CinematicColorConfig, CinematicColorLimits, DecodedExr, ExrAttribute,
    ExrInspectLimits, PixelType, PngColor, PreviewDither, TemporalDenoiseConfig,
    TemporalDenoiseInput, TemporalDenoiseLimits, TemporalDenoisedFrame, TemporalFrameBoundary,
    inspect_exr, read_exr, temporal_denoise_rgb, transform_cinematic_preview, write_png16,
};

const DAILY_CORE_CHANNELS: [&str; 14] = [
    "B",
    "G",
    "R",
    "albedo.B",
    "albedo.G",
    "albedo.R",
    "depth.Z",
    "motion.prev.X",
    "motion.prev.Y",
    "normal.X",
    "normal.Y",
    "normal.Z",
    "primary.coverage",
    "variance.Y",
];
const MAX_4K_PIXELS: u64 = 3_840 * 2_160;
const MAX_DAILY_CORE_DECODED_BYTES: u64 = MAX_4K_PIXELS * DAILY_CORE_CHANNELS.len() as u64 * 4;
const MAX_DAILY_CORE_ENCODED_BYTES: u64 = 512 * 1024 * 1024;
const MAX_EXR_HEADER_BYTES: u64 = 1024 * 1024;
const MAX_EXR_METADATA_BYTES: u64 = 1024 * 1024;
const DAILY_CORE_PROFILE: &str = "daily-core-v1";

#[derive(Debug, Clone, PartialEq, Eq)]
struct Cli {
    input: PathBuf,
    output: PathBuf,
    frame_start: u64,
    frame_count: u64,
}

#[derive(Debug)]
struct DailyCoreFrame {
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
    sequence_identity: SequenceIdentity,
}

/// Header values which must name one coherent raw-render sequence before the
/// history-dependent denoiser is allowed to join frames.
#[derive(Debug, Clone, PartialEq, Eq)]
struct SequenceIdentity {
    source_trajectory: String,
    scene_hash: String,
    composition: String,
    aov_profile: String,
}

impl DailyCoreFrame {
    fn temporal_input(&self, frame_index: u64) -> TemporalDenoiseInput<'_> {
        TemporalDenoiseInput {
            frame_index,
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
            object_ids: None,
            material_ids: None,
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
    let cli = parse_cli(std::env::args().skip(1))?;
    run_cli(&cli, |message| eprintln!("status=progress {message}"))
}

fn run_cli(cli: &Cli, mut progress: impl FnMut(&str)) -> Result<(), String> {
    if !cli.input.is_dir() {
        return Err(format!("input is not a directory: {}", cli.input.display()));
    }
    let frame_end = cli
        .frame_start
        .checked_add(cli.frame_count)
        .ok_or_else(|| "frame-start plus frame-count overflows u64".to_owned())?;
    if cli.frame_count == 0 {
        return Err("frame-count must be positive".to_owned());
    }
    fs::create_dir_all(&cli.output)
        .map_err(|error| format!("create output directory {}: {error}", cli.output.display()))?;
    if !cli.output.is_dir() {
        return Err(format!(
            "output is not a directory: {}",
            cli.output.display()
        ));
    }

    // Refuse an entire range before computing: a rerun can never overwrite a
    // prior preview, and no later frame is rendered after this preflight fails.
    for frame in cli.frame_start..frame_end {
        let output = preview_path(&cli.output, frame);
        match fs::metadata(&output) {
            Ok(_) => {
                return Err(format!(
                    "refusing to overwrite existing output: {}",
                    output.display()
                ));
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(format!("inspect output {}: {error}", output.display())),
        }
    }

    progress(&format!(
        "stage=denoise begin input={} output={} frame_start={} frame_count={}",
        cli.input.display(),
        cli.output.display(),
        cli.frame_start,
        cli.frame_count,
    ));
    let mut color = CinematicColorConfig::reference_srgb_16();
    color.exposure_ev = 1;
    color.dither = PreviewDither::Disabled;
    let denoise_config = TemporalDenoiseConfig::default();
    let mut history: Option<TemporalDenoisedFrame> = None;
    let mut expected_dimensions: Option<(u32, u32)> = None;
    let mut expected_sequence: Option<SequenceIdentity> = None;

    for (ordinal, frame) in (cli.frame_start..frame_end).enumerate() {
        let input_path = raw_path(&cli.input, frame);
        progress(&format!(
            "stage=denoise frame={}/{} absolute_frame={} action=read path={}",
            ordinal + 1,
            cli.frame_count,
            frame,
            input_path.display(),
        ));
        let raw = read_daily_core(&input_path, frame)?;
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
        if let Some(expected) = &expected_sequence {
            if raw.sequence_identity != *expected {
                return Err(format!(
                    "frame continuity violation at absolute frame {frame}: EXR provenance identity differs from the first requested frame"
                ));
            }
        } else {
            expected_sequence = Some(raw.sequence_identity.clone());
        }

        let boundary = if ordinal == 0 {
            TemporalFrameBoundary::Cut
        } else {
            TemporalFrameBoundary::Continuous
        };
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
        let output_path = preview_path(&cli.output, frame);
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
    progress(&format!(
        "stage=denoise complete frame_start={} frame_count={}",
        cli.frame_start, cli.frame_count
    ));
    Ok(())
}

fn read_daily_core(path: &Path, expected_frame: u64) -> Result<DailyCoreFrame, String> {
    let length = fs::metadata(path)
        .map_err(|error| format!("inspect input {}: {error}", path.display()))?
        .len();
    if length > MAX_DAILY_CORE_ENCODED_BYTES {
        return Err(format!(
            "DailyCore EXR {} is {length} bytes; maximum encoded input is {MAX_DAILY_CORE_ENCODED_BYTES} bytes",
            path.display()
        ));
    }
    let bytes =
        fs::read(path).map_err(|error| format!("read input {}: {error}", path.display()))?;
    let inspection = inspect_exr(
        &bytes,
        ExrInspectLimits {
            max_input_bytes: MAX_DAILY_CORE_ENCODED_BYTES,
            max_header_bytes: MAX_EXR_HEADER_BYTES,
            max_decoded_bytes: MAX_DAILY_CORE_DECODED_BYTES,
            max_metadata_bytes: MAX_EXR_METADATA_BYTES,
        },
    )
    .map_err(|error| format!("inspect DailyCore EXR {}: {error}", path.display()))?;
    validate_dimensions(inspection.width, inspection.height, path)?;
    let decoded = read_exr(&bytes)
        .map_err(|error| format!("decode DailyCore EXR {}: {error}", path.display()))?;
    drop(bytes);
    decode_daily_core(decoded, path, expected_frame)
}

fn decode_daily_core(
    exr: DecodedExr,
    path: &Path,
    expected_frame: u64,
) -> Result<DailyCoreFrame, String> {
    validate_dimensions(exr.width, exr.height, path)?;
    let sequence_identity = validate_daily_core_attributes(&exr, path, expected_frame)?;
    let pixels = usize::try_from(u64::from(exr.width) * u64::from(exr.height)).map_err(|_| {
        format!(
            "DailyCore EXR {} pixel count does not fit usize",
            path.display()
        )
    })?;
    let mut planes = BTreeMap::new();
    for Channel { name, ty, data } in exr.channels {
        if ty != PixelType::Float {
            return Err(format!(
                "DailyCore EXR {} channel {name} is not FLOAT",
                path.display()
            ));
        }
        if data.len() != pixels {
            return Err(format!(
                "DailyCore EXR {} channel {name} has {} samples; expected {pixels}",
                path.display(),
                data.len()
            ));
        }
        if planes.insert(name.clone(), data).is_some() {
            return Err(format!(
                "DailyCore EXR {} contains duplicate channel {name}",
                path.display()
            ));
        }
    }
    for name in DAILY_CORE_CHANNELS {
        if !planes.contains_key(name) {
            return Err(format!(
                "DailyCore EXR {} is missing required channel {name}",
                path.display()
            ));
        }
    }
    if let Some((unexpected, _)) = planes
        .iter()
        .find(|(name, _)| !DAILY_CORE_CHANNELS.contains(&name.as_str()))
    {
        return Err(format!(
            "DailyCore EXR {} contains unexpected channel {unexpected}",
            path.display()
        ));
    }
    Ok(DailyCoreFrame {
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
        primary_coverage: take_plane(&mut planes, "primary.coverage", path)?,
        variance_luminance: take_plane(&mut planes, "variance.Y", path)?,
        sequence_identity,
    })
}

fn take_plane(
    planes: &mut BTreeMap<String, Vec<f32>>,
    name: &str,
    path: &Path,
) -> Result<Vec<f32>, String> {
    planes.remove(name).ok_or_else(|| {
        format!(
            "DailyCore EXR {} lost validated required channel {name} during reconstruction",
            path.display()
        )
    })
}

fn validate_daily_core_attributes(
    exr: &DecodedExr,
    path: &Path,
    expected_frame: u64,
) -> Result<SequenceIdentity, String> {
    let frame_index = exr_attribute_u64(exr, "frankensim.frame.index", path)?;
    if frame_index != expected_frame {
        return Err(format!(
            "DailyCore EXR {} declares frame index {frame_index}; expected {expected_frame}",
            path.display()
        ));
    }
    let aov_profile = exr_attribute_string(exr, "frankensim.aov.profile", path)?;
    if aov_profile != DAILY_CORE_PROFILE {
        return Err(format!(
            "DailyCore EXR {} declares AOV profile {aov_profile:?}; expected {DAILY_CORE_PROFILE:?}",
            path.display()
        ));
    }
    // The per-frame AOV configuration hash is mandatory integrity metadata,
    // but it is not a sequence identity: the hash deliberately commits to the
    // absolute frame index and neighbouring frame times.  Requiring equality
    // across frames would reject every valid temporal sequence at frame 1.
    let _frame_config_hash = exr_attribute_string(exr, "frankensim.aov.configHash", path)?;
    Ok(SequenceIdentity {
        source_trajectory: exr_attribute_string(exr, "frankensim.source.trajectory", path)?,
        scene_hash: exr_attribute_string(exr, "frankensim.source.sceneHash", path)?,
        composition: exr_attribute_string(exr, "frankensim.source.composition", path)?,
        aov_profile,
    })
}

fn exr_attribute_string(exr: &DecodedExr, name: &str, path: &Path) -> Result<String, String> {
    let attribute = required_exr_attribute(exr, name, path)?;
    if attribute.ty != "string" {
        return Err(format!(
            "DailyCore EXR {} attribute {name} has type {:?}; expected string",
            path.display(),
            attribute.ty
        ));
    }
    String::from_utf8(attribute.value.clone()).map_err(|_| {
        format!(
            "DailyCore EXR {} attribute {name} is not valid UTF-8 string data",
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
                    "DailyCore EXR {} attribute {name} is not valid UTF-8 string data",
                    path.display()
                )
            })?
            .parse::<u64>()
            .map_err(|_| {
                format!(
                    "DailyCore EXR {} attribute {name} is not a decimal u64",
                    path.display()
                )
            }),
        "uint64" => {
            let bytes: [u8; 8] = attribute.value.as_slice().try_into().map_err(|_| {
                format!(
                    "DailyCore EXR {} attribute {name} uint64 payload is not exactly eight bytes",
                    path.display()
                )
            })?;
            Ok(u64::from_le_bytes(bytes))
        }
        ty => Err(format!(
            "DailyCore EXR {} attribute {name} has type {ty:?}; expected string or uint64",
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
            "DailyCore EXR {} is missing required attribute {name}",
            path.display()
        )
    })?;
    if found.next().is_some() {
        return Err(format!(
            "DailyCore EXR {} has duplicate required attribute {name}",
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
                "DailyCore EXR {} dimension product overflows",
                path.display()
            )
        })?;
    if width == 0 || height == 0 || pixels > MAX_4K_PIXELS {
        return Err(format!(
            "DailyCore EXR {} dimensions {}x{} exceed the bounded 4K envelope",
            path.display(),
            width,
            height,
        ));
    }
    Ok(())
}

fn raw_path(input: &Path, frame: u64) -> PathBuf {
    input.join(format!("frame-{frame:06}.exr"))
}

fn preview_path(output: &Path, frame: u64) -> PathBuf {
    output.join(format!("frame-{frame:06}.png"))
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

fn parse_cli(args: impl Iterator<Item = String>) -> Result<Cli, String> {
    let mut input = None;
    let mut output = None;
    let mut frame_start = None;
    let mut frame_count = None;
    let mut args = args;
    while let Some(argument) = args.next() {
        match argument.as_str() {
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
            "--help" | "-h" => return Err(usage().to_owned()),
            _ => return Err(format!("unknown argument: {argument}\n{}", usage())),
        }
    }
    let cli = Cli {
        input: input.ok_or_else(|| format!("missing --input\n{}", usage()))?,
        output: output.ok_or_else(|| format!("missing --output\n{}", usage()))?,
        frame_start: frame_start.ok_or_else(|| format!("missing --frame-start\n{}", usage()))?,
        frame_count: frame_count.ok_or_else(|| format!("missing --frame-count\n{}", usage()))?,
    };
    if cli.frame_count == 0 {
        return Err("frame-count must be positive".to_owned());
    }
    Ok(cli)
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

const fn usage() -> &'static str {
    "Usage: euler_cinematic_denoise --input DIR --output DIR --frame-start N --frame-count N"
}

#[cfg(test)]
mod tests {
    use super::*;

    fn daily_core() -> DecodedExr {
        let channels = DAILY_CORE_CHANNELS
            .iter()
            .map(|name| Channel {
                name: (*name).to_owned(),
                ty: PixelType::Float,
                data: vec![1.0, 2.0],
            })
            .collect();
        DecodedExr {
            width: 2,
            height: 1,
            channels,
            attributes: [
                ("frankensim.frame.index", "0"),
                ("frankensim.source.trajectory", "trajectory-a"),
                ("frankensim.source.sceneHash", "scene-a"),
                ("frankensim.source.composition", "composition-a"),
                ("frankensim.aov.profile", DAILY_CORE_PROFILE),
                ("frankensim.aov.configHash", "config-a"),
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
    fn daily_core_schema_reconstructs_exact_denoise_planes() {
        let frame =
            decode_daily_core(daily_core(), Path::new("fixture.exr"), 0).expect("DailyCore");
        assert_eq!(frame.temporal_input(42).frame_index, 42);
        assert_eq!(frame.temporal_input(42).red, &[1.0, 2.0]);
        assert_eq!(frame.temporal_input(42).variance_luminance, &[1.0, 2.0]);
    }

    #[test]
    fn missing_daily_core_plane_refuses_before_denoising() {
        let mut exr = daily_core();
        exr.channels
            .retain(|channel| channel.name != "motion.prev.Y");
        let error = decode_daily_core(exr, Path::new("fixture.exr"), 0).expect_err("missing plane");
        assert!(error.contains("missing required channel motion.prev.Y"));
    }

    #[test]
    fn header_frame_index_must_match_the_absolute_requested_frame() {
        let error = decode_daily_core(daily_core(), Path::new("fixture.exr"), 1)
            .expect_err("mismatched metadata must refuse temporal history");
        assert!(error.contains("declares frame index 0; expected 1"));
    }

    #[test]
    fn uint64_frame_index_attribute_is_accepted_when_it_matches() {
        let mut exr = daily_core();
        let attribute = exr
            .attributes
            .iter_mut()
            .find(|attribute| attribute.name == "frankensim.frame.index")
            .expect("fixture frame index");
        attribute.ty = "uint64".to_owned();
        attribute.value = 7_u64.to_le_bytes().to_vec();
        decode_daily_core(exr, Path::new("fixture.exr"), 7)
            .expect("matching uint64 metadata must be accepted");
    }

    #[test]
    fn consecutive_frames_allow_distinct_frame_bound_aov_config_hashes() {
        let first = daily_core();
        let first_identity = validate_daily_core_attributes(&first, Path::new("frame-0.exr"), 0)
            .expect("first frame provenance");

        let mut second = daily_core();
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
        config_hash.value = b"config-frame-1".to_vec();

        let second_identity = validate_daily_core_attributes(&second, Path::new("frame-1.exr"), 1)
            .expect("second frame provenance");
        assert_eq!(first_identity, second_identity);
    }

    #[test]
    fn every_frame_still_requires_its_aov_config_hash() {
        let mut exr = daily_core();
        exr.attributes
            .retain(|attribute| attribute.name != "frankensim.aov.configHash");
        let error = validate_daily_core_attributes(&exr, Path::new("frame-0.exr"), 0)
            .expect_err("missing per-frame integrity metadata must be refused");
        assert!(error.contains("missing required attribute frankensim.aov.configHash"));
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
    }
}
