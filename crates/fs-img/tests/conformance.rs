//! fs-img conformance suite (plan §13.3; the qfx.6 bead). Acceptance:
//! bit-exact deterministic encodes; AOV round-trips lossless; external
//! validation of PNG/EXR outputs (dev-only oracle: macOS `sips`, skipped
//! with a note where absent); the denoiser improves MSE on a fixture
//! while the bias label propagates; fuzzed readers reject structurally.
//! Completed aggregate cases emit canonical fs-obs verdicts. Randomized
//! inputs carry their literal root seeds, while fixed inputs use zero;
//! there is no execution seed in this suite. Assertions and expectations
//! reached before a verdict remain ordinary Rust test diagnostics.

use fs_img::{
    Channel, DenoiseParams, ExrAttribute, LabeledPlane, PixelProvenance, PixelType, PngColor,
    SOURCE_ARTIFACT_HASH_ATTRIBUTE, atrous_denoise, mse, read_exr, read_png, write_exr,
    write_exr_with_attributes, write_png8, write_png16,
};

const SUITE: &str = "fs-img/conformance";
const FIXED_INPUT_SEED: u64 = 0;
const IM_003_INPUT_SEED: u64 = 0x5EED_D401_5E00_0003;
const IM_004_INPUT_SEED: u64 = 0x5EED_F077_0000_0004;

fn print_event(event: fs_obs::Event) {
    fs_obs::lint_failure_record(&event).expect("fs-img event must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("fs-img event must use the fs-obs wire schema");
    println!("{line}");
}

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    print_event(event);
    assert!(pass, "case {case}: {detail}");
}

fn custom_event(scope: &str, severity: fs_obs::Severity, name: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new(SUITE, scope);
    let event = emitter.emit(
        severity,
        fs_obs::EventKind::Custom {
            name: name.to_string(),
            json,
        },
        None,
    );
    print_event(event);
}

fn finite_json_number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6}")
    } else {
        "null".to_string()
    }
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
}

#[test]
fn im_001_encodes_are_bit_exact_and_round_trip() {
    let (w, h) = (16u32, 9u32);
    let px8: Vec<u8> = (0..w * h * 3).map(|i| (i * 31 % 251) as u8).collect();
    let png_a = write_png8(w, h, PngColor::Rgb, &px8).unwrap();
    let png_b = write_png8(w, h, PngColor::Rgb, &px8).unwrap();
    assert_eq!(png_a, png_b, "PNG byte determinism");
    assert_eq!(read_png(&png_a).unwrap().bytes, px8);

    let px16: Vec<u16> = (0..w * h).map(|i| (i * 6151 % 65_521) as u16).collect();
    let png16 = write_png16(w, h, PngColor::Gray, &px16).unwrap();
    assert_eq!(read_png(&png16).unwrap().samples16(), px16);

    let n = (w * h) as usize;
    let chans = vec![
        Channel {
            name: "R".to_string(),
            ty: PixelType::Float,
            data: (0..n).map(|i| i as f32 * 0.5 - 7.0).collect(),
        },
        Channel {
            name: "normal.X".to_string(),
            ty: PixelType::Half,
            data: (0..n).map(|i| (i % 32) as f32 * 0.062_5).collect(),
        },
    ];
    let source_hash = b"0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";
    let attributes = [ExrAttribute {
        name: SOURCE_ARTIFACT_HASH_ATTRIBUTE.to_string(),
        ty: "string".to_string(),
        value: source_hash.to_vec(),
    }];
    let exr_a = write_exr_with_attributes(w, h, &chans, &attributes).unwrap();
    let exr_b = write_exr_with_attributes(w, h, &chans, &attributes).unwrap();
    assert_eq!(exr_a, exr_b, "EXR byte determinism");
    let dec = read_exr(&exr_a).unwrap();
    for c in &dec.channels {
        let orig = chans.iter().find(|o| o.name == c.name).unwrap();
        assert_eq!(c.data, orig.data, "AOV {} round-trip", c.name);
    }
    assert_eq!(dec.attributes.as_slice(), &attributes);
    assert_eq!(
        write_exr(w, h, &chans).unwrap(),
        write_exr_with_attributes(w, h, &chans, &[]).unwrap(),
        "empty metadata must preserve legacy bytes"
    );
    verdict(
        "im-001",
        true,
        "PNG8/PNG16/EXR byte-deterministic; AOV and source-hash metadata round-trip lossless",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn im_002_external_oracle_validates_outputs_when_available() {
    // Dev-only oracle per the bead: macOS `sips` (CoreImage) reads both
    // formats. Skipped with an explicit note when absent (Linux CI).
    if !std::process::Command::new("sips")
        .arg("--help")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
    {
        custom_event(
            "im-002",
            fs_obs::Severity::Warn,
            "image-oracle-skip",
            format!(
                "{{\"oracle\":\"sips\",\"available\":false,\"status\":\"skipped\",\
                 \"reason\":\"sips oracle not present on this machine\",\
                 \"input_seed\":{FIXED_INPUT_SEED}}}"
            ),
        );
        return;
    }
    let dir = std::env::temp_dir();
    let png_path = dir.join(format!("fs-img-oracle-{}.png", std::process::id()));
    let exr_path = dir.join(format!("fs-img-oracle-{}.exr", std::process::id()));
    let px: Vec<u8> = (0..24 * 10 * 3).map(|i| (i % 256) as u8).collect();
    std::fs::write(&png_path, write_png8(24, 10, PngColor::Rgb, &px).unwrap()).unwrap();
    let chan = Channel {
        name: "R".to_string(),
        ty: PixelType::Half,
        data: (0..24 * 10).map(|i| i as f32 / 240.0).collect(),
    };
    std::fs::write(
        &exr_path,
        write_exr(24, 10, std::slice::from_ref(&chan)).unwrap(),
    )
    .unwrap();
    for (path, label) in [(&png_path, "png"), (&exr_path, "exr")] {
        let out = std::process::Command::new("sips")
            .args(["-g", "pixelWidth", "-g", "pixelHeight"])
            .arg(path)
            .output()
            .expect("run sips");
        let text = String::from_utf8_lossy(&out.stdout);
        assert!(
            out.status.success()
                && text.contains("pixelWidth: 24")
                && text.contains("pixelHeight: 10"),
            "sips rejected our {label}: {text}"
        );
    }
    let _ = std::fs::remove_file(&png_path);
    let _ = std::fs::remove_file(&exr_path);
    verdict(
        "im-002",
        true,
        "sips (CoreImage) parsed our PNG and EXR with correct dimensions",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn im_003_denoiser_improves_mse_and_label_propagates() {
    // Fixture: smooth gradient + seeded noise. The denoiser must reduce
    // MSE vs the clean image, and the output must carry the bias tag.
    let (w, h) = (32usize, 32usize);
    let clean: Vec<f32> = (0..w * h)
        .map(|i| f32::midpoint((i % w) as f32 / w as f32, (i / w) as f32 / h as f32))
        .collect();
    let mut seed = IM_003_INPUT_SEED;
    let noisy_data: Vec<f32> = clean
        .iter()
        .map(|&c| c + 0.1 * (lcg(&mut seed) as f32 - 0.5))
        .collect();
    let noisy = LabeledPlane {
        width: w,
        height: h,
        data: noisy_data,
        provenance: PixelProvenance::RawEstimate,
    };
    let out = atrous_denoise(&noisy, None, &DenoiseParams::default()).unwrap();
    let before = mse(&noisy.data, &clean).unwrap();
    let after = mse(&out.data, &clean).unwrap();
    assert!(
        after < before * 0.5,
        "denoiser must clearly improve the fixture: {before:.6} -> {after:.6}"
    );
    assert!(
        matches!(
            out.provenance,
            PixelProvenance::BiasedDenoised { iterations: 3 }
        ),
        "bias label must propagate: {:?}",
        out.provenance
    );
    custom_event(
        "im-003/measurement",
        fs_obs::Severity::Info,
        "image-denoise-mse",
        format!(
            "{{\"before\":{},\"after\":{},\"bias_label\":\"BiasedDenoised\",\
             \"iterations\":3,\"input_seed\":{IM_003_INPUT_SEED}}}",
            finite_json_number(before),
            finite_json_number(after),
        ),
    );
    verdict(
        "im-003",
        true,
        &format!("MSE {before:.5} -> {after:.5}; output labeled biased"),
        IM_003_INPUT_SEED,
    );
}

#[test]
fn im_004_readers_reject_garbage_structurally() {
    let mut seed = IM_004_INPUT_SEED;
    let mut rejected = 0usize;
    for _ in 0..2000 {
        let len = (lcg(&mut seed) * 64.0) as usize;
        let junk: Vec<u8> = (0..len).map(|_| (lcg(&mut seed) * 256.0) as u8).collect();
        if read_png(&junk).is_err() {
            rejected += 1;
        }
        if read_exr(&junk).is_err() {
            rejected += 1;
        }
    }
    assert_eq!(
        rejected, 4000,
        "all 4000 seeded reader attempts must reject random junk"
    );
    // Truncation of a valid file is caught at every prefix length.
    let px = vec![9u8; 27];
    let good = write_png8(3, 3, PngColor::Rgb, &px).unwrap();
    for cut in 1..good.len() {
        assert!(
            read_png(&good[..cut]).is_err(),
            "truncated at {cut} must not decode"
        );
    }
    verdict(
        "im-004",
        true,
        "2000 seeded junk buffers produced 4000 rejected reader attempts; every truncation prefix rejected structurally",
        IM_004_INPUT_SEED,
    );
}
