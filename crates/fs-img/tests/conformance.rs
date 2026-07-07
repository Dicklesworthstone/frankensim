//! fs-img conformance suite (plan §13.3; the qfx.6 bead). Acceptance:
//! bit-exact deterministic encodes; AOV round-trips lossless; external
//! validation of PNG/EXR outputs (dev-only oracle: macOS `sips`, skipped
//! with a note where absent); the denoiser improves MSE on a fixture
//! while the bias label propagates; fuzzed readers reject structurally.

use fs_img::{
    Channel, DenoiseParams, LabeledPlane, PixelProvenance, PixelType, PngColor, atrous_denoise,
    mse, read_exr, read_png, write_exr, write_png8, write_png16,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-img/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
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
    let exr_a = write_exr(w, h, &chans).unwrap();
    let exr_b = write_exr(w, h, &chans).unwrap();
    assert_eq!(exr_a, exr_b, "EXR byte determinism");
    let dec = read_exr(&exr_a).unwrap();
    for c in &dec.channels {
        let orig = chans.iter().find(|o| o.name == c.name).unwrap();
        assert_eq!(c.data, orig.data, "AOV {} round-trip", c.name);
    }
    verdict(
        "im-001",
        "PNG8/PNG16/EXR byte-deterministic; AOV round-trips lossless",
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
        println!(
            "{{\"suite\":\"fs-img/conformance\",\"case\":\"im-002\",\"verdict\":\"skip\",\
             \"detail\":\"sips oracle not present on this machine\"}}"
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
        "sips (CoreImage) parsed our PNG and EXR with correct dimensions",
    );
}

#[test]
fn im_003_denoiser_improves_mse_and_label_propagates() {
    // Fixture: smooth gradient + seeded noise. The denoiser must reduce
    // MSE vs the clean image, and the output must carry the bias tag.
    let (w, h) = (32usize, 32usize);
    let clean: Vec<f32> = (0..w * h)
        .map(|i| ((i % w) as f32 / w as f32 + (i / w) as f32 / h as f32) / 2.0)
        .collect();
    let mut seed = 0x5EED_D401_5E00_0003u64;
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
    println!(
        "{{\"suite\":\"fs-img/conformance\",\"metric\":\"denoise_mse\",\"before\":{before:.6},\
         \"after\":{after:.6},\"bias_label\":\"BiasedDenoised\"}}"
    );
    verdict(
        "im-003",
        &format!("MSE {before:.5} -> {after:.5}; output labeled biased"),
    );
}

#[test]
fn im_004_readers_reject_garbage_structurally() {
    let mut seed = 0x5EED_F077_0000_0004u64;
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
    assert!(
        rejected >= 3999,
        "random junk must essentially never decode: {rejected}/4000"
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
        "4000 junk parses + every truncation prefix rejected structurally",
    );
}
