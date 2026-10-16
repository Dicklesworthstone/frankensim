//! E2E for the music-lane render CLI (bead
//! `frankensim-music-t-out-render-ib15w`): drive the REAL binary; assert
//! determinism (same args → bit-identical WAV + provenance), the
//! never-overwrite refusal, the fixture/argument refusals, RIFF shape,
//! provenance completeness, and the never-peak-normalize law (a hotter
//! full-scale yields QUIETER samples of the SAME physics, and clipping is
//! counted, not hidden).

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_music_render")
}

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("music-render-{}-{name}", std::process::id()));
    if dir.exists() {
        std::fs::remove_dir_all(&dir).expect("clean scratch");
    }
    std::fs::create_dir_all(&dir).expect("mkdir scratch");
    dir
}

fn run(args: &[&str]) -> (bool, String) {
    let output = Command::new(bin()).args(args).output().expect("spawn");
    (
        output.status.success(),
        String::from_utf8_lossy(&output.stdout).into_owned(),
    )
}

#[test]
fn physical_inputs_reach_real_cli_audio_and_provenance() {
    // Keep artifacts; unique creation avoids overwriting an earlier run.
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("physical-inputs-{}-{stamp}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let render = |name: &str, extra: &[&str]| {
        let path = dir.join(format!("{name}.wav"));
        let mut args = vec![
            "reed",
            path.to_str().unwrap(),
            "--seconds",
            "0.1",
            "--full-scale-pa",
            "200000",
        ];
        args.extend_from_slice(extra);
        let (ok, stdout) = run(&args);
        assert!(ok, "{stdout}");
        (
            std::fs::read(&path).unwrap(),
            std::fs::read_to_string(path.with_extension("provenance.json")).unwrap(),
        )
    };
    let (baseline, _) = render("baseline", &[]);
    let (explicit_dry, dry_provenance) = render("dry", &["--relative-humidity", "0"]);
    assert_eq!(baseline, explicit_dry);
    assert!(dry_provenance.contains("\"spec\":\"dry_air_ussa1976\""));
    assert!(dry_provenance.contains("\"relative_humidity\":0e0"));
    let (humid, humid_provenance) =
        render("humid-37", &["--relative-humidity", "0.8", "--block", "37"]);
    let (humid_replay, _) = render(
        "humid-512",
        &["--relative-humidity", "0.8", "--block", "512"],
    );
    assert_ne!(humid, baseline);
    assert_eq!(humid, humid_replay);
    assert!(humid_provenance.contains("\"spec\":\"moist_air_ussa1976_water_vapor_nist\""));
    assert!(humid_provenance.contains("\"relative_humidity\":8e-1"));
    // Interrupt an unfinished ramp, then ramp from the value actually reached.
    // Exercise decoding, control compilation and the retained voice together.
    use fs_scenario::gesture::{GestureEvent, GestureSchedule, GestureValue};
    let mut track = pressure_performance().tracks()[0].clone();
    track.initial = GestureValue::PressurePa(0.0);
    track.events = vec![
        GestureEvent {
            time_s: 0.0,
            transition_s: 0.08,
            value: GestureValue::PressurePa(4000.0),
        },
        GestureEvent {
            time_s: 0.02,
            transition_s: 0.02,
            value: GestureValue::PressurePa(2000.0),
        },
    ];
    let source = dir.join("interrupted.gesture");
    let schedule = GestureSchedule::try_new(137, vec![track.clone()]).unwrap();
    std::fs::write(&source, schedule.to_canonical_bytes()).unwrap();
    let (performed, _) = render(
        "interrupted-37",
        &["--schedule", source.to_str().unwrap(), "--block", "37"],
    );
    let (replayed, _) = render(
        "interrupted-512",
        &["--schedule", source.to_str().unwrap(), "--block", "512"],
    );
    assert_eq!(performed, replayed);
    assert_ne!(performed, baseline);
    // The same physical trajectory with its first ramp explicitly shortened:
    // 4000 Pa * (0.02 / 0.08) = 1000 Pa at interruption.
    track.events[0].transition_s = 0.02;
    track.events[0].value = GestureValue::PressurePa(1000.0);
    let equivalent = GestureSchedule::try_new(137, vec![track]).unwrap();
    let equivalent_source = dir.join("equivalent.gesture");
    std::fs::write(&equivalent_source, equivalent.to_canonical_bytes()).unwrap();
    let (equivalent_audio, _) = render(
        "equivalent",
        &["--schedule", equivalent_source.to_str().unwrap()],
    );
    assert_eq!(performed, equivalent_audio);
    let options = [
        "--temperature-k",
        "313.15",
        "--ambient-pressure-pa",
        "80000",
        "--relative-humidity",
        "0.5",
        "--duct-length-m",
        "0.6",
        "--duct-radius-m",
        "0.002",
        "--duct-outlet-radius-m",
        "0.0015",
    ];
    let mut a_args = options.to_vec();
    a_args.extend_from_slice(&["--block", "37"]);
    let (a, provenance) = render("a", &a_args);
    let mut b_args = options.to_vec();
    b_args.extend_from_slice(&["--block", "512"]);
    let (b, _) = render("b", &b_args);
    assert_eq!(a, b);
    assert_ne!(a, baseline);
    assert_eq!(&a[..4], b"RIFF");
    for field in [
        "\"temperature_k\":3.1315e2",
        "\"pressure_pa\":8e4",
        "\"relative_humidity\":5e-1",
        "\"spec\":\"moist_air_ussa1976_water_vapor_nist\"",
        "\"shape\":\"cone\"",
        "\"length_m\":6e-1",
        "\"inlet_radius_m\":2e-3",
        "\"outlet_radius_m\":1.5e-3",
    ] {
        assert!(provenance.contains(field), "missing {field}: {provenance}");
    }
    for (index, (fixture, option, value)) in [
        ("reed", "--temperature-k", "NaN"),
        ("reed", "--ambient-pressure-pa", "0"),
        ("reed", "--relative-humidity", "NaN"),
        ("reed", "--relative-humidity", "inf"),
        ("reed", "--relative-humidity", "-0.1"),
        ("reed", "--relative-humidity", "1.1"),
        ("reed", "--duct-outlet-radius-m", "-1"),
        ("string", "--duct-outlet-radius-m", "0.002"),
        ("string", "--temperature-k", "313.15"),
        ("string", "--relative-humidity", "0.5"),
    ]
    .into_iter()
    .enumerate()
    {
        let path = dir.join(format!("refused-{index}.wav"));
        let (ok, stdout) = run(&[fixture, path.to_str().unwrap(), option, value]);
        assert!(!ok, "{stdout}");
        assert!(stdout.contains("\"verdict\":\"refused\""));
        assert!(!path.exists());
        assert!(!path.with_extension("provenance.json").exists());
    }
    // Mixture-owner limits must refuse before opening output files, even
    // where the dry-gas constructor alone would admit the temperature/pressure.
    for (name, extra) in [
        ("hot-moist", ["--temperature-k", "400"]),
        ("vapor-rich", ["--ambient-pressure-pa", "1000"]),
    ] {
        let path = dir.join(format!("{name}.wav"));
        let (ok, stdout) = run(&[
            "reed",
            path.to_str().unwrap(),
            "--relative-humidity",
            "1",
            extra[0],
            extra[1],
        ]);
        assert!(!ok, "{stdout}");
        assert!(stdout.contains("ambient gas state refused"), "{stdout}");
        assert!(!path.exists());
        assert!(!path.with_extension("provenance.json").exists());
    }
    println!("retained CLI artifacts: {}", dir.display());
}

#[test]
fn renders_deterministically_with_provenance() {
    let dir = scratch("determinism");
    let a = dir.join("a.wav");
    let b = dir.join("b.wav");
    for path in [&a, &b] {
        let (ok, stdout) = run(&["string", path.to_str().expect("utf8"), "--seconds", "0.25"]);
        assert!(ok, "render must succeed:\n{stdout}");
        assert!(stdout.contains("\"verdict\":\"rendered\""), "{stdout}");
        assert!(stdout.contains("\"wav_blake3\":\""), "{stdout}");
    }
    let wav_a = std::fs::read(&a).expect("read a");
    let wav_b = std::fs::read(&b).expect("read b");
    assert_eq!(wav_a, wav_b, "same args must produce bit-identical WAVs");
    assert_eq!(&wav_a[0..4], b"RIFF", "WAV container shape");
    let prov_a = std::fs::read_to_string(a.with_extension("provenance.json")).expect("sidecar");
    let prov_b = std::fs::read_to_string(b.with_extension("provenance.json")).expect("sidecar");
    // Sidecars differ only in nothing — fully deterministic.
    assert_eq!(prov_a, prov_b);
    for field in [
        "\"schema\":\"frankensim-music-render-provenance-v1\"",
        "\"sample_rate_hz\":48000",
        "\"clipped_samples\":",
        "\"wav_blake3\":\"",
        "never peak-normalized",
    ] {
        assert!(
            prov_a.contains(field),
            "sidecar missing {field:?}:\n{prov_a}"
        );
    }
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn refuses_to_overwrite_evidence() {
    let dir = scratch("overwrite");
    let out = dir.join("once.wav");
    let (ok, _) = run(&["string", out.to_str().expect("utf8"), "--seconds", "0.05"]);
    assert!(ok);
    let before = std::fs::read(&out).expect("read");
    let (ok, stdout) = run(&["string", out.to_str().expect("utf8"), "--seconds", "0.05"]);
    assert!(!ok, "re-render onto an existing path must refuse");
    assert!(stdout.contains("refuses to overwrite evidence"), "{stdout}");
    assert_eq!(
        std::fs::read(&out).expect("read"),
        before,
        "the refusal must leave the artifact untouched"
    );
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn refusal_arms_are_typed() {
    let dir = scratch("refusals");
    let out = dir.join("x.wav");
    let out_str = out.to_str().expect("utf8");
    for (args, needle) in [
        (vec!["kazoo", out_str], "fixture must be"),
        (vec!["string"], "usage:"),
        (
            vec!["string", out_str, "--seconds", "0"],
            "--seconds must be",
        ),
        (vec!["string", out_str, "--block", "0"], "--block must be"),
    ] {
        let (ok, stdout) = run(&args);
        assert!(!ok, "args {args:?} must refuse");
        assert!(
            stdout.contains(needle),
            "args {args:?}: wrong refusal:\n{stdout}"
        );
    }
    assert!(!out.exists(), "refusals must write nothing");
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn full_scale_is_physics_not_normalization() {
    // The never-peak-normalize law, observable: rendering the same
    // physics at 10x the full-scale (200 kPa vs 20 kPa: the reed peaks ~10.7 kPa, so both scales are clip-free) must produce the SAME peak_pa in
    // provenance (physics unchanged) and quieter PCM (samples scale
    // down); a tiny full-scale must CLIP and say so.
    let dir = scratch("fullscale");
    let quiet = dir.join("quiet.wav");
    let loud = dir.join("loud.wav");
    let clipped = dir.join("clipped.wav");
    let (ok, _) = run(&[
        "reed",
        quiet.to_str().expect("utf8"),
        "--seconds",
        "0.1",
        "--full-scale-pa",
        "200000",
    ]);
    assert!(ok);
    let (ok, _) = run(&[
        "reed",
        loud.to_str().expect("utf8"),
        "--seconds",
        "0.1",
        "--full-scale-pa",
        "20000",
    ]);
    assert!(ok);
    let (ok, stdout) = run(&[
        "reed",
        clipped.to_str().expect("utf8"),
        "--seconds",
        "0.1",
        "--full-scale-pa",
        "10",
    ]);
    assert!(ok, "clipping renders (counted, not refused):\n{stdout}");

    let read_prov =
        |p: &Path| std::fs::read_to_string(p.with_extension("provenance.json")).expect("sidecar");
    let peak_of = |prov: &str| -> f64 {
        let start = prov.find("\"peak_pa\":").expect("peak field") + "\"peak_pa\":".len();
        let rest = &prov[start..];
        let end = rest.find(',').expect("comma");
        rest[..end].parse().expect("peak parse")
    };
    let clips_of = |prov: &str| -> u64 {
        let start =
            prov.find("\"clipped_samples\":").expect("clip field") + "\"clipped_samples\":".len();
        let rest = &prov[start..];
        let end = rest.find(',').expect("comma");
        rest[..end].parse().expect("clip parse")
    };
    let quiet_prov = read_prov(&quiet);
    let loud_prov = read_prov(&loud);
    let clipped_prov = read_prov(&clipped);
    // Same physics: identical peak pascals across full-scale choices.
    assert!(
        (peak_of(&quiet_prov) - peak_of(&loud_prov)).abs() < 1.0e-12,
        "full-scale must not change the physics"
    );
    assert_eq!(clips_of(&quiet_prov), 0, "200 kPa full-scale must not clip");
    assert_eq!(
        clips_of(&loud_prov),
        0,
        "20 kPa full-scale must not clip either"
    );
    assert!(
        clips_of(&clipped_prov) > 0,
        "10 Pa full-scale must clip a screaming reed and SAY so"
    );
    // Quieter mapping: the 200-kPa WAV's peak sample magnitude is ~10x
    // smaller than the 20-kPa one's (both clip-free, so the ratio is exact).
    let peak_sample = |path: &Path| -> i32 {
        let bytes = std::fs::read(path).expect("wav");
        bytes[44..]
            .as_chunks::<2>()
            .0
            .iter()
            .map(|c| i32::from(i16::from_le_bytes(*c)).abs())
            .max()
            .unwrap_or(0)
    };
    let ratio = f64::from(peak_sample(&loud)) / f64::from(peak_sample(&quiet)).max(1.0);
    assert!(
        (8.0..=12.0).contains(&ratio),
        "PCM peaks must scale ~10x with a 10x full-scale change (got {ratio})"
    );
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

fn pressure_performance() -> fs_scenario::gesture::GestureSchedule {
    use fs_scenario::gesture::{
        GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue,
    };
    GestureSchedule::try_new(
        137,
        vec![GestureTrack {
            id: "blow\"quoted".to_string(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(2800.0),
            events: vec![
                GestureEvent {
                    time_s: 0.02,
                    transition_s: 0.01,
                    value: GestureValue::PressurePa(0.0),
                },
                GestureEvent {
                    time_s: 0.06,
                    transition_s: 0.01,
                    value: GestureValue::PressurePa(3500.0),
                },
            ],
        }],
    )
    .expect("performance admits")
}

#[test]
fn gesture_file_changes_real_physics_and_replays_across_audio_blocks() {
    let dir = scratch("gesture-replay");
    let schedule = pressure_performance();
    let first_source = dir.join("first.gesture");
    let moved_source = dir.join("relocated.gesture");
    let source_bytes = schedule.to_canonical_bytes();
    std::fs::write(&first_source, &source_bytes).unwrap();
    std::fs::write(&moved_source, &source_bytes).unwrap();
    let a = dir.join("a.wav");
    let b = dir.join("b.wav");
    let replay = dir.join("replay.wav");
    for (output, source, block) in [
        (&a, &first_source, "37"),
        (&b, &first_source, "512"),
        (&replay, &moved_source, "37"),
    ] {
        let (ok, stdout) = run(&[
            "reed",
            output.to_str().unwrap(),
            "--seconds",
            "0.1",
            "--full-scale-pa",
            "200000",
            "--block",
            block,
            "--schedule",
            source.to_str().unwrap(),
        ]);
        assert!(ok, "scheduled render must succeed: {stdout}");
        assert!(stdout.contains("\"verdict\":\"rendered\""));
    }
    let audio = std::fs::read(&a).unwrap();
    assert_eq!(
        audio,
        std::fs::read(&b).unwrap(),
        "host block size changed the performed audio"
    );
    assert_eq!(
        audio,
        std::fs::read(&replay).unwrap(),
        "source relocation changed the performed audio"
    );
    let provenance = std::fs::read_to_string(a.with_extension("provenance.json")).unwrap();
    assert_eq!(
        provenance,
        std::fs::read_to_string(replay.with_extension("provenance.json")).unwrap()
    );
    #[allow(clippy::format_collect)]
    let source_hash: String = schedule
        .content_hash()
        .0
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    assert!(provenance.contains(&format!("\"blake3\":\"{source_hash}\"")));
    assert!(provenance.contains("\"control_rate_hz\":137"));
    assert!(provenance.contains("\"track\":\"blow\\\"quoted\""));
    assert!(provenance.contains("ceil-control-tick-to-audio-sample"));

    // A vacuous file-loader test could produce the old fixed fixture unchanged.
    // Require the same real voice WITHOUT the authored release/re-entry to differ.
    let held = dir.join("held.wav");
    let (ok, stdout) = run(&[
        "reed",
        held.to_str().unwrap(),
        "--seconds",
        "0.1",
        "--full-scale-pa",
        "200000",
        "--block",
        "37",
    ]);
    assert!(ok, "unscheduled reference failed: {stdout}");
    assert_ne!(
        audio,
        std::fs::read(&held).unwrap(),
        "gesture input had no physical effect"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn malformed_or_unsupported_schedules_refuse_before_writing_artifacts() {
    use fs_scenario::gesture::{
        GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue,
    };
    let dir = scratch("gesture-refusals");
    let valid = pressure_performance().to_canonical_bytes();
    let mut trailing = valid.clone();
    trailing.extend_from_slice(b"ignored\n");
    let unsupported = GestureSchedule::try_new(
        137,
        vec![GestureTrack {
            id: "pedal".to_string(),
            target: GestureTarget::SustainPedal,
            initial: GestureValue::Fraction(0.0),
            events: Vec::new(),
        }],
    )
    .unwrap()
    .to_canonical_bytes();
    let mut overlapping_track = pressure_performance().tracks()[0].clone();
    overlapping_track.events = vec![
        GestureEvent {
            time_s: 0.0,
            transition_s: 0.1,
            value: GestureValue::PressurePa(3000.0),
        },
        GestureEvent {
            time_s: 0.01,
            transition_s: 0.0,
            value: GestureValue::PressurePa(0.0),
        },
    ];
    let overlapping = GestureSchedule::try_new(137, vec![overlapping_track])
        .unwrap()
        .to_canonical_bytes();
    let excessive =
        b"frankensim-gesture-schedule-v1\ncontrol_rate_hz\t137\ntracks\t18446744073709551615\n"
            .to_vec();
    for (index, (bytes, needle)) in [
        (trailing, "canonical"),
        (unsupported, "not a blowing-pressure input"),
        (excessive, "count"),
    ]
    .into_iter()
    .enumerate()
    {
        let source = dir.join(format!("bad-{index}.gesture"));
        let output = dir.join(format!("bad-{index}.wav"));
        std::fs::write(&source, bytes).unwrap();
        let (ok, stdout) = run(&[
            "reed",
            output.to_str().unwrap(),
            "--schedule",
            source.to_str().unwrap(),
        ]);
        assert!(!ok, "invalid schedule rendered: {stdout}");
        assert!(stdout.contains(needle), "wrong refusal: {stdout}");
        assert_eq!(
            stdout.lines().count(),
            1,
            "refusal must be one escaped JSON record"
        );
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
    let source = dir.join("interrupted.gesture");
    let output = dir.join("interrupted.wav");
    std::fs::write(&source, overlapping).unwrap();
    let (ok, stdout) = run(&[
        "reed",
        output.to_str().unwrap(),
        "--seconds",
        "0.1",
        "--schedule",
        source.to_str().unwrap(),
    ]);
    assert!(ok, "supported interrupted ramp refused: {stdout}");
    let source = dir.join("valid.gesture");
    std::fs::write(&source, valid).unwrap();
    let output = dir.join("string.wav");
    let (ok, stdout) = run(&[
        "string",
        output.to_str().unwrap(),
        "--schedule",
        source.to_str().unwrap(),
    ]);
    assert!(!ok);
    assert!(stdout.contains("requires the reed fixture"));
    assert!(!output.exists());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn argument_and_sidecar_refusals_preserve_existing_evidence() {
    let dir = scratch("gesture-output-refusals");
    let out = dir.join("x.wav");
    for args in [
        vec!["string", out.to_str().unwrap(), "--full-scale-pa", "NaN"],
        vec!["string", out.to_str().unwrap(), "--full-scale-pa", "inf"],
        vec!["string", out.to_str().unwrap(), "--seconds", "1e-20"],
        vec!["string", out.to_str().unwrap(), "--bad\"\noption"],
    ] {
        let (ok, stdout) = run(&args);
        assert!(!ok, "invalid arguments must refuse");
        assert_eq!(
            stdout.lines().count(),
            1,
            "diagnostic contains a raw newline: {stdout}"
        );
        assert!(stdout.contains("\"verdict\":\"refused\""));
        assert!(!out.exists());
    }
    let sidecar = out.with_extension("provenance.json");
    std::fs::write(&sidecar, b"retained evidence").unwrap();
    let (ok, stdout) = run(&["string", out.to_str().unwrap()]);
    assert!(!ok);
    assert!(stdout.contains("refuses to overwrite evidence"));
    assert_eq!(std::fs::read(&sidecar).unwrap(), b"retained evidence");
    assert!(!out.exists());
    std::fs::remove_dir_all(&dir).unwrap();
}
