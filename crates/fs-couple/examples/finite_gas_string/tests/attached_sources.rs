//! The attached runner uses the original source decoders and physical trajectory.
use super::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT: AtomicU64 = AtomicU64::new(0);
const CORE: &str = "materials/core.pack";
const LOCAL: &str = "gas/local.gas";
const RESERVOIR: &str = "gas/reservoir.gas";
const SAMPLES: usize = 24;

struct Fixture {
    directory: PathBuf,
    sources: Vec<String>,
}

impl Fixture {
    fn new(rate: u32) -> Self {
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "fs-attached-finite-gas-{}-{stamp}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&directory).unwrap();
        for name in ["materials", "gas"] {
            std::fs::create_dir(directory.join(name)).unwrap();
        }
        let pack = NormalizedPack::new(
            "attached-live-string",
            "integration-v1",
            fs_blake3::ContentHash([17; 32]),
            "CC0-1.0",
            source_law(1.0, None).claims().clone(),
            vec![],
            vec![],
        )
        .unwrap();
        std::fs::write(directory.join(CORE), pack.to_bytes()).unwrap();
        let mut sources = vec![CORE.to_owned(), LOCAL.to_owned(), RESERVOIR.to_owned()];
        for (name, mass, cp_r) in [("heavy", 0.028, 3.5), ("light", 0.004, 2.5)] {
            let packs = gas_support::packs(name, mass, cp_r, 200.0);
            for (suffix, bytes) in [
                ("species", packs.species.to_bytes()),
                ("nasa9", packs.nasa9.to_bytes()),
                ("transport", packs.transport.to_bytes()),
            ] {
                let path = format!("gas/{name}.{suffix}");
                std::fs::write(directory.join(&path), bytes).unwrap();
                sources.push(path);
            }
        }
        std::fs::write(
            directory.join(LOCAL),
            "thermochemical 300 85000 mole 1 heavy 1",
        )
        .unwrap();
        std::fs::write(
            directory.join(RESERVOIR),
            "thermochemical 300 95000 mole 1 light 1",
        )
        .unwrap();
        let scenario = include_str!("../feedback-flow.string")
            .replace(
                "audio 16000 3200 0.01",
                &format!("audio {rate} {SAMPLES} 0.01"),
            )
            .replace("string 0.65 0.00025 20 4", "string 0.65 0.00025 20 1")
            .replace("release 0.4 0.001 16", "release 0.4 0.001 2")
            .replace("nitrogen-reservoir.gas", RESERVOIR);
        std::fs::write(directory.join("input.string"), scenario).unwrap();
        Self { directory, sources }
    }

    fn output(&self, stem: &str, extension: &str) -> PathBuf {
        self.directory.join(stem).with_extension(extension)
    }

    fn args(&self, stem: &str, attached: Option<&[String]>) -> Vec<String> {
        let name = |name: &str| {
            if attached.is_some() {
                name.to_owned()
            } else {
                self.directory.join(name).to_str().unwrap().to_owned()
            }
        };
        let mut args = vec![
            self.directory
                .join("input.string")
                .to_str()
                .unwrap()
                .to_owned(),
            name(CORE),
            self.output(stem, "wav").to_str().unwrap().to_owned(),
            "--history-output".into(),
            self.output(stem, "csv").to_str().unwrap().to_owned(),
            "--gas-record".into(),
            name(LOCAL),
            "--gas-water".into(),
            "none".into(),
            "--gas-phase".into(),
            "frozen".into(),
        ];
        for component in ["heavy", "light"] {
            args.push("--gas-component".into());
            for suffix in ["species", "nasa9", "transport"] {
                args.push(name(&format!("gas/{component}.{suffix}")));
            }
        }
        if let Some(sources) = attached {
            args.push("--self-contained".into());
            for source in sources {
                args.extend(["--attached-source".into(), source.clone()]);
            }
        }
        args
    }

    fn no_output(&self, stem: &str) {
        for extension in ["wav", "csv"] {
            assert!(!self.output(stem, extension).exists());
        }
    }
}

#[test]
fn attached_feedback_and_nozzle_preserve_native_audio_and_complete_physical_history() {
    let fixture = Fixture::new(48_000);
    for (stem, sources) in [
        ("ordinary", None),
        ("attached", Some(fixture.sources.as_slice())),
    ] {
        let mut args = fixture.args(stem, sources);
        args.extend(["--require-sample-rate".into(), "48000".into()]);
        run(&args).unwrap();
    }
    let wav = std::fs::read(fixture.output("attached", "wav")).unwrap();
    assert_eq!(
        wav,
        std::fs::read(fixture.output("ordinary", "wav")).unwrap()
    );
    assert_eq!(u32::from_le_bytes(wav[24..28].try_into().unwrap()), 48_000);
    assert_eq!(wav.len(), 44 + 2 * SAMPLES);
    let csv = std::fs::read_to_string(fixture.output("attached", "csv")).unwrap();
    assert_eq!(
        csv,
        std::fs::read_to_string(fixture.output("ordinary", "csv")).unwrap()
    );
    let mut lines = csv.lines();
    let headers = lines.next().unwrap().split(',').collect::<Vec<_>>();
    let rows = lines
        .map(|line| {
            line.split(',')
                .map(|v| v.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), SAMPLES);
    assert!(rows.iter().all(|row| row.len() == headers.len()));
    assert!(rows.iter().all(|row| row[15].abs() <= row[16]));
    assert!(rows.iter().any(|row| row[1] != 0.0));
    let last = rows.last().unwrap();
    let value = |name| last[headers.iter().position(|header| *header == name).unwrap()];
    assert!(value("solid_temperature_k") > 300.0);
    assert!(value("cumulative_mass_into_surroundings_kg") > 0.0);
    assert!(value("cumulative_electrical_source_work_j") > 0.0);
    assert!(value("electrical_feedback_trials") >= 1.0);
    assert!(headers.contains(&"cumulative_electrical_source_residual_j"));
    assert!(
        headers
            .iter()
            .any(|name| name.contains("gas_mole_fraction_light"))
    );
}

#[test]
fn unlisted_solid_gas_and_reservoir_sources_never_fall_back_to_existing_files() {
    let fixture = Fixture::new(48_000);
    for (i, missing) in [
        CORE,
        LOCAL,
        RESERVOIR,
        "gas/heavy.species",
        "gas/heavy.nasa9",
        "gas/heavy.transport",
    ]
    .iter()
    .enumerate()
    {
        assert!(fixture.directory.join(missing).is_file());
        let names = fixture
            .sources
            .iter()
            .filter(|name| name.as_str() != *missing)
            .cloned()
            .collect::<Vec<_>>();
        let stem = format!("unlisted-{i}");
        let error = run(&fixture.args(&stem, Some(&names)))
            .unwrap_err()
            .to_string();
        assert!(
            error.contains(*missing) && error.contains("unlisted source"),
            "{error}"
        );
        fixture.no_output(&stem);
    }

    // Bonded references use the same admitted map, even when the sleeve file
    // already exists beside its scenario and would be readable ordinarily.
    let sleeve = "materials/sleeve.pack";
    std::fs::write(
        fixture.directory.join(sleeve),
        std::fs::read(fixture.directory.join(CORE)).unwrap(),
    )
    .unwrap();
    let bonded = include_str!("../bonded.string")
        .replace("audio 16000 3200 0.01", "audio 48000 24 0.01")
        .replace("sleeve.pack", sleeve);
    std::fs::write(fixture.directory.join("input.string"), bonded).unwrap();
    let error = run(&fixture.args("unlisted-sleeve", Some(&fixture.sources)))
        .unwrap_err()
        .to_string();
    assert!(
        error.contains(sleeve) && error.contains("unlisted source"),
        "{error}"
    );
    fixture.no_output("unlisted-sleeve");
}

#[test]
fn existing_audio_or_history_is_preserved_without_creating_its_companion() {
    let fixture = Fixture::new(48_000);
    let previous = b"previous accepted output";
    for (existing, absent) in [("wav", "csv"), ("csv", "wav")] {
        let stem = format!("existing-{existing}");
        let path = fixture.output(&stem, existing);
        std::fs::write(&path, previous).unwrap();
        assert!(run(&fixture.args(&stem, Some(&fixture.sources))).is_err());
        assert_eq!(std::fs::read(path).unwrap(), previous);
        assert!(!fixture.output(&stem, absent).exists());
    }
}

#[test]
fn required_clock_refuses_before_attached_source_io_and_never_resamples() {
    let fixture = Fixture::new(16_000);
    let mut args = fixture.args("wrong-clock", Some(&["missing-source.pack".into()]));
    args.extend(["--require-sample-rate".into(), "48000".into()]);
    let error = run(&args).unwrap_err().to_string();
    assert!(
        error.contains("16000") && error.contains("48000"),
        "{error}"
    );
    fixture.no_output("wrong-clock");
}

#[test]
fn csv_byte_budget_admits_the_complete_original_history_or_refuses() {
    scoped(|cx, _| {
        let mut spec = scenario();
        spec.samples = 4;
        spec.release_sample = 1;
        let gas = parcel("history-source", 0.028, 3.5, spec.gas_volume_m3);
        let history = simulate(cx, &spec, &source(1.0), gas, None).unwrap();
        let mut streamed = Vec::new();
        output::write_csv(&history, &mut streamed).unwrap();
        assert_eq!(
            output::csv_bytes(&history, streamed.len()).unwrap(),
            streamed
        );
        let error = output::csv_bytes(&history, streamed.len() - 1)
            .unwrap_err()
            .to_string();
        assert!(error.contains("history byte limit"), "{error}");
    });
}

#[test]
fn incompatible_clock_and_repeated_output_flags_refuse_without_artifacts() {
    let fixture = Fixture::new(48_000);
    for (i, tail) in [
        vec!["--decimate"],
        vec!["--history-output", "another.csv"],
        vec![
            "--require-sample-rate",
            "48000",
            "--require-sample-rate",
            "48000",
        ],
    ]
    .iter()
    .enumerate()
    {
        let stem = format!("invalid-option-{i}");
        let mut args = fixture.args(&stem, Some(&fixture.sources));
        args.extend(tail.iter().map(|value| (*value).to_owned()));
        let error = run(&args).unwrap_err().to_string();
        assert!(
            error.contains(if i == 0 {
                "--decimate"
            } else {
                "cannot be repeated"
            }),
            "{error}"
        );
        fixture.no_output(&stem);
    }
}
