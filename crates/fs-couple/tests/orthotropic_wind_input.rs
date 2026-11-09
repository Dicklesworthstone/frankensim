//! Original material bytes through the shared wind loader, native physics and PCM.
//! All cards below are manufactured integration data, not instrument specimens.
use fs_couple::bernoulli_aperture::performance::file::material_modal::dg::MaterialModalDgValvePerformance;
use fs_couple::dg_port::input::wall::DgMaterialWallInput;
use fs_couple::flue_loop::dg::performance::file::FlueDgFilePerformance;
use fs_couple::pcm_wav::{encode_pcm16_wav, observation::{DecimatedRenderer, PressureRenderer}};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack,
    NormalizedPack, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use std::{path::PathBuf, process::Command, sync::Arc};

const TIME: &str = "normal_plane_viscosity_time";
const FLUE: &str = include_str!("../examples/flue-dg.performance");
const PAIRED: &str = include_str!("../examples/material_modal_valve/paired-dg.performance");
const WALL: &str = "dg_wall orthotropic wall.pack 293.15 1000 0.00002 strain-fixed 50 2000 1 2 proportional-kelvin-voigt normal_plane_viscosity_time\n";

fn pack(values: &[(&str, Dims, f64)]) -> NormalizedMaterialCardPack {
    pack_with_omega_limit(values, 100000.0)
}
fn pack_with_omega_limit(values: &[(&str, Dims, f64)], omega_max: f64) -> NormalizedMaterialCardPack {
    let artifact = fs_blake3::hash_domain("org.frankensim.test.orthotropic-wind", b"manufactured source");
    let mut claims = ClaimSet::new();
    for &(name, dims, value) in values {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, omega_max),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance {
                source: "manufactured directional wall/valve regression; no measured qualification".into(),
                license: "CC0-1.0".into(), artifact: Some(artifact),
            },
        }).unwrap();
    }
    let claims = NormalizedPack::new("manufactured-wind", "synthetic-v1", artifact,
        "CC0-1.0", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: "manufactured source".into(), phase: "solid".into(),
        process: "integration test".into(), revision: 0,
    }, claims).unwrap()
}
fn wall_pack(young: f64, density: f64, time: f64) -> NormalizedMaterialCardPack {
    pack(&[
        ("density", Density::DIMS, density),
        ("young_modulus_1", Pressure::DIMS, young),
        ("young_modulus_2", Pressure::DIMS, 4e6),
        ("poisson_ratio_12", Dims::NONE, 0.1),
        (TIME, Dims([0, 0, 1, 0, 0, 0]), time),
    ])
}
fn blade_pack() -> NormalizedMaterialCardPack {
    pack(&[
        ("density", Density::DIMS, 900.0),
        ("young_modulus", Pressure::DIMS, 4e9),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("eta11", DynViscosity::DIMS, 22000.0),
        ("eta22", DynViscosity::DIMS, 35000.0),
        ("eta66", DynViscosity::DIMS, 11000.0),
        ("eta12", DynViscosity::DIMS, 4500.0),
        ("eta16", DynViscosity::DIMS, 1000.0),
        ("eta26", DynViscosity::DIMS, -1500.0),
    ])
}
fn flue_source(declaration: &str) -> String {
    FLUE.replace("output 2048", "output 1024")
        .replace("\nflue ", &format!("\n{declaration}flue "))
}
fn flue(text: &str, material: &NormalizedMaterialCardPack, block: usize) -> FlueDgFilePerformance {
    FlueDgFilePerformance::from_bytes_with_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "wall.pack"); Ok(material.clone())
        }).unwrap()
}
fn trace(mut renderer: impl PressureRenderer, count: usize, block: usize) -> Vec<f64> {
    let mut pressure = vec![0.0; count];
    for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    pressure
}

#[test]
fn directional_grammar_requires_explicit_axes_loss_and_original_source_support() {
    let fields = WALL.trim_start_matches("dg_wall ");
    let gate = CancelGate::new_clock_free();
    let material = Arc::new(wall_pack(1e6, 1000.0, 2e-5));
    let read = |text: &str| DgMaterialWallInput::read(7, &mut text.split_ascii_whitespace());
    let source = read(fields).unwrap().unwrap().resolve(material.clone(), &gate).unwrap();
    assert_eq!(source.pack(), material.as_ref());
    assert_eq!(source.wall().orthotropic_axes(), Some([0, 1]));
    assert_eq!(source.wall().material().properties().len(), 5);
    assert_eq!(source.wall().proportional_viscosity_time_property(), Some(TIME));
    for invalid in [
        fields.replace("1 2 proportional", "0 2 proportional"),
        fields.replace("1 2 proportional", "1 4 proportional"),
        fields.replace("1 2 proportional", "2 2 proportional"),
        fields.replace("1 2 proportional", "1.0 2 proportional"),
        fields.replace("1 2 proportional", "18446744073709551615 2 proportional"),
        fields.replace("proportional-kelvin-voigt", "kelvin-voigt"),
        fields.replace("proportional-kelvin-voigt", "maxwell"),
        fields.replace(TIME, "density"),
        fields.replace(TIME, ""),
        fields.replace("wall.pack", "../wall.pack"),
        format!("{fields} extra"),
    ] { assert!(read(&invalid).is_err(), "{invalid}"); }
    let missing = fields.replace(TIME, "missing_time");
    assert!(read(&missing).unwrap().unwrap().resolve(material.clone(), &gate).is_err());
    let absent_axes = fields.replace("1 2 proportional", "1 3 proportional");
    assert!(read(&absent_axes).unwrap().unwrap().resolve(material.clone(), &gate).is_err());
    gate.request();
    assert!(read(fields).unwrap().unwrap().resolve(material, &gate).is_err());
}

#[test]
fn changing_material_direction_loss_or_gas_alters_actual_flue_pressure_not_only_metadata() {
    let text = flue_source(WALL);
    let material = wall_pack(1e6, 1000.0, 2e-5);
    let source = flue(&text, &material, 37);
    let duct = source.renderer().voice().dg().unwrap().load().duct();
    let pin = duct.compliant_wall().unwrap();
    let expected = 1e6 / (1.0 - 0.1 * 0.1 * 4.0) * 0.00002 / 0.001_f64.powi(2);
    assert!((pin.stiffness_per_area / expected - 1.0).abs() < 1e-13);
    assert!((pin.resistance / (2e-5 * expected) - 1.0).abs() < 1e-13);
    assert!(duct.has_viscothermal_memory());
    let reference = trace(source.into_renderer(), 1024, 37);
    assert_eq!(reference, trace(flue(&text, &material, 1).into_renderer(), 1024, 1));
    for changed in [wall_pack(2e6, 1000.0, 2e-5), wall_pack(1e6, 2000.0, 2e-5),
        wall_pack(1e6, 1000.0, 5e-5)] {
        assert_ne!(reference, trace(flue(&text, &changed, 37).into_renderer(), 1024, 37));
    }
    for changed in [text.replace("1 2 proportional", "2 1 proportional"),
        text.replace("N2 0.8 O2 0.2", "N2 0.7 CO2 0.3")] {
        assert_ne!(reference, trace(flue(&changed, &material, 37).into_renderer(), 1024, 37));
    }
}

#[test]
fn one_shared_loader_binds_the_same_directional_wall_to_multimode_reeds() {
    let material = wall_pack(1e6, 1000.0, 2e-5);
    let blade = blade_pack();
    let text = PAIRED.replace("audio 96000 96000 4800 20", "audio 96000 96000 512 20")
        .replace("slit_length", &format!("{WALL}slit_length"));
    let build = |text: &str, block| {
        let mut wall_loads = 0;
        let source = MaterialModalDgValvePerformance::from_bytes_with_sources(text.as_bytes(), block,
            &CancelGate::new_clock_free(), None, |path| {
                if path == "wall.pack" { wall_loads += 1; Ok(material.clone()) }
                else { assert_eq!(path, "blade.matpack"); Ok(blade.clone()) }
            }).unwrap();
        assert_eq!(wall_loads, 1);
        assert_eq!(source.wall_record().unwrap().pack(), &material);
        assert_eq!(source.material_records()[0].pack(), &blade);
        assert_eq!(source.info().modes, 6);
        source
    };
    let reference = trace(build(&text, 37).into_renderer(), 512, 37);
    assert_eq!(reference, trace(build(&text, 1).into_renderer(), 512, 1));
    let changed = text.replace("1 2 proportional", "2 1 proportional");
    assert_ne!(reference, trace(build(&changed, 37).into_renderer(), 512, 37));
}

#[test]
fn physical_regions_keep_local_direction_without_changing_the_source_artifact() {
    let material = wall_pack(1e6, 1000.0, 2e-5);
    let regions = "dg_wall regions 2\ndg_wall_region 0 0.025 orthotropic wall.pack 293.15 1000 0.00002 strain-fixed 50 2000 1 2 elastic\ndg_wall_region 0.025 0.05 orthotropic wall.pack 293.15 1000 0.00002 strain-fixed 50 2000 2 1 elastic\n";
    let text = flue_source(regions)
        .replace("open 1 2 8", "open 2 2 8");
    let mut loads = 0;
    let source = FlueDgFilePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "wall.pack"); loads += 1; Ok(material.clone())
        }).unwrap();
    assert_eq!(loads, 1, "one original pack serves both orientations");
    let layout = source.wall_layout();
    assert_eq!(layout.regions().len(), 2);
    assert_eq!(layout.regions()[0].source().wall().orthotropic_axes(), Some([0, 1]));
    assert_eq!(layout.regions()[1].source().wall().orthotropic_axes(), Some([1, 0]));
    assert_eq!(layout.regions()[0].source().pack(), layout.regions()[1].source().pack());
    let reference = trace(source.into_renderer(), 1024, 37);
    let changed = text.replace("2000 2 1 elastic", "2000 1 2 elastic");
    assert_ne!(reference, trace(flue(&changed, &material, 37).into_renderer(), 1024, 37));
    let misaligned = text.replace("0.025", "0.02");
    assert!(FlueDgFilePerformance::from_bytes_with_sources(misaligned.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(material.clone())).is_err());
}

struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) { let _ = std::fs::remove_dir_all(&self.0); }
}
#[test]
fn native_cli_and_attached_source_preserve_physical_pcm_and_directional_metadata() {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = Scratch(std::env::temp_dir().join(format!("fs-orthotropic-wind-{}-{stamp}", std::process::id())));
    std::fs::create_dir(&dir.0).unwrap();
    let input = dir.0.join("source.performance");
    let text = flue_source(WALL);
    std::fs::write(&input, &text).unwrap();
    let material = wall_pack(1e6, 1000.0, 2e-5);
    std::fs::write(dir.0.join("wall.pack"), material.to_bytes()).unwrap();
    let source = flue(&text, &material, 37);
    let info = source.info();
    let renderer = DecimatedRenderer::new(source.into_renderer(), info.sample_rate_hz, 48000, 37).unwrap();
    let count = renderer.output_samples_for(info.samples).unwrap() as usize;
    let expected = encode_pcm16_wav(&trace(renderer, count, 37), 48000, info.full_scale_pa).unwrap().0;
    assert!(expected[44..].iter().any(|byte| *byte != 0));
    let mut receipts = Vec::new();
    for attached in [false, true] {
        let output = dir.0.join(format!("out-{attached}.wav"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_music_render"));
        command.arg("wind").arg(&input).arg(&output).args(["--block", "37", "--decimate"]);
        if attached { command.args(["--self-contained", "--attached-source", "wall.pack"]); }
        let result = command.output().unwrap();
        assert!(result.status.success(), "{}{}", String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["orthotropic-principal-local-ring-membrane", "\"hoop\":1,\"axial\":2",
            "proportional-kelvin-voigt", "C_viscous=tau*C_elastic", TIME,
            "\"propagator\":\"original-exact-dg\"",
            &material.content_hash().to_hex()] {
            assert!(metadata.contains(token), "missing {token}: {metadata}");
        }
        assert!(!metadata.contains("\"model\":\"isotropic-local-ring-membrane\""));
        receipts.push(metadata);
    }
    assert_eq!(receipts[0], receipts[1]);
}

#[test]
fn reed_wind_and_ensemble_receipts_retain_the_original_wall_source() {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = Scratch(std::env::temp_dir().join(format!("fs-reed-wall-receipt-{}-{stamp}", std::process::id())));
    std::fs::create_dir(&dir.0).unwrap();
    let material = pack_with_omega_limit(&[
        ("density", Density::DIMS, 1000.0),
        ("young_modulus_1", Pressure::DIMS, 2e8),
        ("young_modulus_2", Pressure::DIMS, 4e8),
        ("poisson_ratio_12", Dims::NONE, 0.1),
        (TIME, Dims([0, 0, 1, 0, 0, 0]), 2e-5),
    ], 200000.0);
    let wall = "bore_wall orthotropic wall.pack 300 1000 0.00001 strain-fixed 0.000001 25000 1 2 proportional-kelvin-voigt normal_plane_viscosity_time\n";
    let source = DgMaterialWallInput::read(1,
        &mut wall.trim_start_matches("bore_wall ").split_ascii_whitespace())
        .unwrap().unwrap().resolve(Arc::new(material.clone()), &CancelGate::new_clock_free()).unwrap();
    let text = include_str!("../examples/reed-baffled-mixture.performance")
        .replace("audio 48000 4801 2", "audio 48000 4801 20000")
        .replace("baffled_outlet 0 0 0.5 4 16 1000", "listener_m 1")
        .replace("observer_law stokes-kirchhoff", "observer_law none")
        .replace("compile_limits", &format!("{wall}compile_limits"));
    let input = dir.0.join("source.performance");
    std::fs::write(&input, text).unwrap();
    std::fs::write(dir.0.join("wall.pack"), material.to_bytes()).unwrap();
    for ensemble in [false, true] {
        let output = dir.0.join(format!("out-{ensemble}.wav"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_music_render"));
        if ensemble {
            command.arg("ensemble").arg(&output).arg("--reed").arg(&input)
                .args(["--full-scale-pa", "20000"]);
        } else {
            command.arg("wind").arg(&input).arg(&output);
        }
        let result = command.args(["--block", "37"]).output().unwrap();
        assert!(result.status.success(), "{}{}", String::from_utf8_lossy(&result.stdout),
            String::from_utf8_lossy(&result.stderr));
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in [
            "\"bore_wall\":{\"model\":\"orthotropic-principal-local-ring-membrane\"",
            "\"pack_path\":\"wall.pack\"", "\"thickness_m\":1e-5",
            "\"axial_constraint\":\"strain-fixed\"", "\"band_hz\":[1e-6,2.5e4]",
            "\"hoop\":1,\"axial\":2", "proportional-kelvin-voigt", TIME,
            "\"propagator\":\"regional-tmm-characteristic-line\"",
            &material.content_hash().to_hex(), &source.wall().material().identity().to_hex(),
        ] {
            assert!(metadata.contains(token), "missing {token}: {metadata}");
        }
        assert!(!metadata.contains("\"propagator\":\"original-exact-dg\""), "{metadata}");
    }
}

#[path = "orthotropic_wind_input/normal.rs"]
mod normal;
