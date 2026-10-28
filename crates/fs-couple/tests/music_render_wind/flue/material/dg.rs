//! Source cards drive the same native DG mechanics for every wind source.
use super::{CancelGate, NormalizedMaterialCardPack, directory, pack, run, success};
use fs_couple::dg_port::input::wall::DgMaterialWallInput;
use fs_couple::flue_loop::dg::performance::file::FlueDgFilePerformance;
use fs_couple::pcm_wav::{
    encode_pcm16_wav,
    observation::{DecimatedRenderer, PressureRenderer},
};

const WALL: &str = "dg_wall isotropic wall.fsmcdpk 293.15 1000 0.00002 stress-free 50 2000 ";
const TUBE: &str = "dg_tube 0.05 0.001 open 1 2 8 50 2000 0.01 1e-11 256";

fn input(loss: &str) -> String {
    include_str!("../../../../examples/flue-dg.performance")
        .replace("output 2048", "output 256")
        .replace("\nflue ", &format!("\n{WALL}{loss}\nflue "))
}

fn loaded(
    text: &str,
    material: &NormalizedMaterialCardPack,
    block: usize,
) -> FlueDgFilePerformance {
    FlueDgFilePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            assert_eq!(path, "wall.fsmcdpk");
            Ok(material.clone())
        },
    )
    .unwrap()
}

fn trace(source: FlueDgFilePerformance) -> (Vec<f64>, fs_duct::dg::DgCompliantWallState) {
    let samples = source.info().samples as usize;
    let mut renderer = source.into_renderer();
    let mut pressure = vec![0.0; samples];
    for chunk in pressure.chunks_mut(renderer.max_block_len()) {
        renderer.block(chunk).unwrap();
    }
    let load = renderer.voice().dg().unwrap().load();
    assert_eq!(load.state().step_count(), samples as u64);
    let wall = load.duct().compliant_wall_state(load.state()).unwrap();
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-12));
    assert!(wall.stored_energy_j > 0.0);
    assert!(wall.displacement_m.iter().any(|x| *x != 0.0));
    (pressure, wall)
}

#[test]
fn source_density_stiffness_and_loss_drive_the_actual_dg_wall_and_jet_pressure() {
    let material = pack(1000.0, 1e6, None, 100_000.0);
    for loss in [
        "elastic",
        "kelvin-voigt",
        "maxwell 1 1 0.5 boundary_modulus boundary_time",
    ] {
        let text = input(loss);
        let source = loaded(&text, &material, 37);
        let record = source.wall_record().unwrap();
        assert_eq!(record.pack(), &material);
        assert_eq!(
            record.wall().material().card_identity(),
            material.card().content_hash()
        );
        let duct = source.renderer().voice().dg().unwrap().load().duct();
        let wall = duct.compliant_wall().unwrap();
        // Independent thin-cylinder membrane reduction: sigma=rho*h and K=E*h/r².
        assert!((wall.surface_density - 0.02).abs() < 1e-15);
        assert!((wall.stiffness_per_area - 2e7).abs() < 1e-7);
        assert_eq!(
            wall.resistance,
            if loss == "kelvin-voigt" { 40.0 } else { 0.0 }
        );
        assert!(
            duct.has_viscothermal_memory(),
            "gas and solid losses are separate owners"
        );
        let reference = trace(source);
        assert_eq!(reference, trace(loaded(&text, &material, 1)));
        if loss.starts_with("maxwell") {
            assert_eq!(reference.1.relaxation_pressure_pa.len(), 1);
            assert!(reference.1.relaxation_stored_energy_j > 0.0);
            assert!(reference.1.dissipated_power_w > 0.0);
        }
        for changed in [
            pack(2000.0, 1e6, None, 100_000.0),
            pack(1000.0, 4e6, None, 100_000.0),
        ] {
            let changed = trace(loaded(&text, &changed, 37));
            assert_ne!(
                reference.0, changed.0,
                "source-only change must alter returning bore pressure"
            );
            assert_ne!(reference.1.displacement_m, changed.1.displacement_m);
        }
    }
}

#[test]
fn nonuniform_composition_and_the_solid_source_share_one_original_operator() {
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let text = input("maxwell 1 1 0.5 boundary_modulus boundary_time").replace(TUBE,
        "dg_bore 1 open 2 1 8 50 2000 0.01 1e-11 256 2 gas_field 2 mole 2 N2 CO2\ndg_segment cylinder 0.001 0.05\ngas_station 0 293.15 1 0\ngas_station 0.05 293.15 1 0");
    let changed = text.replace(
        "gas_station 0.05 293.15 1 0",
        "gas_station 0.05 293.15 0.7 0.3",
    );
    let a = loaded(&text, &material, 37);
    let b = loaded(&changed, &material, 37);
    let a_duct = a.renderer().voice().dg().unwrap().load().duct();
    let b_duct = b.renderer().voice().dg().unwrap().load().duct();
    assert_eq!(a_duct.inlet_gas(), b_duct.inlet_gas());
    assert_ne!(a_duct.outlet_gas(), b_duct.outlet_gas());
    assert_eq!(a_duct.compliant_wall(), b_duct.compliant_wall());
    assert!(b_duct.gas_field().is_some());
    assert_ne!(a_duct.operator_identity(), b_duct.operator_identity());
    let a = trace(a);
    let b = trace(b);
    assert_ne!(a.0, b.0);
    assert_ne!(a.1.displacement_m, b.1.displacement_m);
    assert_ne!(a.1.relaxation_pressure_pa, b.1.relaxation_pressure_pa);
}

#[test]
fn sourced_dg_wall_reaches_the_actual_cli_microphone_and_pcm_owner() {
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let text = include_str!("../../../../examples/flue-dg-baffled.performance")
        .replace("output 2048", "output 1024")
        .replace("\nflue ", &format!("\n{WALL}kelvin-voigt\nflue "));
    let source = loaded(&text, &material, 37);
    let info = source.info();
    let mut renderer =
        DecimatedRenderer::new(source.into_renderer(), info.sample_rate_hz, 48000, 37).unwrap();
    let mut pressure = vec![0.0; renderer.output_samples_for(info.samples).unwrap() as usize];
    for chunk in pressure.chunks_mut(37) {
        renderer.block(chunk).unwrap();
    }
    let expected = encode_pcm16_wav(&pressure, 48000, info.full_scale_pa)
        .unwrap()
        .0;
    assert!(expected[44..].iter().any(|byte| *byte != 0));
    let dir = directory();
    let path = dir.join("flue-dg-wall.performance");
    std::fs::write(&path, &text).unwrap();
    std::fs::write(dir.join("wall.fsmcdpk"), material.to_bytes()).unwrap();
    let output = dir.join("flue-dg-wall.wav");
    success(run(&path, &output, 37, true));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let restricted = dir.join("flue-dg-wall-attached.wav");
    success(std::process::Command::new(env!("CARGO_BIN_EXE_music_render"))
        .arg("wind").arg(&path).arg(&restricted)
        .args(["--block", "37", "--decimate", "--self-contained",
            "--attached-source", "wall.fsmcdpk"]).output().unwrap());
    assert_eq!(std::fs::read(&restricted).unwrap(), expected,
        "explicit source snapshots must retain the original physical PCM");
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert_eq!(std::fs::read_to_string(restricted.with_extension("provenance.json")).unwrap(), metadata);
    for token in [
        "\"solid_wall\":{",
        "\"loss\":\"kelvin-voigt\"",
        "\"pressure_point\":\"baffled-outlet\"",
        "\"propagator\":\"original-exact-dg\"",
        &material.content_hash().to_hex(),
    ] {
        assert!(metadata.contains(token), "{token}: {metadata}");
    }
}

#[test]
fn rigid_is_exactly_the_original_path_and_invalid_sources_never_default_to_it() {
    let original = include_str!("../../../../examples/flue-dg.performance");
    let rigid = original.replace("\nflue ", "\ndg_wall rigid\nflue ");
    let a =
        FlueDgFilePerformance::from_bytes(original.as_bytes(), 37, &CancelGate::new_clock_free())
            .unwrap();
    let b = FlueDgFilePerformance::from_bytes_with_sources(
        rigid.as_bytes(),
        37,
        &CancelGate::new_clock_free(),
        None,
        |_| panic!("rigid does not resolve a solid"),
    )
    .unwrap();
    assert!(b.wall_record().is_none());
    assert_eq!(
        a.renderer()
            .voice()
            .dg()
            .unwrap()
            .load()
            .duct()
            .operator_identity(),
        b.renderer()
            .voice()
            .dg()
            .unwrap()
            .load()
            .duct()
            .operator_identity()
    );
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let text = input("elastic");
    assert!(
        FlueDgFilePerformance::from_bytes(text.as_bytes(), 37, &CancelGate::new_clock_free())
            .is_err()
    );
    for invalid in [
        text.replace("0.00002 stress-free", "0.00006 stress-free"),
        text.replace("stress-free 50 2000", "stress-free 100 2000"),
        text.replace("1e-11 256", "1e-11 54"),
        text.replace(
            TUBE,
            "dg_bore 1 open 1 2 12 50 2000 0.01 1e-11 256 3\ndg_segment cone 0.001 0.002 0.05",
        ),
    ] {
        assert!(
            FlueDgFilePerformance::from_bytes_with_sources(
                invalid.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(material.clone())
            )
            .is_err()
        );
    }
    let missing = pack(1000.0, 1e6, Some("young_modulus"), 100_000.0);
    assert!(
        FlueDgFilePerformance::from_bytes_with_sources(
            text.as_bytes(),
            37,
            &CancelGate::new_clock_free(),
            None,
            |_| Ok(missing.clone())
        )
        .is_err()
    );
    for row in [
        "rigid extra",
        "isotropic ../wall.fsmcdpk 293.15 0 0.001 stress-free 50 2000 elastic",
        "isotropic wall.fsmcdpk 293.15 0 0.001 stress-free 50 2000 unknown",
        "isotropic wall.fsmcdpk 293.15 NaN 0.001 stress-free 50 2000 elastic",
    ] {
        assert!(DgMaterialWallInput::read(1, &mut row.split_ascii_whitespace()).is_err());
    }
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(
        FlueDgFilePerformance::from_bytes_with_sources(
            text.as_bytes(),
            37,
            &gate,
            None,
            |_| panic!("cancelled preparation must not load material")
        )
        .is_err()
    );
}

#[path = "dg/regions.rs"]
mod regions;

#[path = "dg/chimneys.rs"]
mod chimneys;
