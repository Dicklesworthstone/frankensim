//! Actual finite branch materials share the main-bore owner and native microphone.
use super::{CancelGate, FlueDgFilePerformance, NormalizedMaterialCardPack,
    PressureRenderer, DecimatedRenderer, directory, encode_pcm16_wav, pack, run, success};

const WALL: &str = "isotropic wall.fsmcdpk 293.15 1000 0.00002 stress-free 50 2000 maxwell 1 1 0.5 boundary_modulus boundary_time";

fn input(main: &str, chimney: &str) -> String {
    include_str!("../../../../../examples/flue-dg-chimney.performance")
        .replace("output 2048", "output 1024")
        .replace("\nflue ", &format!("\ndg_wall assembly 2\ndg_wall_main {main}\ndg_wall_chimney 7 {chimney}\nflue "))
}

fn build(text: &str, material: &NormalizedMaterialCardPack, block: usize) -> FlueDgFilePerformance {
    FlueDgFilePerformance::from_bytes_with_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "wall.fsmcdpk"); Ok(material.clone())
        }).unwrap()
}

fn trace(source: FlueDgFilePerformance) -> (Vec<f64>, fs_duct::dg::DgCompliantWallState) {
    let samples = source.info().samples as usize;
    let block = source.renderer().max_block_len();
    let mut renderer = source.into_renderer();
    let mut pressure = vec![0.0; samples];
    for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
    let branch = &renderer.voice().dg().unwrap().load().branches()[0];
    let wall = branch.duct().compliant_wall_state(branch.state()).unwrap();
    assert!(wall.stored_energy_j > 0.0);
    assert!(wall.relaxation_stored_energy_j > 0.0);
    assert!(wall.displacement_m.iter().any(|x| *x != 0.0));
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-12));
    (pressure, wall)
}

#[test]
fn finite_branch_reuses_the_main_artifact_and_keeps_its_own_gas_and_material_memory() {
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let text = input(WALL, WALL);
    let mut loads = 0;
    let source = FlueDgFilePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "wall.fsmcdpk"); loads += 1; Ok(material.clone())
        }).unwrap();
    assert_eq!(loads, 1);
    let assignment = &source.wall_layout().chimneys()[0];
    assert_eq!(assignment.port(), fs_duct::dg::PortId(7));
    assert!(std::ptr::eq(source.wall_record().unwrap().pack(),
        assignment.layout().uniform_record().unwrap().pack()));
    let load = source.renderer().voice().dg().unwrap().load();
    let duct = load.branches()[0].duct();
    assert_ne!(duct.inlet_gas(), load.duct().inlet_gas());
    assert!(duct.has_viscothermal_memory());
    let wall = duct.compliant_wall().unwrap();
    // Actual branch radius, never the main bore's radius or an invented scale.
    assert!((wall.surface_density - 1000.0 * 0.00002).abs() < 1e-15);
    assert!((wall.stiffness_per_area / (1e6 * 0.00002 / 0.0006_f64.powi(2)) - 1.0).abs() < 1e-14);
    let reference = trace(source);
    assert_eq!(reference, trace(build(&text, &material, 1)));

    // Hold the main bore rigid so only this branch's material changes.
    let text = input("rigid", WALL);
    let reference = trace(build(&text, &material, 37));
    for altered in [pack(2000.0, 1e6, None, 100_000.0), pack(1000.0, 4e6, None, 100_000.0)] {
        let other = trace(build(&text, &altered, 37));
        assert_ne!(reference.0, other.0);
        assert_ne!(reference.1.displacement_m, other.1.displacement_m);
    }
    let changed_gas = text.replace("N2 0.6 CO2 0.4", "N2 0.9 CO2 0.1");
    let other = trace(build(&changed_gas, &material, 37));
    assert_ne!(reference.0, other.0);
    assert_ne!(reference.1.relaxation_pressure_pa, other.1.relaxation_pressure_pa);
}

#[test]
fn branch_regions_use_local_coordinates_and_unknown_or_unbounded_assemblies_refuse() {
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let declaration = format!("regions 2\ndg_wall_region 0 0.006 {WALL}\ndg_wall_region 0.006 0.012 {}",
        WALL.replace("0.00002 stress-free", "0.00001 stress-free"));
    let text = input("rigid", &declaration).replace("0.0006 0.012 1 2", "0.0006 0.012 2 2");
    let source = build(&text, &material, 37);
    let branch = source.renderer().voice().dg().unwrap().load().branches()[0].duct();
    assert!(branch.compliant_wall().is_none());
    assert_eq!(branch.compliant_wall_at_element(1).unwrap().surface_density,
        0.5 * branch.compliant_wall_at_element(0).unwrap().surface_density);
    assert_eq!(source.wall_layout().chimneys()[0].layout().regions()[1].interval_m(), [0.006, 0.012]);
    assert_eq!(trace(source), trace(build(&text, &material, 1)));

    for invalid in [
        text.replace("assembly 2", "assembly 0"),
        text.replace("dg_wall_main rigid", "dg_wall_chimney 7 rigid"),
        text.replace("dg_wall_chimney 7", "dg_wall_chimney 99"),
        text.replace("dg_wall_main rigid", "dg_wall_main assembly 1"),
        text.replace("0.006", "0.007"),
        text.replace("0.012 isotropic", "0.013 isotropic"),
        text.replace("1e-11 256", "1e-11 220"),
        text.replacen("stress-free 50 2000", "stress-free 100 2000", 1),
    ] {
        assert!(FlueDgFilePerformance::from_bytes_with_sources(invalid.as_bytes(), 37,
            &CancelGate::new_clock_free(), None, |_| Ok(material.clone())).is_err(),
            "admitted invalid branch wall: {invalid}");
    }
}

#[test]
fn finite_chimney_material_reaches_cli_microphone_pcm_and_original_source_metadata() {
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let text = input("rigid", WALL);
    let source = build(&text, &material, 37);
    let info = source.info();
    let mut renderer = DecimatedRenderer::new(source.into_renderer(),
        info.sample_rate_hz, 48000, 37).unwrap();
    let mut pressure = vec![0.0; renderer.output_samples_for(info.samples).unwrap() as usize];
    for chunk in pressure.chunks_mut(37) { renderer.block(chunk).unwrap(); }
    let expected = encode_pcm16_wav(&pressure, 48000, info.full_scale_pa).unwrap().0;
    let dir = directory();
    let input = dir.join("branch-wall.performance");
    std::fs::write(&input, text).unwrap();
    std::fs::write(dir.join("wall.fsmcdpk"), material.to_bytes()).unwrap();
    let output = dir.join("branch-wall.wav");
    success(run(&input, &output, 37, true));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    for token in ["\"solid_wall\":null", "\"chimney_walls\":[{\"port\":7,\"source\":{",
        "\"pack_path\":\"wall.fsmcdpk\"", "\"loss\":\"maxwell\"", "\"propagator\":\"original-exact-dg\""] {
        assert!(metadata.contains(token), "{token}: {metadata}");
    }
}
