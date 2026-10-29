//! Physical material partitions survive mesh changes and reach native pressure.
use super::{CancelGate, FlueDgFilePerformance, NormalizedMaterialCardPack,
    PressureRenderer, DecimatedRenderer, TUBE, WALL, ThermochemicalPhase, ThermochemicalSources,
    directory, encode_pcm16_wav, gas_packs, gas_support, input, loaded, pack, run, success, trace};

fn declaration(first: &str, first_loss: &str, second: &str, second_loss: &str) -> String {
    format!("dg_wall regions 2\n\
dg_wall_region 0 0.025 isotropic {first} 293.15 1000 0.00002 stress-free 50 2000 {first_loss}\n\
dg_wall_region 0.025 0.05 isotropic {second} 293.15 1000 0.00002 stress-free 50 2000 {second_loss}")
}

fn regional(first_loss: &str, second_loss: &str) -> String {
    input("elastic")
        .replace(TUBE, &TUBE.replace("open 1 2", "open 2 2"))
        .replace(&format!("{WALL}elastic"), &declaration(
            "wall.fsmcdpk", first_loss, "second.fsmcdpk", second_loss))
}

fn build(
    text: &str, first: &NormalizedMaterialCardPack,
    second: &NormalizedMaterialCardPack, block: usize,
) -> FlueDgFilePerformance {
    build_in(text, first, second, block, None)
}

fn build_in(
    text: &str, first: &NormalizedMaterialCardPack,
    second: &NormalizedMaterialCardPack, block: usize, gas: Option<&ThermochemicalSources<'_>>,
) -> FlueDgFilePerformance {
    FlueDgFilePerformance::from_bytes_with_sources(
        text.as_bytes(), block, &CancelGate::new_clock_free(), gas,
        |path| match path {
            "wall.fsmcdpk" => Ok(first.clone()),
            "second.fsmcdpk" => Ok(second.clone()),
            _ => Err("missing explicit source".into()),
        },
    ).unwrap()
}

#[test]
fn equal_regions_reuse_one_artifact_and_preserve_the_uniform_operator_and_output() {
    let material = pack(1000.0, 1e6, None, 100_000.0);
    for loss in ["kelvin-voigt", "maxwell 1 1 0.5 boundary_modulus boundary_time"] {
        let uniform = input(loss).replace(TUBE, &TUBE.replace("open 1 2", "open 2 2"));
        let regional = uniform.replace(&format!("{WALL}{loss}"),
            &declaration("wall.fsmcdpk", loss, "wall.fsmcdpk", loss));
        let mut loads = 0;
        let source = FlueDgFilePerformance::from_bytes_with_sources(
            regional.as_bytes(), 37, &CancelGate::new_clock_free(), None, |path| {
                assert_eq!(path, "wall.fsmcdpk"); loads += 1; Ok(material.clone())
            },
        ).unwrap();
        assert_eq!(loads, 1);
        assert!(source.wall_record().is_none(), "regional sources have explicit physical support");
        let records = source.wall_layout().regions();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].interval_m(), [0.0, 0.025]);
        assert_eq!(records[1].interval_m(), [0.025, 0.05]);
        assert!(std::ptr::eq(records[0].source().pack(), records[1].source().pack()));
        let original = loaded(&uniform, &material, 37);
        assert_eq!(source.renderer().voice().dg().unwrap().load().duct().operator_identity(),
            original.renderer().voice().dg().unwrap().load().duct().operator_identity());
        assert_eq!(trace(source), trace(original));
    }
}

#[test]
fn physical_regions_preserve_local_sources_under_refinement_and_change_native_pressure() {
    let first = pack(1000.0, 1e6, None, 100_000.0);
    let second = pack(2000.0, 4e6, None, 100_000.0);
    let text = regional("kelvin-voigt", "kelvin-voigt");
    for elements in [2, 4] {
        let refined = text.replace("open 2 2", &format!("open {elements} 2"));
        let source = build(&refined, &first, &second, 37);
        let duct = source.renderer().voice().dg().unwrap().load().duct();
        assert!(duct.compliant_wall().is_none(), "different sources must not become an average wall");
        assert!(duct.has_viscothermal_memory());
        for element in 0..elements {
            let wall = duct.compliant_wall_at_element(element).unwrap();
            let (density, stiffness) = if element < elements / 2 {
                (0.02, 2e7)
            } else { (0.04, 8e7) };
            assert!((wall.surface_density / density - 1.0).abs() < 1e-14);
            assert!((wall.stiffness_per_area / stiffness - 1.0).abs() < 1e-14);
        }
        assert_eq!(source.wall_layout().regions()[1].interval_m(), [0.025, 0.05]);
    }
    let reference = trace(build(&text, &first, &second, 37));
    assert_eq!(reference, trace(build(&text, &first, &second, 1)));
    let changed = trace(build(&text, &first, &first, 37));
    assert_ne!(reference.0, changed.0, "a change to one original material must reach native pressure");
    assert_ne!(reference.1.displacement_m, changed.1.displacement_m);
}

#[test]
fn distinct_regional_histories_and_a_gas_field_reach_the_same_renderer() {
    let first = pack(1000.0, 1e6, None, 100_000.0);
    let second = pack(2000.0, 4e6, None, 100_000.0);
    let text = regional("maxwell 1 1 0.5 boundary_modulus boundary_time", "kelvin-voigt")
        .replace("dg_tube 0.05 0.001 open 2 2 8 50 2000 0.01 1e-11 256", concat!(
            "dg_bore 1 open 2 2 8 50 2000 0.01 1e-11 256 3 gas_field 2 mole 2 A B\n",
            "dg_segment cylinder 0.001 0.05\n",
            "gas_station 0 293.15 1 0\n",
            "gas_station 0.05 293.15 0.7 0.3",
        ));
    let packs = gas_packs();
    let bindings = packs.each_ref().map(gas_support::Packs::binding);
    let gas = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let source = build_in(&text, &first, &second, 37, Some(&gas));
    let duct = source.renderer().voice().dg().unwrap().load().duct();
    assert!(duct.gas_field().is_some());
    assert!(duct.spatial_viscothermal_admission().is_some());
    let reference = trace(source);
    assert_eq!(reference, trace(build_in(&text, &first, &second, 1, Some(&gas))));
    assert_eq!(reference.1.relaxation_pressure_pa.len(), 1);
    let pressure = &reference.1.relaxation_pressure_pa[0];
    assert!(pressure[..3].iter().any(|p| *p != 0.0));
    assert!(pressure[3..].iter().all(|p| *p == 0.0),
        "a material arm must not acquire support in another source region");
    assert!(reference.1.relaxation_stored_energy_j > 0.0);
}

#[test]
fn incomplete_partitions_cut_elements_missing_sources_and_complete_state_overflow_refuse() {
    let text = regional("elastic", "elastic");
    let material = pack(1000.0, 1e6, None, 100_000.0);
    for invalid in [
        text.replace("dg_wall regions 2", "dg_wall regions 0"),
        text.replace("dg_wall regions 2", "dg_wall regions 257"),
        text.replace("dg_wall_region 0 ", "dg_wall_region 0.001 "),
        text.replace("dg_wall_region 0.025 0.05", "dg_wall_region 0.026 0.05"),
        text.replace("dg_wall_region 0.025 0.05", "dg_wall_region 0.024 0.05"),
        text.replace("0.05 isotropic second", "0.0501 isotropic second"),
        text.replace("0.025", "0.02"),
        text.replace("1e-11 256", "1e-11 110"),
        text.replacen("stress-free 50 2000", "stress-free 100 2000", 1),
    ] {
        assert!(FlueDgFilePerformance::from_bytes_with_sources(
            invalid.as_bytes(), 37, &CancelGate::new_clock_free(), None,
            |_| Ok(material.clone()),
        ).is_err(), "admitted invalid source partition: {invalid}");
    }
    assert!(FlueDgFilePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None,
        |path| if path == "wall.fsmcdpk" { Ok(material.clone()) }
            else { Err("unavailable source".into()) },
    ).is_err());
}

#[test]
fn physical_source_regions_reach_the_cli_baffled_microphone_and_pcm() {
    let first = pack(1000.0, 1e6, None, 100_000.0);
    let second = pack(2000.0, 4e6, None, 100_000.0);
    let text = include_str!("../../../../../examples/flue-dg-baffled.performance")
        .replace("output 2048", "output 1024")
        .replace("baffled 1 2 8", "baffled 2 2 8")
        .replace("\nflue ", &format!("\n{}\nflue ",
            declaration("wall.fsmcdpk", "kelvin-voigt", "second.fsmcdpk", "kelvin-voigt")));
    let source = build(&text, &first, &second, 37);
    let info = source.info();
    let mut renderer = DecimatedRenderer::new(
        source.into_renderer(), info.sample_rate_hz, 48000, 37).unwrap();
    let mut pressure = vec![0.0; renderer.output_samples_for(info.samples).unwrap() as usize];
    for chunk in pressure.chunks_mut(37) { renderer.block(chunk).unwrap(); }
    let expected = encode_pcm16_wav(&pressure, 48000, info.full_scale_pa).unwrap().0;
    assert!(expected[44..].iter().any(|byte| *byte != 0));
    let dir = directory();
    let performance = dir.join("regional-wall.performance");
    std::fs::write(&performance, text).unwrap();
    std::fs::write(dir.join("wall.fsmcdpk"), first.to_bytes()).unwrap();
    std::fs::write(dir.join("second.fsmcdpk"), second.to_bytes()).unwrap();
    let output = dir.join("regional-wall.wav");
    success(run(&performance, &output, 37, true));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    for token in [
        "regional-isotropic-local-ring-membranes", "\"regions\":[", "\"interval_m\":[",
        "\"pack_path\":\"wall.fsmcdpk\"", "\"pack_path\":\"second.fsmcdpk\"",
        "\"pressure_point\":\"baffled-outlet\"", &second.content_hash().to_hex(),
    ] { assert!(metadata.contains(token), "{token}: {metadata}"); }
}
