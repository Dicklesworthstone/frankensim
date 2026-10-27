//! Independent normal viscosity through the ORIGINAL wind loaders and PCM path.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::WallMaterialLoss;

const NAMES: [&str; 3] = ["wall_eta_hh", "wall_eta_ha", "wall_eta_aa"];
const NORMAL: &str = "dg_wall orthotropic wall.pack 293.15 1000 0.00002 stress-free 50 2000 1 2 normal-kelvin-voigt 1 0.8 wall_eta_hh wall_eta_ha wall_eta_aa\n";
fn normal_pack(cross: f64) -> NormalizedMaterialCardPack {
    pack(&[
        ("density", Density::DIMS, 1000.0),
        ("young_modulus_1", Pressure::DIMS, 1e6),
        ("young_modulus_2", Pressure::DIMS, 4e6),
        ("poisson_ratio_12", Dims::NONE, 0.1),
        (NAMES[0], DynViscosity::DIMS, 400.0),
        (NAMES[1], DynViscosity::DIMS, cross),
        (NAMES[2], DynViscosity::DIMS, 90.0),
    ])
}

#[test]
fn normal_viscosity_grammar_retains_original_source_names_and_refuses_missing_data() {
    let read = |s: &str| DgMaterialWallInput::read(3, &mut s.split_ascii_whitespace());
    let fields = NORMAL.trim_start_matches("dg_wall ");
    let material = Arc::new(normal_pack(-60.0));
    let source = read(fields).unwrap().unwrap()
        .resolve(material.clone(), &CancelGate::new_clock_free()).unwrap();
    assert_eq!(source.normal_viscosity_properties(), Some(NAMES));
    assert_eq!(source.pack(), material.as_ref());
    assert_eq!(source.wall().material().properties().len(), 7);
    assert_eq!(source.wall().proportional_viscosity_time_property(), None);
    let WallMaterialLoss::Maxwell(memory) = source.wall().loss() else {
        panic!("zero axial stress needs the exact induced material history");
    };
    let c_aa = 4e6 / 0.96;
    let (modulus, time) = memory.material().terms[0];
    assert!((modulus / (c_aa * (0.1 + 60.0 / 90.0_f64).powi(2)) - 1.0).abs() < 1e-13);
    assert!((time / (90.0 / c_aa) - 1.0).abs() < 1e-13);
    assert!(memory.provenance().contains("exact normal-plane"));
    assert!(memory.provenance().contains(NAMES[1]));
    for bad in [
        fields.replace("normal-kelvin-voigt 1", "normal-kelvin-voigt 0"),
        fields.replace("normal-kelvin-voigt 1", "normal-kelvin-voigt 2.1"),
        fields.replace("1 0.8 wall_eta", "1 1.1 wall_eta"),
        fields.replace("1 0.8 wall_eta", "1 NaN wall_eta"),
        fields.replace(NAMES[2], ""),
        fields.replace(NAMES[0], "density"),
        format!("{fields} extra"),
    ] { assert!(read(&bad).is_err(), "{bad}"); }
    for missing in NAMES {
        let unresolved = fields.replace(missing, "absent_source");
        assert!(read(&unresolved).unwrap().unwrap()
            .resolve(material.clone(), &CancelGate::new_clock_free()).is_err());
    }
    assert!(read(fields).unwrap().unwrap()
        .resolve(Arc::new(normal_pack(300.0)), &CancelGate::new_clock_free()).is_err());
    let gate = CancelGate::new_clock_free(); gate.request();
    assert!(read(fields).unwrap().unwrap().resolve(material, &gate).is_err());
}

#[test]
fn cross_viscosity_and_gas_changes_reach_actual_flue_pressure_and_retained_memory() {
    let text = flue_source(NORMAL);
    let material = normal_pack(-60.0);
    let source = flue(&text, &material, 37);
    let duct = source.renderer().voice().dg().unwrap().load().duct();
    assert!(duct.has_viscothermal_memory());
    let base_pin = duct.compliant_wall().unwrap();
    let positive = flue(&text, &normal_pack(60.0), 37);
    assert_eq!(base_pin, positive.renderer().voice().dg().unwrap().load().duct().compliant_wall().unwrap(),
        "changing only the sign of D_ha changes memory, not residual Kelvin resistance");
    let mut renderer = source.into_renderer();
    let mut reference = vec![0.0; 1024];
    for chunk in reference.chunks_mut(37) { renderer.block(chunk).unwrap(); }
    let load = renderer.voice().dg().unwrap().load();
    let motion = load.duct().compliant_wall_state(load.state()).unwrap();
    assert_eq!(motion.relaxation_pressure_pa.len(), 1);
    assert!(motion.relaxation_stored_energy_j > 0.0);
    assert!(motion.dissipated_power_w > 0.0);
    assert_eq!(reference, trace(flue(&text, &material, 1).into_renderer(), 1024, 1));
    assert_ne!(reference, trace(positive.into_renderer(), 1024, 37));
    for changed in [text.replace("stress-free", "strain-fixed"),
        text.replace("N2 0.8 O2 0.2", "N2 0.7 CO2 0.3")] {
        assert_ne!(reference, trace(flue(&changed, &material, 37).into_renderer(), 1024, 37));
    }
}

#[test]
fn regional_axial_constraints_share_one_pack_and_keep_memory_on_its_physical_support() {
    let properties = NORMAL.trim_start_matches("dg_wall ").trim();
    let regions = format!("dg_wall regions 2\ndg_wall_region 0 0.025 {properties}\ndg_wall_region 0.025 0.05 {}\n",
        properties.replace("stress-free", "strain-fixed"));
    let text = flue_source(&regions).replace("open 1 2 8", "open 2 2 8");
    let material = normal_pack(-60.0);
    let mut loads = 0;
    let source = FlueDgFilePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path,"wall.pack"); loads += 1; Ok(material.clone())
        }).unwrap();
    assert_eq!(loads, 1);
    for region in source.wall_layout().regions() {
        assert_eq!(region.source().normal_viscosity_properties(), Some(NAMES));
        assert_eq!(region.source().pack(), &material);
    }
    let mut renderer = source.into_renderer();
    renderer.block(&mut [0.0;37]).unwrap();
    let load = renderer.voice().dg().unwrap().load();
    let motion = load.duct().compliant_wall_state(load.state()).unwrap();
    assert_eq!(motion.relaxation_pressure_pa.len(),1);
    assert!(motion.relaxation_pressure_pa[0][..3].iter().any(|p| *p != 0.0));
    assert!(motion.relaxation_pressure_pa[0][3..].iter().all(|p| *p == 0.0),
        "no artificial hereditary state outside the free-axial region");
    let misaligned = text.replace("0.025", "0.02");
    assert!(FlueDgFilePerformance::from_bytes_with_sources(misaligned.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(material.clone())).is_err());
}

#[test]
fn multimode_reed_uses_identical_material_wall_elimination_without_replacing_its_blades() {
    let blade = blade_pack();
    let text = PAIRED.replace("audio 96000 96000 4800 20", "audio 96000 96000 512 20")
        .replace("slit_length", &format!("{NORMAL}slit_length"));
    let build = |cross, block| {
        let material = normal_pack(cross);
        let source = MaterialModalDgValvePerformance::from_bytes_with_sources(text.as_bytes(), block,
            &CancelGate::new_clock_free(), None, |path| {
                if path == "wall.pack" { Ok(material.clone()) }
                else { assert_eq!(path, "blade.matpack"); Ok(blade.clone()) }
            }).unwrap();
        assert_eq!(source.wall_record().unwrap().normal_viscosity_properties(), Some(NAMES));
        assert_eq!(source.wall_record().unwrap().pack(), &material);
        assert_eq!(source.material_records()[0].pack(), &blade);
        assert_eq!(source.info().modes, 6);
        assert!(source.renderer().system().load().duct().has_viscothermal_memory());
        source
    };
    let reference = trace(build(-60.0,37).into_renderer(),512,37);
    assert_eq!(reference,trace(build(-60.0,1).into_renderer(),512,1));
    assert_ne!(reference,trace(build(60.0,37).into_renderer(),512,37));
}

#[test]
fn native_cli_snapshot_and_direct_renderer_produce_identical_normal_viscosity_pcm() {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = Scratch(std::env::temp_dir().join(format!("fs-normal-viscous-wind-{}-{stamp}",std::process::id())));
    std::fs::create_dir(&dir.0).unwrap();
    let text = flue_source(NORMAL);
    let input = dir.0.join("normal.performance");
    std::fs::write(&input,&text).unwrap();
    let material = normal_pack(-60.0);
    std::fs::write(dir.0.join("wall.pack"),material.to_bytes()).unwrap();
    let source = flue(&text,&material,37);
    let info = source.info();
    let renderer = DecimatedRenderer::new(source.into_renderer(),info.sample_rate_hz,48000,37).unwrap();
    let count = renderer.output_samples_for(info.samples).unwrap() as usize;
    let expected = encode_pcm16_wav(&trace(renderer,count,37),48000,info.full_scale_pa).unwrap().0;
    let mut receipts = Vec::new();
    for attached in [false,true] {
        let output = dir.0.join(format!("out-{attached}.wav"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_music_render"));
        command.arg("wind").arg(&input).arg(&output).args(["--block","37","--decimate"]);
        if attached { command.args(["--self-contained","--attached-source","wall.pack"]); }
        let result = command.output().unwrap();
        assert!(result.status.success(),"{}{}",String::from_utf8_lossy(&result.stdout),String::from_utf8_lossy(&result.stderr));
        assert_eq!(std::fs::read(&output).unwrap(),expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert!(metadata.contains("orthotropic-principal-local-ring-membrane"));
        assert!(metadata.contains("\"loss\":\"maxwell\",\"propagator\":\"original-exact-dg\",\"maxwell_terms\":[{"));
        assert!(metadata.contains(&material.content_hash().to_hex()));
        assert!(!metadata.contains("C_viscous=tau*C_elastic"));
        receipts.push(metadata);
    }
    assert_eq!(receipts[0],receipts[1]);
}
