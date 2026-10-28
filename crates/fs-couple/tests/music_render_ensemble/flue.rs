//! Full native flue phrases through the common physical-pressure ensemble.
use super::{directory, modal, run, succeeded};
use fs_couple::flue_loop::dg::performance::file::FlueDgFilePerformance;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::{encode_pcm16_wav, observation::{DecimatedRenderer, PressureRenderer}};
use fs_couple::render::schedule::force::file::ModalPerformance;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack,
    NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    UncertaintyModel};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

fn graph() -> String {
    let head = include_str!("../../examples/flue.performance")
        .split_once("\ngesture\n").unwrap().0.replace("output 24000 5", "output 1024 5");
    format!("{head}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t3e2\nevents\t2\nevent\t5e-3\t2e-3\tpressure-pa\t4e2\nevent\t1.5e-2\t2e-3\tpressure-pa\t0e0\n")
}

fn wall_graph() -> String {
    graph().replace("tube 0.3 0.006 1e-10 -0.9",
        "network 3 2\nduct_node inlet\nduct_node junction\nduct_node reflection -0.9\nduct_section 0 1 0.15 0.006 5e-11\nduct_section 1 2 0.15 0.006 5e-11")
        .replace("\nflue ", "\nwalls 1 20 1000\nwall 1 isotropic wall.fsmcdpk 293.15 1000 0.0002 stress-free kelvin-voigt\nradiation compact-baffled\nflue ")
        .replace("viscothermal 8 300 500",
            "viscothermal sections 2\nloss_section 0 8 300 500\nloss_section 1 8 300 500")
        .replace("\ncontrols ", "\noutlets 2 0\ncontrols ")
}

// Manufactured compliant-wall fixture; no claim of measured instrument data.
fn material(density: f64, young: f64) -> NormalizedMaterialCardPack {
    let artifact = fs_blake3::hash_domain("org.frankensim.test.flue-ensemble", b"manufactured-wall");
    let provenance = Provenance { source: "manufactured flue ensemble wall; not measured material".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact) };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic compliant cylinder".into(), method: "manufactured integration fixture".into(),
        artifact, caveats: "not a measured material or instrument".into(), provenance: provenance.clone(),
    }).unwrap();
    for (name, value, dims) in [("density", density, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS), ("poisson_ratio", 0.3, Dims::NONE),
        ("kelvin_voigt_young_viscosity", 2.0, DynViscosity::DIMS)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: fs_evidence::ValidityDomain::unconstrained().with("T", 280.0, 330.0)
                .with("omega", 0.0, 100_000.0),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated, observations: vec![observation.clone()],
            provenance: provenance.clone(),
        }).unwrap();
    }
    let pack = NormalizedPack::new("synthetic-wall", "ensemble-test-v1", artifact,
        "synthetic redistribution permitted", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: "synthetic isotropic wall".into(), phase: "solid".into(),
        process: "manufactured test state".into(), revision: 0,
    }, pack).unwrap()
}

fn trace(mut renderer: impl PressureRenderer, samples: usize) -> Vec<f64> {
    let mut pressure = vec![0.0; samples];
    for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    pressure
}

fn add_delayed_modal(pressure: &mut [f64], delay: f64) {
    assert!(delay >= 0.0 && delay.fract() == 0.0);
    let modal = ModalPerformance::from_bytes(modal(48000, pressure.len() as u64, 200.0).as_bytes(), 37).unwrap();
    let dry = trace(modal.into_renderer(), pressure.len());
    for (sum, value) in pressure[delay as usize..].iter_mut().zip(dry) { *sum += value; }
}

#[test]
fn graph_phrase_keeps_material_control_history_and_original_mechanics_clock() {
    let dir = directory();
    let source = wall_graph();
    let input = dir.join("flue.performance");
    let other = dir.join("modal.performance");
    std::fs::write(&input, &source).unwrap();
    std::fs::write(&other, modal(48000, 1024, 200.0)).unwrap();
    let mut waves = Vec::new();
    for (index, pack) in [material(1000.0, 1e6), material(2000.0, 4e6)].into_iter().enumerate() {
        // The source file and its entire phrase stay fixed; only material data changes.
        std::fs::write(dir.join("wall.fsmcdpk"), pack.to_bytes()).unwrap();
        let p = FlueFilePerformance::from_bytes_with_sources(source.as_bytes(), 37,
            &CancelGate::new_clock_free(), None, |name| {
                assert_eq!(name, "wall.fsmcdpk");
                Ok(pack.clone())
            }).unwrap();
        let info = p.info();
        assert!(info.compiled_controls > 2);
        let delay = p.renderer().receiver().filter_delay_samples();
        let mut pressure = trace(p.into_renderer(), info.samples as usize);
        add_delayed_modal(&mut pressure, delay);
        let expected = encode_pcm16_wav(&pressure, 48000, 1.0).unwrap().0;
        for block in [1, 37] {
            let output = dir.join(format!("graph-{index}-{block}.wav"));
            succeeded(run(&output, &[("--flue", &input), ("--modal", &other)], block, false));
            assert_eq!(std::fs::read(&output).unwrap(), expected);
            let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            for token in ["\"kind\":\"flue\"", "\"mechanics_samples\":2048",
                "\"mechanics_sample_rate_hz\":96000", "\"pressure_sample_rate_hz\":48000",
                "\"source_pcm_scale_applied\":false", "\"composition_authority\":\"EstimateOnly\"",
                &pack.content_hash().to_hex(), &format!("\"compiled_controls\":{}", info.compiled_controls),
                &format!("\"source_filter_delay_output_samples\":{delay}")] {
                assert!(metadata.contains(token), "{token}: {metadata}");
            }
        }
        waves.push(expected);
    }
    assert_ne!(waves[0][44..], waves[1][44..], "material data must alter mixed physical pressure");
}

#[test]
fn dg_flue_keeps_its_microphone_decimator_and_continuous_pressure_phrase() {
    let dir = directory();
    let source = include_str!("../../examples/flue-dg-baffled.performance")
        .replace("output 2048", "output 1024");
    let input = dir.join("flue-dg.performance");
    let other = dir.join("modal.performance");
    std::fs::write(&input, &source).unwrap();
    std::fs::write(&other, modal(48000, 512, 200.0)).unwrap();
    let p = FlueDgFilePerformance::from_bytes(source.as_bytes(), 37,
        &CancelGate::new_clock_free()).unwrap();
    let info = p.info();
    assert!(info.compiled_controls > 1);
    let receiver_delay = p.renderer().baffled_receiver().unwrap().filter_delay_samples() / 2.0;
    let exterior = *p.renderer().voice().dg().unwrap().load().radiation_record().unwrap().gas_record().state();
    let observed = DecimatedRenderer::new(p.into_renderer(), info.sample_rate_hz, 48000, 37).unwrap();
    let delay = receiver_delay + observed.info().delay_output_samples;
    let mut pressure = trace(observed, 512);
    add_delayed_modal(&mut pressure, delay);
    let expected = encode_pcm16_wav(&pressure, 48000, 1.0).unwrap().0;
    for block in [1, 37] {
        let output = dir.join(format!("dg-{block}.wav"));
        succeeded(run(&output, &[("--flue", &input), ("--modal", &other)], block, true));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["frankensim-flue-dg-performance-v1", "original-exact-dg",
            "\"pressure_point\":\"baffled-outlet\"", "\"mechanics_samples\":1024",
            "\"ratio\":2", "\"state_authority\":\"complete_gas_state\"",
            &format!("\"density_kg_m3\":{:e}", exterior.density),
            &format!("\"source_filter_delay_output_samples\":{receiver_delay}"),
            &format!("\"common_delay_output_samples\":{delay}")] {
            assert!(metadata.contains(token), "{token}: {metadata}");
        }
    }
    let refused = dir.join("implicit-decimation.wav");
    assert!(!run(&refused, &[("--flue", &input)], 37, false).status.success());
    assert!(!refused.exists() && !refused.with_extension("provenance.json").exists());
}

#[test]
fn explicitly_internal_dg_trace_stays_unconstrained_in_the_ensemble() {
    let dir = directory();
    let source = include_str!("../../examples/flue-dg.performance")
        .replace("output 2048", "output 1024");
    let input = dir.join("internal.performance");
    let output = dir.join("internal.wav");
    std::fs::write(&input, &source).unwrap();
    let p = FlueDgFilePerformance::from_bytes(source.as_bytes(), 37,
        &CancelGate::new_clock_free()).unwrap();
    assert!(p.renderer().baffled_receiver().is_none());
    let observed = DecimatedRenderer::new(p.into_renderer(), 96000, 48000, 37).unwrap();
    let pressure = trace(observed, 512);
    succeeded(run(&output, &[("--flue", &input)], 37, true));
    assert_eq!(std::fs::read(&output).unwrap(), encode_pcm16_wav(&pressure, 48000, 1.0).unwrap().0);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(metadata.contains("\"pressure_point\":\"mean-inlet\""));
    assert!(metadata.contains("internal DG inlet trace; no exterior microphone"));
    assert!(!metadata.contains("\"state_authority\":\"complete_gas_state\""));
}

#[test]
fn flue_source_and_listener_refusals_precede_both_output_artifacts() {
    let dir = directory();
    let input = dir.join("flue.performance");
    let other = dir.join("other.performance");
    let source = graph();
    std::fs::write(&input, &source).unwrap();
    for (index, changed) in [
        source.replace("exterior same", "exterior 323.15 101325 0"),
        source.replace("exterior same", "exterior mixture 293.15 101325 mole 2 N2 0.6 O2 0.4"),
        wall_graph(), // The explicitly named material pack is missing.
        source.replace("frankensim-flue-performance-v1", "unrecognized-performance-v1"),
    ].into_iter().enumerate() {
        std::fs::write(&other, changed).unwrap();
        let output = dir.join(format!("refused-{index}.wav"));
        let result = run(&output, &[("--flue", &input), ("--flue", &other)], 37, false);
        assert!(!result.status.success());
        if index < 2 {
            assert!(String::from_utf8_lossy(&result.stdout).contains("one ensemble has one listener gas"));
        }
        assert!(!output.exists() && !output.with_extension("provenance.json").exists());
    }
}
