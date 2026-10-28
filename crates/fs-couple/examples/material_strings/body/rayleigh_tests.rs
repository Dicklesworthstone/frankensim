//! Source data -> full modal radiation -> original performer and PCM.
use super::*;

fn selected(pressure_pa: f64) -> input::Spec {
    let text = include_str!("../rayleigh-pluck.material-strings")
        .replace("directional-body.fsmatpk", "body.fsmatpk")
        .replace("293.15 101325 mole", &format!("293.15 {pressure_pa} mole"));
    let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
    spec.samples = 48;
    spec
}

fn ready(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, profile: Profile) -> Ready {
    spec.compile_with_body(cx, gate, Some(profile), |name| Ok(pack(name, 1.0, None))).unwrap()
}

fn report(ready: &Ready) -> &PlateRayleighReport {
    let BodyRadiationReport::Rayleigh(report) = ready.body_radiation.as_ref().unwrap() else {
        panic!("Rayleigh source selected a different radiation model");
    };
    report
}

#[test]
fn computed_air_load_changes_motion_and_streams_the_same_accepted_pressure() {
    scoped(|cx, gate| {
        let spec = selected(101_325.0);
        let mut loaded = ready(cx, gate, spec.clone(), profile());
        let report = report(&loaded);
        assert_eq!(report.body_modes, 1);
        assert_eq!(report.surface_cells, 96);
        assert!(report.acoustic_poles > 0);
        assert!(report.peak_error <= 0.05 && report.rms_error <= 0.02);
        assert_eq!(report.gas, spec.gas);
        let mut bare_spec = spec.clone();
        bare_spec.body_radiation = None;
        let mut bare = ready(cx, gate, bare_spec, profile());
        assert_eq!(loaded.performance.bank().system().mode_count(), bare.performance.bank().system().mode_count());
        let p = pressure(&mut loaded, 7);
        let q = pressure(&mut bare, 48);
        assert!(p.iter().any(|v| v.abs() > 1e-12));
        assert!(p.iter().zip(&q).any(|(a,b)| (a-b).abs() > 1e-15), "air never reacted on the instrument");
        assert_ne!(loaded.performance.bank().system().components()[0].states(),
            bare.performance.bank().system().components()[0].states(), "radiation never reached the string");
        let mut streamed = ready(cx, gate, spec, profile());
        let encoded = fs_couple::pcm_wav::encode_pcm16_wav(&p, streamed.spec.rate, streamed.spec.full_scale).unwrap().0;
        let (wav, summary) = streamed.write(std::io::Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), encoded);
        assert_eq!(summary.samples, 48);
        for image in [&loaded, &bare, &streamed] {
            let frame = image.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn gas_and_grain_data_rebuild_the_actual_same_basis_radiation_matrix() {
    scoped(|cx, gate| {
        let a = ready(cx, gate, selected(101_325.0), profile());
        let b = ready(cx, gate, selected(202_650.0), profile());
        let turned = Profile::parse(include_str!("../orthotropic.body")
            .replace("grain 0.17", "grain 0.87").as_bytes()).unwrap();
        let c = ready(cx, gate, selected(101_325.0), turned);
        let ar = report(&a);
        let br = report(&b);
        let cr = report(&c);
        assert_eq!(ar.specimen_identity, br.specimen_identity);
        assert_ne!(ar.specimen_identity, cr.specimen_identity);
        assert_eq!(ar.gas.sound_speed, br.gas.sound_speed);
        assert_eq!(2.0*ar.gas.density, br.gas.density);
        let impedance = |ready: &Ready| ready.performance.bank().source().bank()
            .body_radiation().unwrap().model().impedance(core::f64::consts::TAU*100.0).unwrap();
        let za = impedance(&a); let zb = impedance(&b); let zc = impedance(&c);
        assert!((zb[0]-za[0].scale(2.0)).abs() <= 1e-8 * za[0].abs().max(1e-12));
        assert!((zc[0]-za[0]).abs() > 1e-10 * za[0].abs().max(1e-12));
    });
}

#[test]
fn incompatible_or_unbounded_radiation_is_refused_before_pack_io() {
    scoped(|cx, gate| {
        let mut calls = 0;
        assert!(selected(101_325.0).compile_with_body(cx, gate, None, |name| {
            calls += 1; Ok(pack(name, 1.0, None))
        }).is_err());
        assert_eq!(calls, 0);
        let text = include_str!("../rayleigh-pluck.material-strings");
        for selected in ["rayleigh 1 5 300 64 256 10000000 48",
            "rayleigh 1 5 3000 65 256 10000000 48",
            "rayleigh 1 5 300 65 999999 10000000 48",
            "rayleigh 3 5 300 65 256 10000000 48"] {
            assert!(input::Spec::parse(text.replace("rayleigh 1 5 300 65 256 10000000 48", selected).as_bytes()).is_err());
        }
    });
}

#[test]
fn bem_source_grammar_retains_its_complete_physical_resolution() {
    let original = include_str!("../rayleigh-pluck.material-strings");
    let selected = original.replace("rayleigh 1 5 300 65 256 10000000 48",
        "bem 10 1000 65 80 100000000 64");
    let spec = input::Spec::parse(selected.as_bytes()).unwrap();
    let Some(input::BodyRadiation::ShellBem(options)) = spec.body_radiation else {
        panic!("BEM source was silently converted to a different acoustic problem");
    };
    assert_eq!(options.minimum_hz, 10.0);
    assert_eq!(options.maximum_hz, 1000.0);
    assert_eq!(options.frequency_samples, 65);
    assert_eq!(options.maximum_surface_panels, 80);
    assert_eq!(options.maximum_sweep_terms, 100_000_000);
    assert_eq!(options.maximum_poles, 64);
    for invalid in [
        "bem 10 1000 64 80 100000000 64",
        "bem 10 1000 65 3 100000000 64",
        "bem 10 1000 65 2049 100000000 64",
        "bem 10 1000 65 80 0 64",
        "bem 10 1000 65 80 100000000 0",
        "bem 10 1000 65 80 100000000 64 ignored",
    ] {
        assert!(input::Spec::parse(original.replace(
            "rayleigh 1 5 300 65 256 10000000 48", invalid).as_bytes()).is_err(),
            "accepted {invalid}");
    }
}

#[test]
fn bem_on_a_missing_or_planar_profile_refuses_before_material_io() {
    scoped(|cx, gate| {
        for body in [None, Some(profile())] {
            let selected = include_str!("../rayleigh-pluck.material-strings")
                .replace("rayleigh 1 5 300 65 256 10000000 48",
                    "bem 10 1000 65 80 100000000 64");
            let spec = input::Spec::parse(selected.as_bytes()).unwrap();
            let mut loads = 0;
            let result = spec.compile_with_body(cx, gate, body, |name| {
                loads += 1;
                Ok(pack(name, 1.0, None))
            });
            assert!(result.is_err());
            assert_eq!(loads, 0);
        }
    });
}

fn observed() -> input::Spec {
    let text = include_str!("../rayleigh-observed.material-strings")
        .replace("directional-body.fsmatpk", "body.fsmatpk");
    let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
    spec.samples = 1024;
    spec
}

#[test]
fn material_plate_observer_keeps_loaded_motion_and_streams_physical_pcm() {
    scoped(|cx, gate| {
        let spec = observed();
        let mut observed = ready(cx, gate, spec.clone(), profile());
        let report = observed.performance.material_plate_observer().unwrap();
        assert_eq!(report.surface_cells, 96);
        assert_eq!(report.body_modes, 1);
        assert_eq!(report.gas.gas, spec.gas);
        assert_eq!(report.gas.surface_filter_indices.len(), report.surface_cells);
        assert!(report.gas.filters.len() > 1, "surface paths were collapsed to one distance");
        assert!(report.filter_delay_samples > 0.0);
        assert!(report.replay_terms <= report.options.maximum_replay_terms);
        let mut compact_spec = spec.clone();
        compact_spec.plate_observer = None;
        let mut compact = ready(cx, gate, compact_spec, profile());
        let p = pressure(&mut observed, 37);
        let q = pressure(&mut compact, 128);
        assert!(p.iter().any(|v| v.abs() > 1e-12));
        assert_ne!(p, q, "the finite-distance receiver left compact pressure in place");
        for (a, b) in observed.performance.bank().system().components().iter()
            .zip(compact.performance.bank().system().components()) {
            assert_eq!(a.states(), b.states(), "an output observer changed accepted mechanics");
        }
        let frame = observed.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let mut streamed = ready(cx, gate, spec, profile());
        let expected = fs_couple::pcm_wav::encode_pcm16_wav(&p, streamed.spec.rate, streamed.spec.full_scale).unwrap().0;
        let (wav, summary) = streamed.write(std::io::Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), expected);
        assert_eq!(summary.samples, 1024);
    });
}

#[test]
fn material_plate_observer_uses_transport_and_preserves_both_mechanical_histories() {
    scoped(|cx, gate| {
        let a = observed();
        let mut b = a.clone();
        b.plate_observer.as_mut().unwrap().propagation.law = fs_couple::air_path::GasAbsorptionModel::Lossless;
        let mut classical = ready(cx, gate, a, profile());
        let mut lossless = ready(cx, gate, b, profile());
        let p = pressure(&mut classical, 37);
        let q = pressure(&mut lossless, 64);
        let scale = p.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!(scale > 0.0);
        assert!(p.iter().zip(q).any(|(a, b)| (a-b).abs() > 1e-10*scale),
            "the complete material gas did not reach actual pressure");
        for (a, b) in classical.performance.bank().system().components().iter()
            .zip(lossless.performance.bank().system().components()) {
            assert_eq!(a.states(), b.states(), "propagation changed physical reaction or motion");
        }
    });
}

#[test]
fn material_plate_observer_refuses_incomplete_or_unbounded_native_data() {
    let text = include_str!("../rayleigh-observed.material-strings");
    for (from, to) in [
        ("plate_observer 0.7 0.1 1", "plate_observer 0.7 0.1 0.01"),
        ("300 256 100000", "300 2 100000"),
        ("300 256 100000", "300 256 0"),
        ("stokes-kirchhoff 256", "lossless 0"),
        ("stokes-kirchhoff 256", "molecular-relaxation 256"),
        ("body-radiation rayleigh 1 5 300 65 256 10000000 48\n", ""),
        ("cavity none", "cavity neck 0.01 0.015 0.008 1000 1"),
    ] {
        assert!(input::Spec::parse(text.replace(from, to).as_bytes()).is_err(), "accepted {to}");
    }
    scoped(|cx, gate| {
        let mut spec = observed();
        spec.plate_observer.as_mut().unwrap().maximum_surface_cells = 95;
        assert!(spec.compile_with_body(cx, gate, Some(profile()), |name| Ok(pack(name, 1.0, None))).is_err(),
            "the receiver silently dropped a surface quadrature point");
        let mut spec = observed();
        spec.plate_observer.as_mut().unwrap().maximum_replay_terms = 1;
        assert!(spec.compile_with_body(cx, gate, Some(profile()), |name| Ok(pack(name, 1.0, None))).is_err(),
            "the receiver ignored its declared runtime work");
    });
}

#[test]
fn material_plate_observer_molecular_data_binds_the_original_gas_species() {
    let text = include_str!("../rayleigh-observed.material-strings").replace(
        "observer_gas stokes-kirchhoff 256 1000000 50000000",
        "observer_gas molecular 256 1000000 50000000 declared 2 N2 1 vibration 1 3390 1 power 0.02 293.15 0 280 320 O2 1 vibration 1 2270 1 power 0.01 293.15 0 280 320",
    );
    let spec = input::Spec::parse(text.as_bytes()).unwrap();
    let propagation = &spec.plate_observer.as_ref().unwrap().propagation;
    assert_eq!(propagation.law, fs_couple::air_path::GasAbsorptionModel::MolecularRelaxation);
    let relaxation = propagation.relaxation.as_ref().unwrap();
    assert_eq!(relaxation.gas(), &spec.gas);
    assert_eq!(relaxation.processes().len(), 2);
    assert!(input::Spec::parse(text.replace("N2 0.8 O2 0.2", "He 0.8 O2 0.2").as_bytes()).is_err(),
        "a new atmosphere silently borrowed another species' molecular data");
}
