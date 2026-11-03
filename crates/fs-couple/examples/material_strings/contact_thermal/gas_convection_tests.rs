//! Native profile -> original finite heat -> pressure/PCM, manufactured sources.
use super::*;
use fs_couple::air_path::input::convection::{
    ForcedGasBoundary, ForcedGasGeometry, GasConvectionBoundary,
    NaturalConvectionGeometry, NaturalGasBoundary,
};
use fs_material::gas_species::input::temperature::IsobaricGasLaw;

fn version_five(text: String, routes: &str) -> String {
    text.replace("frankensim-contact-thermal-v2", "frankensim-contact-thermal-v5")
        .replace("boundary adiabatic", "boundary gas-reservoir")
        .replace("dry-contacts 1", &format!("{routes}dry-contacts 1"))
        .replace("series 0", "film-surfaces 0\nseries 0\nmaterial-heat-budget 256 1048576 10000000\nfiber-heat 0\nstring-heat 0")
}

#[test]
fn gas_boundary_grammar_preserves_old_profiles_and_refuses_false_exposure() {
    assert!(Profile::parse(profile_text(false).as_bytes()).unwrap().convection.is_empty());
    let route = "convection 1\nconvection-cell 10 cylinder-crossflow 1e-6 0.0003 1\n";
    let text = version_five(profile_text(false), route);
    let parsed = Profile::parse(text.as_bytes()).unwrap();
    assert_eq!(parsed.convection.len(), 1);
    assert_eq!(parsed.convection[0].cell, CellId(10));
    let sphere_text = text.replace("cylinder-crossflow 1e-6", "sphere 1e-7");
    let spherical = Profile::parse(sphere_text.as_bytes()).unwrap();
    assert!(matches!(spherical.convection[0].boundary,
        GasConvectionBoundary::Forced(ForcedGasBoundary { geometry: ForcedGasGeometry::Sphere, .. })));
    assert!(Profile::parse(text.replace("cylinder-crossflow", "sphere").as_bytes()).is_err());
    for geometry in ["natural-vertical-plate", "natural-horizontal-cylinder",
        "natural-horizontal-plate-up", "natural-horizontal-plate-down", "natural-sphere"]
    {
        let natural = text.replace("cylinder-crossflow 1e-6 0.0003 1",
            &format!("{geometry} 1e-7 0.0003 9.81"));
        let parsed = Profile::parse(natural.as_bytes()).unwrap();
        assert!(matches!(parsed.convection[0].boundary, GasConvectionBoundary::Natural(_)));
        assert!(Profile::parse(natural.replace("9.81", "-9.81").as_bytes()).is_err());
    }
    assert!(Profile::parse(text.replace("cylinder-crossflow", "turbulent-flat-plate").as_bytes()).is_ok());
    for (old, new) in [
        ("convection 1", "convection 0"),
        ("convection-cell 10", "convection-cell 99"),
        ("1e-6 0.0003 1", "1 0.0003 1"),
        ("1e-6 0.0003 1", "1e-6 0.0003 0"),
        ("cylinder-crossflow", "instrument-name"),
        ("frankensim-contact-thermal-v5", "frankensim-contact-thermal-v4"),
        ("boundary gas-reservoir", "boundary adiabatic"),
    ] {
        assert!(Profile::parse(text.replace(old, new).as_bytes()).is_err(), "accepted {old} -> {new}");
    }
    let duplicate = text.replace("convection 1", "convection 2")
        .replace("dry-contacts 1", "convection-cell 10 cylinder-crossflow 1e-6 0.0003 1\ndry-contacts 1");
    assert!(Profile::parse(duplicate.as_bytes()).is_err());
}

#[test]
fn source_gas_cooling_enters_the_original_backward_euler_external_heat_once() {
    for (route, boundary) in [
        ("cylinder-crossflow 1e-6 0.0003 1",
            GasConvectionBoundary::Forced(ForcedGasBoundary {
                geometry: ForcedGasGeometry::CylinderCrossflow,
                characteristic_length_m: 0.0003, exposed_area_m2: 1e-6, speed_m_s: 1.0,
            })),
        ("natural-sphere 1e-6 0.001 9.81",
            GasConvectionBoundary::Natural(NaturalGasBoundary {
                geometry: NaturalConvectionGeometry::Sphere,
                characteristic_length_m: 0.001, exposed_area_m2: 1e-6, gravity_m_s2: 9.81,
            })),
    ] {
    scoped(|cx, gate| {
        let text = profile_text(false).replace("cells 8", "cells 9")
            .replace("dry-contacts 1", "cell 30 wire.fsmatpk 280 310 320 0 0 2 0.001 0 2 1e-9 1e-4\ndry-contacts 1");
        let profile = Profile::parse(version_five(text,
            &format!("convection 1\nconvection-cell 30 {route}\n")).as_bytes()).unwrap();
        let source = ready(cx, gate, &profile, false, 20.0, 37, 48000, None);
        let prepared = profile.prepare(&source, gate).unwrap();
        let original = prepared.network(gate).unwrap();
        let before = original.cells().clone();
        let law = IsobaricGasLaw::try_from_record(source.spec.gas_record.clone()).unwrap();
        let initial = 310.0;
        let capacity = 7800.0*1e-9*20.0;
        let dt = prepared.sample_period_s;
        let (mut lo, mut hi) = (source.spec.gas.temperature, initial);
        assert!(lo < hi);
        // Independent scalar implicit heat balance, not the native solver.
        for _ in 0..80 {
            let mid = 0.5*(lo+hi);
            let q = boundary.evaluate(&law, mid).unwrap().heat_into_wall_w;
            if capacity*(mid-initial)-dt*q > 0.0 { hi=mid; } else { lo=mid; }
        }
        let expected = 0.5*(lo+hi);
        let mut increments = vec![0.0; 9];
        increments[8] = capacity*(expected-initial);
        let step = fs_thermal::enthalpy::backward_euler::BackwardEulerStep::new(&original, dt, || Ok(())).unwrap();
        let trial = step.evaluate(&increments, |state| prepared.resolve(state), || Ok(())).unwrap();
        assert!(trial.maximum_cell_residual_j() < 1e-12);
        let (next, (), audit) = trial.accept(1e-12, 1e-12, || Ok(())).unwrap();
        assert!((next.temperature_k(CellId(30)).unwrap()-expected).abs() < 1e-9);
        assert!((audit.external_heat_j-increments[8]).abs() < 1e-12);
        assert!(audit.external_heat_j < 0.0);
        assert_eq!(next.cells()[&CellId(30)].mass_kg, original.cells()[&CellId(30)].mass_kg);
        assert_eq!(original.cells(), &before);
        assert_eq!(source.performance.samples_rendered(), 0);
        gate.request();
        assert!(prepared.resolve(&next).is_err());
        assert_eq!(original.cells(), &before);
    });
    }
}

#[test]
fn gas_boundary_changes_temperature_dependent_pressure_and_preserves_pcm_blocking() {
    scoped(|cx, gate| {
        let make_profile = |speed: f64| {
            let mut routes = String::from("convection 8\n");
            for id in [10,11,12,13,20,21,22,23] {
                routes.push_str(&format!("convection-cell {id} cylinder-crossflow 1e-6 0.0003 {speed}\n"));
            }
            Profile::parse(version_five(profile_text(false), &routes).as_bytes()).unwrap()
        };
        let mut runs = Vec::new();
        for speed in [1.0, 4.0] {
            let profile = make_profile(speed);
            // Only boundary speed changes: material, geometry, excitation and gas stay fixed.
            let mut source = temperature_ready(cx, gate, &profile, 0.02, 0.005, 37);
            let prepared = profile.prepare(&source, gate).unwrap();
            let routes = prepared.dry_routes(&source).unwrap();
            let mut renderer = source.performance.with_dry_contact_thermal_network(
                prepared.network(gate).unwrap(), routes, vec![], profile.controls,
                gate, |endpoint| prepared.resolve(endpoint),
            ).unwrap();
            let mut pressure = vec![0.0; 128];
            for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
            let temperature = renderer.last_frame().unwrap().dry_contacts[0]
                .interface_temperature_k.unwrap();
            runs.push((pressure, temperature));
        }
        assert!((runs[0].1-runs[1].1).abs() > 1e-10);
        assert!(runs[0].0.iter().zip(&runs[1].0).any(|(a,b)| a.to_bits()!=b.to_bits()));
        let profile = make_profile(1.0);
        let mut source = temperature_ready(cx, gate, &profile, 0.02, 0.005, 128);
        let prepared = profile.prepare(&source, gate).unwrap();
        let mut trace = Vec::new();
        let (wav, _) = prepared.write(&mut source, Cursor::new(Vec::new()), Some(&mut trace), gate).unwrap();
        assert_eq!(wav.into_inner(), fs_couple::pcm_wav::encode_pcm16_wav(
            &runs[0].0, 48000, 0.01).unwrap().0);
        assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 129);
    });
}
