//! Original material packs -> finite heat -> actual CLI pressure/PCM.
//! These are manufactured integration fixtures, never identified rosin data.
use super::*;
use crate::contact_thermal::Profile;
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_qty::Temperature;

fn thermal_pack(cp: f64, omit: Option<&str>, band: f64) -> NormalizedPack {
    let source = layer_pack(4.0, band, Some(PLASTIC_SHEAR_VISCOSITY_PROPERTY));
    let mut claims = source.claims().clone();
    let artifact = fs_blake3::ContentHash([19; 32]);
    let provenance = Provenance {
        source: "MANUFACTURED finite contact-heat fixture".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "manufactured-thermal-layer".into(),
            method: "analytic integration fixture, not measured rosin".into(),
            artifact,
            caveats: "manufactured Cp contrast exercises heat feedback".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, cp),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.2),
        (PLASTIC_SHEAR_VISCOSITY_PROPERTY, DynViscosity::DIMS, 4.0),
    ] {
        if omit == Some(name) {
            continue;
        }
        let variable = name == PLASTIC_SHEAR_VISCOSITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if variable {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(280.0, 5.315), (320.0, 1.315)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, band),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if variable {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    NormalizedPack::new(
        "manufactured-shear-layer",
        "manufactured-thermal-test-v1",
        artifact,
        "manufactured CC0 integration fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

fn thermal_spec() -> input::Spec {
    let text = include_str!("../../shear-layer-bow.material-strings")
        .replace("48000 4096 0.01 128", "48000 512 0.01 128")
        .replace("bow_ribbon 0.01 3\n", "")
        .replace("2e-2\t5e-3", "2e-3\t1e-3")
        .replace("5e-2\t5e-3", "6e-3\t1e-3");
    let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
    spec.contact_thermal = true;
    spec
}

const PROFILE: &str = "frankensim-contact-thermal-v1
boundary adiabatic
solver 256 24 1e-12 1e-13 1e-13 2e-12 1e-15 1e-12 1 1e-8
cells 1
cell 10 manufactured-shear-layer.fsmatpk 280 293.15 320 0 0 0 0.001 0 0 1e-10 2.0006e-6
films 1
film 0 10
series 0
links 0
";

fn compiled(cx: &Cx<'_>, gate: &CancelGate, cp: f64, omit: Option<&str>) -> Result<Ready> {
    thermal_spec().compile(cx, gate, |name| {
        Ok(if name == "manufactured-shear-layer.fsmatpk" {
            thermal_pack(cp, omit, 200_000.0)
        } else {
            pack(name, 7800.0, 200_000.0, None)
        })
    })
}

#[test]
fn high_rate_thermal_contact_keeps_heat_feedback_and_mechanical_trace_under_decimation() {
    scoped(|cx, gate| {
        let profile = Profile::parse(PROFILE.as_bytes()).unwrap();
        let build = |cp, block| {
            let source = String::from_utf8(thermal_spec().bytes).unwrap()
                .replace("48000 512 0.01 128", &format!("96000 1024 0.01 {block}"));
            let mut spec = input::Spec::parse(source.as_bytes()).unwrap();
            spec.contact_thermal = true;
            spec.compile(cx, gate, |name| Ok(if name == "manufactured-shear-layer.fsmatpk" {
                thermal_pack(cp, None, 400_000.0)
            } else { pack(name, 7800.0, 400_000.0, None) })).unwrap()
        };
        let mut raw = build(10.0, 128);
        let mut observed = build(10.0, 37);
        let mut cooler = build(1000.0, 128);
        let native = profile.prepare(&raw, gate).unwrap();
        let decimated = profile.prepare(&observed, gate).unwrap();
        let changed = profile.prepare(&cooler, gate).unwrap();
        native.write(&mut raw, Cursor::new(Vec::new()), None, gate).unwrap();
        let mut trace = Vec::new();
        let (wav, summary) = decimated.write_observed(&mut observed, Cursor::new(Vec::new()),
            Some(&mut trace), true, gate).unwrap();
        changed.write_observed(&mut cooler, Cursor::new(Vec::new()),
            None, true, gate).unwrap();
        assert_eq!(summary.sample_rate_hz, 48_000);
        assert_eq!(summary.samples, 512);
        assert_eq!(observed.performance.samples_rendered(), 1024);
        let text = String::from_utf8(trace).unwrap();
        assert_eq!(text.lines().count(), 1025);
        for (sample, row) in text.lines().skip(1).enumerate() {
            assert_eq!(row.split('\t').next().unwrap().parse::<u64>().unwrap(), sample as u64 + 1);
        }
        assert_eq!(raw.performance.bank().system().last_frame(), observed.performance.bank().system().last_frame());
        assert_eq!(raw.performance.bank().system().material_shear_state(0),
            observed.performance.bank().system().material_shear_state(0));
        let native_layer = raw.performance.bank().bow_shear_layer(0).unwrap();
        let warm = observed.performance.bank().bow_shear_layer(0).unwrap();
        let cool = cooler.performance.bank().bow_shear_layer(0).unwrap();
        assert_eq!(native_layer.temperature_kelvin().to_bits(), warm.temperature_kelvin().to_bits());
        assert!(warm.temperature_kelvin() > cool.temperature_kelvin());
        assert!(cool.temperature_kelvin() > 293.15);
        assert!(warm.plastic_resistance_n_s_m() < cool.plastic_resistance_n_s_m());
        assert_ne!(observed.performance.bank().system().last_frame().unwrap().observer_pressure_pa,
            cooler.performance.bank().system().last_frame().unwrap().observer_pressure_pa);
        assert_eq!(wav.into_inner().len(), 44 + 512 * 2);
    });
}

#[test]
fn cli_thermal_writer_and_trace_share_pcm_and_heat_history_and_material_cp_changes_traction() {
    scoped(|cx, gate| {
        let profile = Profile::parse(PROFILE.as_bytes()).unwrap();
        let mut plain = compiled(cx, gate, 10.0, None).unwrap();
        let mut traced = compiled(cx, gate, 10.0, None).unwrap();
        let mut high_capacity = compiled(cx, gate, 1000.0, None).unwrap();
        let physical_bytes = plain.spec.bytes.clone();
        let original = plain.performance.bank().bow_shear_layer(0).unwrap().clone();
        let a = profile.prepare(&plain, gate).unwrap();
        let b = profile.prepare(&traced, gate).unwrap();
        let c = profile.prepare(&high_capacity, gate).unwrap();
        let mut trace = Vec::new();
        let (wav, report) = a
            .write(&mut plain, Cursor::new(Vec::new()), None, gate)
            .unwrap();
        let (traced_wav, traced_report) = b
            .write(&mut traced, Cursor::new(Vec::new()), Some(&mut trace), gate)
            .unwrap();
        c.write(&mut high_capacity, Cursor::new(Vec::new()), None, gate)
            .unwrap();
        assert_eq!(wav.into_inner(), traced_wav.into_inner());
        assert_eq!(report.samples, 512);
        assert_eq!(traced_report.samples, 512);
        assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 513);
        assert_eq!(
            plain.performance.bank().system().material_shear_state(0),
            traced.performance.bank().system().material_shear_state(0)
        );
        let warm = plain.performance.bank().bow_shear_layer(0).unwrap();
        let cooler = high_capacity.performance.bank().bow_shear_layer(0).unwrap();
        assert!(warm.temperature_kelvin() > cooler.temperature_kelvin());
        assert!(cooler.temperature_kelvin() > original.temperature_kelvin());
        assert!(warm.plastic_resistance_n_s_m() < cooler.plastic_resistance_n_s_m());
        assert_eq!(physical_bytes, high_capacity.spec.bytes);
        assert_eq!(plain.catalog.len(), 3);
        assert_ne!(
            warm.resolved().card_identity(),
            cooler.resolved().card_identity()
        );
        assert_ne!(
            plain
                .performance
                .bank()
                .system()
                .last_frame()
                .unwrap()
                .observer_pressure_pa,
            high_capacity
                .performance
                .bank()
                .system()
                .last_frame()
                .unwrap()
                .observer_pressure_pa
        );
        assert_eq!(
            warm.resolved().card_identity(),
            original.resolved().card_identity()
        );
        for (current, initial) in warm
            .resolved()
            .properties()
            .iter()
            .zip(original.resolved().properties())
        {
            assert_eq!(
                current.answer().receipt.selected,
                initial.answer().receipt.selected
            );
        }
    });
}

#[test]
fn thermal_profiles_refuse_missing_sources_wrong_film_volume_and_unbound_materials() {
    for (from, to) in [
        ("boundary adiabatic", "boundary inferred-bath"),
        ("cells 1", "cells 257"),
        ("film 0 10", "film 0 99"),
        ("films 1\nfilm 0 10", "films 2\nfilm 0 10\nfilm 0 10"),
        (
            "manufactured-shear-layer.fsmatpk",
            "../manufactured-shear-layer.fsmatpk",
        ),
    ] {
        assert!(
            Profile::parse(PROFILE.replace(from, to).as_bytes()).is_err(),
            "{to}"
        );
    }
    scoped(|cx, gate| {
        assert!(compiled(cx, gate, 1000.0, Some("density")).is_err());
        let profile = Profile::parse(PROFILE.as_bytes()).unwrap();
        for absent in ["specific_heat_capacity", "thermal_conductivity"] {
            let ready = compiled(cx, gate, 1000.0, Some(absent)).unwrap();
            assert!(profile.prepare(&ready, gate).is_err(), "{absent}");
            assert_eq!(ready.performance.samples_rendered(), 0);
        }
        let ready = compiled(cx, gate, 1000.0, None).unwrap();
        for (from, to) in [
            ("1e-10 2.0006e-6", "2e-10 2.0006e-6"),
            ("280 293.15 320", "280 293.16 320"),
            ("manufactured-shear-layer.fsmatpk", "unbound.fsmatpk"),
        ] {
            let invalid = Profile::parse(PROFILE.replace(from, to).as_bytes()).unwrap();
            assert!(invalid.prepare(&ready, gate).is_err(), "{to}");
            assert_eq!(ready.performance.samples_rendered(), 0);
        }
    });
}
