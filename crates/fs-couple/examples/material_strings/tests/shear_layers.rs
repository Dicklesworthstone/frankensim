//! Exercise pack -> layer geometry -> original coupled pressure through the CLI.
use super::super::*;
use super::{pack, pressure, scoped};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    PLASTIC_SHEAR_VISCOSITY_PROPERTY, SHEAR_LAYER_MODULUS_PROPERTY, SHEAR_LAYER_YIELD_PROPERTY,
};
use std::io::Cursor;

fn layer_pack(viscosity: f64, band: f64, omit: Option<&str>) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let artifact = fs_blake3::ContentHash([17; 32]);
    let provenance = Provenance {
        source: "MANUFACTURED isothermal plastic shear fixture".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "manufactured-shear-layer".into(),
            method: "analytic integration fixture, not measured rosin".into(),
            artifact,
            caveats: "plastic viscosity is not oscillatory viscosity".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, dims, value) in [
        (SHEAR_LAYER_MODULUS_PROPERTY, Pressure::DIMS, 100_000.0),
        (SHEAR_LAYER_YIELD_PROPERTY, Pressure::DIMS, 100.0),
        (
            PLASTIC_SHEAR_VISCOSITY_PROPERTY,
            DynViscosity::DIMS,
            viscosity,
        ),
    ] {
        if omit == Some(name) {
            continue;
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, band),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let source = NormalizedPack::new(
        "manufactured-shear-layer",
        "manufactured-test-v1",
        artifact,
        "manufactured CC0 integration fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedPack::from_bytes(&source.to_bytes()).unwrap()
}
fn specification() -> input::Spec {
    // Include both velocity reversals in a short 1024-sample run.
    let text = include_str!("../shear-layer-bow.material-strings")
        .replace("48000 4096 0.01 128", "48000 1024 0.01 128")
        .replace("2e-2\t5e-3", "5e-3\t2e-3")
        .replace("5e-2\t5e-3", "1.2e-2\t2e-3");
    input::Spec::parse(text.as_bytes()).unwrap()
}
fn ready(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, viscosity: f64) -> Ready {
    spec.compile(cx, gate, |name| {
        Ok(if name == "manufactured-shear-layer.fsmatpk" {
            layer_pack(viscosity, 200_000.0, None)
        } else {
            pack(name, 7800.0, 200_000.0, None)
        })
    })
    .unwrap()
}

#[test]
fn source_only_layer_swap_reaches_pressure_with_one_catalog_and_one_ribbon_area() {
    scoped(|cx, gate| {
        let mut loads = BTreeMap::<String, usize>::new();
        let mut a = specification()
            .compile(cx, gate, |name| {
                *loads.entry(name.into()).or_default() += 1;
                Ok(if name == "manufactured-shear-layer.fsmatpk" {
                    layer_pack(4.0, 200_000.0, None)
                } else {
                    pack(name, 7800.0, 200_000.0, None)
                })
            })
            .unwrap();
        let mut b = ready(cx, gate, specification(), 12.0);
        assert_eq!(a.spec.bytes, b.spec.bytes);
        assert_eq!(a.catalog.len(), 3);
        assert!(loads.values().all(|count| *count == 1));
        assert_ne!(
            a.catalog["manufactured-shear-layer.fsmatpk"]
                .pack
                .content_hash(),
            b.catalog["manufactured-shear-layer.fsmatpk"]
                .pack
                .content_hash()
        );
        let mut total_area = 0.0;
        let mut total_stiffness = 0.0;
        let mut total_yield = 0.0;
        let mut total_resistance = 0.0;
        for port in 0..3 {
            let layer = a.performance.bank().bow_shear_layer(port).unwrap();
            total_area += layer.area_m2();
            total_stiffness += layer.stiffness_n_m();
            total_yield += layer.yield_force_n();
            total_resistance += layer.plastic_resistance_n_s_m();
            assert_eq!(layer.thickness_m(), 1e-4);
            assert_eq!(layer.temperature_kelvin(), 293.15);
        }
        for (actual, expected) in [
            (total_area, 1e-6),
            (total_stiffness, 1000.0),
            (total_yield, 1e-4),
            (total_resistance, 0.04),
        ] {
            assert!((actual / expected - 1.0).abs() < 1e-14);
        }
        assert!(a.performance.bank().bow_shear_layer(3).is_none());
        let x = pressure(&mut a, 37);
        let y = pressure(&mut b, 128);
        assert!(x.iter().any(|sample| sample.abs() > 1e-12));
        assert!(
            x.iter()
                .zip(y)
                .any(|(left, right)| (left - right).abs() > 1e-12)
        );
        let report = a.performance.bank().system().last_frame().unwrap();
        assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
        assert!(report.ports.iter().all(|port| port.shear_layer.is_some()));
    });
}
#[test]
fn layer_history_and_pcm_are_identical_across_callback_boundaries() {
    scoped(|cx, gate| {
        let mut a = ready(cx, gate, specification(), 4.0);
        let mut b = ready(cx, gate, specification(), 4.0);
        let x = pressure(&mut a, 1);
        assert_eq!(x, pressure(&mut b, 37));
        for port in 0..3 {
            assert_eq!(
                a.performance.bank().system().material_shear_state(port),
                b.performance.bank().system().material_shear_state(port)
            );
        }
        let (expected, clipped) = fs_couple::pcm_wav::encode_pcm16_wav(&x, 48_000, 0.01).unwrap();
        let mut output = ready(cx, gate, specification(), 4.0);
        let (wav, summary) = output.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), expected);
        assert_eq!(summary.clipped_samples, clipped as u64);
        assert_eq!(summary.samples, 1024);
    });
}

#[test]
fn layer_material_swaps_reach_nonlinear_string_pressure_without_resetting_interface_history() {
    scoped(|cx, gate| {
        for record in ["geometric-stretching 0.2 64 100000 1e-11 1e-10",
            "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2"] {
            let text = include_str!("../shear-layer-bow.material-strings")
                .replace("48000 4096 0.01 128", "48000 256 0.01 128")
                .replace("cavity neck", &format!("{record}\ncavity neck"));
            let spec = || input::Spec::parse(text.as_bytes()).unwrap();
            let mut direct = ready(cx, gate, spec(), 4.0);
            let mut replay = ready(cx, gate, spec(), 4.0);
            let mut changed = ready(cx, gate, spec(), 12.0);
            assert_eq!(direct.spec.bytes, changed.spec.bytes);
            assert_eq!(direct.performance.bank().system().components().iter().map(|m| m.modes()).collect::<Vec<_>>(),
                changed.performance.bank().system().components().iter().map(|m| m.modes()).collect::<Vec<_>>());
            let expected = pressure(&mut direct, 1);
            assert_eq!(expected, pressure(&mut replay, 37));
            let replacement = pressure(&mut changed, 128);
            assert!(expected.iter().any(|p| p.abs() > 1e-12));
            assert!(expected.iter().zip(replacement).any(|(a, b)| (a - b).abs() > 1e-12));
            for port in 0..3 {
                assert_eq!(direct.performance.bank().system().material_shear_state(port),
                    replay.performance.bank().system().material_shear_state(port));
            }
            let frame = direct.performance.bank().system().last_frame().unwrap();
            assert!(frame.nonlinear.is_some());
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert!(frame.ports.iter().all(|port| port.shear_layer.is_some()));
            assert!(frame.shear_layer_energy_j > 0.0);
            let mut output = ready(cx, gate, spec(), 4.0);
            let (wav, report) = output.write(Cursor::new(Vec::new()), gate).unwrap();
            assert_eq!(wav.into_inner(), fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48_000, 0.01).unwrap().0);
            assert_eq!(report.samples, 256);

            // The thermal endpoint owner remains a separate admitted rung.
            let mut unsupported = spec(); unsupported.contact_thermal = true;
            let mut reads = 0;
            assert!(unsupported.compile(cx, gate, |name| {
                reads += 1; Ok(pack(name, 7800.0, 200_000.0, None))
            }).err().unwrap().to_string().contains("contact-thermal"));
            assert_eq!(reads, 0);
        }
    });
}
#[test]
fn incomplete_layer_declarations_sources_bands_and_opening_do_not_fall_back_to_dry_friction() {
    let text = include_str!("../shear-layer-bow.material-strings");
    for (from, to) in [
        ("friction layers 1", "friction layers 0"),
        ("shear-layer bow ", "shear-layer absent "),
        ("0.000001 0.0001 0.05", "0.000001 0.0001 0.2"),
        (
            "manufactured-shear-layer.fsmatpk",
            "../manufactured-shear-layer.fsmatpk",
        ),
    ] {
        assert!(
            input::Spec::parse(text.replace(from, to).as_bytes()).is_err(),
            "accepted {to}"
        );
    }
    scoped(|cx, gate| {
        for (band, omit) in [
            (200_000.0, Some(PLASTIC_SHEAR_VISCOSITY_PROPERTY)),
            (1000.0, None),
        ] {
            assert!(
                specification()
                    .compile(cx, gate, |name| {
                        Ok(if name == "manufactured-shear-layer.fsmatpk" {
                            layer_pack(4.0, band, omit)
                        } else {
                            pack(name, 7800.0, 200_000.0, None)
                        })
                    })
                    .is_err()
            );
        }
        let opening = text.replace("initial\tbow\t1e-2\t1e-2", "initial\tbow\t1e-2\t0e0");
        let spec = input::Spec::parse(opening.as_bytes()).unwrap();
        assert!(
            spec.compile(cx, gate, |name| {
                Ok(if name == "manufactured-shear-layer.fsmatpk" {
                    layer_pack(4.0, 200_000.0, None)
                } else {
                    pack(name, 7800.0, 200_000.0, None)
                })
            })
            .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        let mut calls = 0;
        assert!(
            specification()
                .compile(cx, &cancelled, |name| {
                    calls += 1;
                    Ok(pack(name, 7800.0, 200_000.0, None))
                })
                .is_err()
        );
        assert_eq!(calls, 0);
    });
}

#[path = "shear_layers/contact_thermal.rs"]
mod contact_thermal;
