//! Full source tensors change physical bridge reactions and emitted pressure.
use super::*;

fn tensor_pack(c33: f64, time: f64, limit: f64, missing: &str, upper: f64) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 900.0, Density::DIMS),
        ("young_modulus", 4e9, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE),
        ("strain_limit", limit, Dims::NONE),
        ("relaxation_time", time, Dims([0, 0, 1, 0, 0, 0])),
        ("branch_11", 1e8, Pressure::DIMS),
        ("branch_22", 2e8, Pressure::DIMS),
        ("branch_33", c33, Pressure::DIMS),
        ("branch_44", 1.5e8, Pressure::DIMS),
        ("branch_55", 1.2e8, Pressure::DIMS),
        ("branch_66", 0.8e8, Pressure::DIMS),
        ("branch_13", 4e7, Pressure::DIMS),
        ("branch_24", 2e7, Pressure::DIMS),
        ("branch_15", -1e7, Pressure::DIMS),
        ("zero_modulus", 0.0, Pressure::DIMS),
    ] {
        if name == missing {
            continue;
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, upper),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured complete tensor body".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    NormalizedPack::new(
        "tensor-body",
        "synthetic-v1",
        fs_blake3::ContentHash([43; 32]),
        "CC0-1.0",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}
fn tensor_spec() -> input::Spec {
    let mut spec = input::Spec::parse(include_bytes!("../tensor-body.material-strings")).unwrap();
    spec.samples = 256;
    spec
}
fn tensor_profile() -> Profile {
    Profile::parse(include_bytes!("../tensor-body.body")).unwrap()
}
fn build(cx: &Cx<'_>, gate: &CancelGate, c33: f64, time: f64) -> Ready {
    tensor_spec()
        .compile_with_body(cx, gate, Some(tensor_profile()), |name| {
            Ok(if name == "tensor-blade.fsmatpk" {
                tensor_pack(c33, time, 0.02, "", 100000.0)
            } else {
                pack(name, 1.0, None)
            })
        })
        .unwrap()
}

#[test]
fn g3_complete_tensor_body_source_changes_normal_history_bridge_feedback_and_pcm() {
    scoped(|cx, gate| {
        let mut original = build(cx, gate, 3e8, 0.001);
        let source = original.performance.bank().source().bank();
        let image = source.material_chart_tensor_relaxation().unwrap();
        assert!(source.material_chart_viscosity().is_none());
        assert!(source.material_chart_relaxation().is_none());
        assert_eq!(image.mode_count(), 2);
        assert_eq!(image.material_points().len(), 6);
        assert_eq!(image.coupling_rows_per_s().len(), 36);
        assert_eq!(
            original.performance.bank().system().material_memory().len(),
            36
        );
        let modes = original.performance.bank().system().components()[1]
            .modes()
            .to_vec();
        let p = pressure(&mut original, 37);
        assert!(p.iter().any(|p| p.abs() > 1e-15));
        for (c33, time) in [(8e9, 0.001), (3e8, 0.005)] {
            let mut changed = build(cx, gate, c33, time);
            assert_eq!(
                changed.performance.bank().system().components()[1].modes(),
                modes
            );
            let other = pressure(&mut changed, 128);
            assert!(p.iter().zip(&other).any(|(a, b)| (a - b).abs() > 1e-15));
            assert_ne!(
                original.performance.bank().system().components()[0].states(),
                changed.performance.bank().system().components()[0].states()
            );
        }
        let mut replay = build(cx, gate, 3e8, 0.001);
        assert_eq!(p, pressure(&mut replay, 1));
        let expected =
            fs_couple::pcm_wav::encode_pcm16_wav(&p, original.spec.rate, original.spec.full_scale)
                .unwrap()
                .0;
        let (wav, summary) = build(cx, gate, 3e8, 0.001)
            .write(std::io::Cursor::new(Vec::new()), gate)
            .unwrap();
        assert_eq!(wav.into_inner(), expected);
        assert_eq!(summary.samples, 256);
        let system = original.performance.bank().system();
        let q: Vec<_> = system.components()[1]
            .states()
            .iter()
            .map(|s| 3.0_f64.sqrt() * s.displacement_m_sqrt_kg)
            .collect();
        for point in original
            .performance
            .bank()
            .source()
            .bank()
            .material_chart_tensor_relaxation()
            .unwrap()
            .material_points()
        {
            let memory: Vec<_> = system.material_memory()[point.memory_range()]
                .iter()
                .map(|z| 3.0_f64.sqrt() * z)
                .collect();
            assert!(point.image().strain(&q, &memory).is_ok());
            assert!(
                point
                    .image()
                    .normal_stress_pa(&q, &memory)
                    .unwrap()
                    .iter()
                    .all(|s| s.abs() < 1e-5)
            );
        }
        let frame = system.last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn g0_native_tensor_body_preserves_complete_source_strain_band_and_capacity_admission() {
    scoped(|cx, gate| {
        for (missing, upper, time, limit) in [
            ("branch_33", 100000.0, 0.001, 0.02),
            ("zero_modulus", 100000.0, 0.001, 0.02),
            ("", 100.0, 0.001, 0.02),
            ("", 100000.0, 1e-9, 0.02),
            ("", 100000.0, 0.001, 1e-20),
        ] {
            assert!(
                tensor_spec()
                    .compile_with_body(cx, gate, Some(tensor_profile()), |name| Ok(
                        if name == "tensor-blade.fsmatpk" {
                            tensor_pack(3e8, time, limit, missing, upper)
                        } else {
                            pack(name, 1.0, None)
                        }
                    ))
                    .is_err()
            );
        }
        for replacement in ["1 6 35 67108864", "1 5 36 67108864", "1 6 36 1"] {
            let text = include_str!("../tensor-body.body").replace("1 6 36 67108864", replacement);
            assert!(
                tensor_spec()
                    .compile_with_body(
                        cx,
                        gate,
                        Some(Profile::parse(text.as_bytes()).unwrap()),
                        |name| Ok(if name == "tensor-blade.fsmatpk" {
                            tensor_pack(3e8, 0.001, 0.02, "", 100000.0)
                        } else {
                            pack(name, 1.0, None)
                        })
                    )
                    .is_err()
            );
        }
        // Orthotropic source directions 1/2 do not authorize inventing axis 3.
        let orthotropic = include_str!("../tensor-body.body")
            .replace("elastic isotropic", "elastic orthotropic12");
        assert!(
            tensor_spec()
                .compile_with_body(
                    cx,
                    gate,
                    Some(Profile::parse(orthotropic.as_bytes()).unwrap()),
                    |_| Ok(pack("body.fsmatpk", 1.0, None))
                )
                .is_err()
        );
    });
}

#[test]
fn g0_cantilever_geometry_and_tensor_grammar_never_invent_missing_rows() {
    let source = include_str!("../tensor-body.body");
    assert!(Profile::parse(source.as_bytes()).is_ok());
    for text in [
        source.replace("clamped-x-min", "clamped"),
        source.replace("bridge 1 1", "bridge 0 1"),
        source.replace("bridge 1 1", "bridge 2 1"),
        source.replace(" 36 67108864", " 257 67108864"),
        source.replace(" 1 6 36", " 9 6 36"),
        source.replace("tensor_branch relaxation_time", "tensor_branch"),
        source.replace(" branch_66\n", "\n"),
        format!("{source}ignored\n"),
    ] {
        assert!(Profile::parse(text.as_bytes()).is_err());
    }
}
