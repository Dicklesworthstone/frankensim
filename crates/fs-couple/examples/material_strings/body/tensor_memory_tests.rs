//! Full source tensors change physical bridge reactions and emitted pressure.
use super::*;
use fs_material::visco_tensor::source::TensorElasticity;

fn tensor_pack(c33: f64, time: f64, limit: f64, missing: &str, upper: f64) -> NormalizedPack {
    tensor_pack_with_elasticity(c33, time, limit, missing, upper, TensorElasticity::Isotropic)
}

fn tensor_pack_with_elasticity(
    c33: f64, time: f64, limit: f64, missing: &str, upper: f64, elasticity: TensorElasticity,
) -> NormalizedPack {
    let mut values = vec![
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
    ];
    match elasticity {
        TensorElasticity::Isotropic => {}
        TensorElasticity::IsotropicYoungShear | TensorElasticity::IsotropicBulkShear => {
            values.retain(|(name, _, _)| *name != "poisson_ratio"
                && (elasticity != TensorElasticity::IsotropicBulkShear || *name != "young_modulus"));
            values.push(("shear_modulus", 4e9 / (2.0 * 1.3), Pressure::DIMS));
            if elasticity == TensorElasticity::IsotropicBulkShear {
                values.push(("bulk_modulus", 4e9 / (3.0 * (1.0 - 2.0 * 0.3)), Pressure::DIMS));
            }
        }
        TensorElasticity::Orthotropic => unreachable!("isotropic body fixture"),
    }
    let raw = format!(
        "manufactured tensor-body fixture; T=[280,320] K; omega=[0,{upper:.17e}] rad/s\n{}",
        values.iter().filter(|(name, _, _)| *name != missing)
            .map(|(name, value, dims)| format!("{name}\t{value:.17e}\t{dims:?}\n"))
            .collect::<String>()
    );
    let artifact = fs_blake3::hash_bytes(raw.as_bytes());
    let provenance = Provenance {
        source: "manufactured complete tensor body".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic tensor-body integration specimen".into(),
        method: "explicit manufactured coefficients; no measurement or fit".into(),
        artifact,
        caveats: format!("Not measured material data. Retained raw fixture:\n{raw}"),
        provenance: provenance.clone(),
    }).unwrap();
    for (name, value, dims) in values {
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
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    NormalizedPack::new(
        "tensor-body",
        "synthetic-v1",
        artifact,
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
    build_profile(cx,gate,c33,time,tensor_profile())
}
fn build_profile(cx:&Cx<'_>,gate:&CancelGate,c33:f64,time:f64,profile:Profile)->Ready {
    tensor_spec()
        .compile_with_body(cx, gate, Some(profile), |name| {
            Ok(if name == "tensor-blade.fsmatpk" {
                tensor_pack(c33, time, 0.02, "", 100000.0)
            } else {
                pack(name, 1.0, None)
            })
        })
        .unwrap()
}

#[test]
fn g1_independent_equilibrium_pairs_reach_native_tensor_body_feedback_and_pressure() {
    scoped(|cx, gate| {
        let mut reference = build(cx, gate, 3e8, 0.001);
        let expected = pressure(&mut reference, 37);
        let scale = expected.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!(scale > 1e-15);
        for (token, elasticity, model) in [
            ("isotropic-young-shear", TensorElasticity::IsotropicYoungShear,
                PlateMaterialModel::IsotropicYoungShear),
            ("isotropic-bulk-shear", TensorElasticity::IsotropicBulkShear,
                PlateMaterialModel::IsotropicBulkShear),
        ] {
            let profile = include_str!("../tensor-body.body")
                .replace("elastic isotropic", &format!("elastic {token}"));
            let mut ready = tensor_spec().compile_with_body(
                cx, gate, Some(Profile::parse(profile.as_bytes()).unwrap()), |name| {
                    Ok(if name == "tensor-blade.fsmatpk" {
                        tensor_pack_with_elasticity(3e8, 0.001, 0.02, "", 100000.0, elasticity)
                    } else {
                        pack(name, 1.0, None)
                    })
                },
            ).unwrap();
            let source = ready.performance.bank().source().bank();
            let chart = source.material_chart().unwrap();
            for region in chart.specimen.regions() {
                assert_eq!(region.input.model, model);
                assert!(region.input.material.property("poisson_ratio").is_none());
                assert!(region.input.material.property("shear_modulus").is_some());
                if elasticity == TensorElasticity::IsotropicBulkShear {
                    assert!(region.input.material.property("young_modulus").is_none());
                    assert!(region.input.material.property("bulk_modulus").is_some());
                }
            }
            let history = source.material_chart_tensor_relaxation().unwrap();
            assert_eq!(history.mode_count(), 2);
            assert_eq!(history.material_points().len(), 6);
            assert_eq!(history.coupling_rows_per_s().len(), 36);
            assert_eq!(ready.performance.bank().system().material_memory().len(), 36);
            let actual = pressure(&mut ready, 128);
            assert!(expected.iter().zip(&actual).all(|(a, b)| (a - b).abs() <= 1e-8 * scale));
            let before = reference.performance.bank().system().last_frame().unwrap();
            let after = ready.performance.bank().system().last_frame().unwrap();
            for (a, b) in [
                (before.network_energy_j, after.network_energy_j),
                (before.network_dissipation_j, after.network_dissipation_j),
            ] {
                assert!((a - b).abs() <= 1e-8 * a.abs().max(b.abs()).max(1e-12));
            }
            assert!(after.energy_residual_j.abs() <= after.energy_tolerance_j);
        }
    });
}

#[test]
fn g1_shared_body_coordinates_match_complete_pressure_and_original_point_history() {
    scoped(|cx,gate| {
        let shared_text=include_str!("../tensor-body.body")
            .replace("loss tensor-maxwell ","loss tensor-maxwell-shared ")
            .replace("1 6 36 67108864","1 6 36 67108864 1048576");
        let mut full=build(cx,gate,3e8,0.001);
        let mut shared=build_profile(cx,gate,3e8,0.001,Profile::parse(shared_text.as_bytes()).unwrap());
        let image=shared.performance.bank().source().bank().material_chart_tensor_relaxation().unwrap();
        assert_eq!(image.shared_coordinates().unwrap().original_coordinates(),36);
        assert!(image.coupling_rows_per_s().len()<36);
        assert_eq!(full.performance.bank().system().components()[1].modes(),
            shared.performance.bank().system().components()[1].modes());
        let p=pressure(&mut full,37);let actual=pressure(&mut shared,128);
        let scale=p.iter().map(|p|p.abs()).fold(0.0_f64,f64::max);
        assert!(scale>1e-15);
        assert!(p.iter().zip(&actual).all(|(a,b)|(a-b).abs()<=1e-8*scale));
        let source=shared.performance.bank().source().bank().material_chart_tensor_relaxation().unwrap();
        let full_image=full.performance.bank().source().bank().material_chart_tensor_relaxation().unwrap();
        let full_system=full.performance.bank().system();let system=shared.performance.bank().system();
        let q:Vec<_>=system.components()[1].states().iter()
            .map(|s|3.0_f64.sqrt()*s.displacement_m_sqrt_kg).collect();
        for (i,point) in source.material_points().iter().enumerate() {
            let complete=full_image.point_memory(i,full_system.material_memory()).unwrap();
            let reconstructed=source.point_memory(i,system.material_memory()).unwrap();
            for (a,b) in complete.iter().zip(&reconstructed) {
                assert!((a-b).abs()<1e-10*(1.0+a.abs()));
            }
            let surface:Vec<_>=reconstructed.iter().map(|z|3.0_f64.sqrt()*z).collect();
            assert!(point.image().strain(&q,&surface).is_ok());
            assert!(point.image().normal_stress_pa(&q,&surface).unwrap().iter().all(|s|s.abs()<1e-5));
        }
        let before=full_system.last_frame().unwrap();let after=system.last_frame().unwrap();
        assert!((before.network_energy_j-after.network_energy_j).abs()<1e-10);
        assert!((before.network_dissipation_j-after.network_dissipation_j).abs()<1e-10);
        assert!(after.energy_residual_j.abs()<=after.energy_tolerance_j);
        let mut replay=build_profile(cx,gate,3e8,0.001,Profile::parse(shared_text.as_bytes()).unwrap());
        assert_eq!(actual,pressure(&mut replay,1));
    });
}

#[test]
fn g3_shared_tensor_coordinates_admit_finer_original_mesh_at_unchanged_state_cap() {
    scoped(|cx,gate| {
        let text=include_str!("../tensor-body-shared.body");
        let mut ready=build_profile(cx,gate,3e8,0.001,Profile::parse(text.as_bytes()).unwrap());
        let image=ready.performance.bank().source().bank().material_chart_tensor_relaxation().unwrap();
        assert_eq!(image.material_points().len(),24);
        assert_eq!(image.shared_coordinates().unwrap().original_coordinates(),144);
        assert_eq!(image.limits().max_memory_states,36);
        assert!(image.coupling_rows_per_s().len()<=36);
        assert_eq!(image.mode_count(),2);
        let p=pressure(&mut ready,37);assert!(p.iter().any(|p|p.abs()>1e-15));
        for refused in [
            text.replace("tensor-maxwell-shared","tensor-maxwell").replace(" 1048576",""),
            text.replace("1 24 36 67108864 1048576","1 23 36 67108864 1048576"),
            text.replace("1 24 36 67108864 1048576","1 24 0 67108864 1048576"),
            text.replace("1 24 36 67108864 1048576","1 24 36 67108864 1"),
        ] {
            let result=tensor_spec().compile_with_body(cx,gate,
                Some(Profile::parse(refused.as_bytes()).unwrap()),|name|Ok(if name=="tensor-blade.fsmatpk" {
                    tensor_pack(3e8,0.001,0.02,"",100000.0)
                } else {pack(name,1.0,None)}));
            assert!(result.is_err());
        }
        let mut changed=build_profile(cx,gate,8e9,0.001,Profile::parse(text.as_bytes()).unwrap());
        assert!(p.iter().zip(pressure(&mut changed,128)).any(|(a,b)|(a-b).abs()>1e-15));
        for (upper,time,limit) in [(100.0,0.001,0.02),(100000.0,0.0,0.02),(100000.0,0.001,1e-20)] {
            assert!(tensor_spec().compile_with_body(cx,gate,Some(Profile::parse(text.as_bytes()).unwrap()),
                |name|Ok(if name=="tensor-blade.fsmatpk" {tensor_pack(3e8,time,limit,"",upper)}
                    else {pack(name,1.0,None)})).is_err());
        }
        for malformed in [text.replace(" 1048576",""),text.replace(" 1048576"," 0"),
            text.replace(" 1048576"," 67108865"),format!("{text}ignored\n")] {
            assert!(Profile::parse(malformed.as_bytes()).is_err());
        }
    });
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
        for (c33, time) in [(8e9, 0.001), (3e8, 0.005), (3e8, 1e-7)] {
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
            ("", 100000.0, 0.0, 0.02),
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
