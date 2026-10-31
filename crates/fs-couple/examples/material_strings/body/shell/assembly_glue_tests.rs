//! G0/G3: original adhesive source data reach the assembled instrument owner.
//! All geometry and coefficients here are manufactured integration fixtures.
use super::*;
use fs_couple::plate_viscosity::{
    GLUE_NORMAL_STIFFNESS_PROPERTY, GLUE_NORMAL_VISCOSITY_PROPERTY,
    GLUE_SHEAR_STIFFNESS_1_PROPERTY, GLUE_SHEAR_STIFFNESS_2_PROPERTY,
    GLUE_SHEAR_VISCOSITY_1_PROPERTY, GLUE_SHEAR_VISCOSITY_2_PROPERTY,
};
use fs_couple::shear_shell::connections::{stiffness_per_area_dims, viscosity_per_area_dims};
use fs_matdb::{NormalizedInterfacePack, SurfaceSpec, SystemContext};
use fs_plate::shell::junction::TieKind;

fn text() -> String {
    ASSEMBLED.replace(
        "slave-node 6\n",
        "slave-node 6\nglue joint.fsintpk left right 0.012\n",
    )
}
fn face_text() -> String {
    // Coincident midsurfaces are a zero-offset manufactured patch, not a corpus.
    let (left, right) = ASSEMBLED.split_once("component right\n").unwrap();
    let right = right
        .replace("node 0.5 ", "node 0 ")
        .replace("node 0.75 ", "node 0.25 ")
        .replace("node 1 ", "node 0.5 ");
    let combined = format!("{left}component right\n{right}");
    let (components, _) = combined.split_once("edge-ties 1\n").unwrap();
    let mut result =
        format!("{components}edge-ties 0\nface-ties 1\nface-tie left right midsurface full 8 8\n");
    for side in ["master", "slave"] {
        for triangle in 0..8 {
            result.push_str(&format!("{side}-triangle {triangle}\n"));
        }
    }
    result.push_str("glue joint.fsintpk left right 0.012\nbridges 1\nbridge right 4 0 0 1\n");
    result
}
fn performance() -> input::Spec {
    let mut s = spec();
    s.omega = 1.0;
    s
}
fn glue_pack(
    k: f64,
    c: f64,
    roles: [&str; 2],
    upper: f64,
    omit: Option<&str>,
    frequency_dependent: bool,
) -> NormalizedInterfacePack {
    let artifact = fs_blake3::ContentHash([17; 32]);
    let provenance = Provenance {
        source: "manufactured glue joint fixture".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic joined shell pair".into(),
            method: "manufactured Kelvin-Voigt interface".into(),
            artifact,
            caveats: "not measured instrument adhesive data".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, value, dims) in [
        (GLUE_NORMAL_STIFFNESS_PROPERTY, k, stiffness_per_area_dims()),
        (
            GLUE_SHEAR_STIFFNESS_1_PROPERTY,
            1.5 * k,
            stiffness_per_area_dims(),
        ),
        (
            GLUE_SHEAR_STIFFNESS_2_PROPERTY,
            0.4 * k,
            stiffness_per_area_dims(),
        ),
        (GLUE_NORMAL_VISCOSITY_PROPERTY, c, viscosity_per_area_dims()),
        (
            GLUE_SHEAR_VISCOSITY_1_PROPERTY,
            0.3 * c,
            viscosity_per_area_dims(),
        ),
        (
            GLUE_SHEAR_VISCOSITY_2_PROPERTY,
            1.7 * c,
            viscosity_per_area_dims(),
        ),
    ] {
        if omit == Some(name) {
            continue;
        }
        let curve = frequency_dependent && name == GLUE_NORMAL_STIFFNESS_PROPERTY;
        let value = if curve {
            PropertyValue::Curve {
                abscissa: "omega".into(),
                abscissa_dims: fs_qty::Frequency::DIMS,
                knots: vec![(0.0, value), (upper, 2.0 * value)],
                dims,
            }
        } else {
            PropertyValue::Scalar { value, dims }
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, upper),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if curve {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let claims = NormalizedPack::new(
        "synthetic-glue",
        "synthetic-interface-v1",
        artifact,
        "CC0-1.0 synthetic redistribution permitted",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: format!("source-pack:{name}"),
            phase: "declared solid model".into(),
            process: "frozen source-pack realization".into(),
            revision: 0,
        },
        texture_frame: "authored mortar frame".into(),
    };
    let pack = NormalizedInterfacePack::new(
        surface(roles[0]),
        surface(roles[1]),
        SystemContext {
            medium: "bonded dry interface".into(),
            third_body: Some("synthetic adhesive".into()),
            environment: "fixed solid temperature".into(),
            history: "cured unloaded reference".into(),
        },
        claims,
    )
    .unwrap();
    NormalizedInterfacePack::from_bytes(&pack.to_bytes()).unwrap()
}
fn source(k: f64, c: f64) -> NormalizedInterfacePack {
    glue_pack(k, c, ["body.fsmatpk"; 2], 200_000.0, None, false)
}
fn loaded(
    s: &input::Spec,
    gate: &CancelGate,
    text: &str,
    pack: NormalizedInterfacePack,
) -> Profile {
    let mut profile = Profile::parse(text.as_bytes()).unwrap();
    let mut reads = 0;
    let bytes = pack.to_bytes().len();
    profile
        .load_interfaces(s, gate, |name| {
            assert_eq!(name, "joint.fsintpk");
            reads += 1;
            Ok(pack.clone())
        })
        .unwrap();
    assert_eq!(reads, 1);
    assert_eq!(profile.interface_pack_bytes().unwrap(), bytes);
    profile
}
fn glue_ready(cx: &Cx<'_>, gate: &CancelGate, k: f64, c: f64) -> Ready {
    let s = performance();
    let profile = loaded(&s, gate, &text(), source(k, c));
    s.compile_with_body(cx, gate, Some(profile), |name| Ok(pack(name, 0.0, None)))
        .unwrap()
}

#[test]
fn g3_glue_source_changes_stiffness_loss_and_pressure_on_the_same_material_strings() {
    scoped(|cx, gate| {
        let mut dry = glue_ready(cx, gate, 2e7, 0.0);
        let mut wet = glue_ready(cx, gate, 2e7, 2000.0);
        let mut stiff = glue_ready(cx, gate, 8e7, 2000.0);
        let modes = |r: &Ready| {
            r.performance.bank().system().components()[1]
                .modes()
                .to_vec()
        };
        assert_eq!(modes(&dry), modes(&wet));
        assert_ne!(modes(&wet), modes(&stiff));
        let matrix = |r: &Ready| {
            r.performance
                .bank()
                .source()
                .bank()
                .material_shell_assembly_image()
                .unwrap()
                .damping()
                .matrix_per_s()
                .to_vec()
        };
        assert!(matrix(&dry).iter().all(|v| *v == 0.0));
        assert!(matrix(&wet).iter().any(|v| v.abs() > 1e-12));
        let tie = &wet
            .performance
            .bank()
            .source()
            .bank()
            .material_shell_assembly()
            .unwrap()
            .specimen
            .edge_ties[0];
        let TieKind::Compliant(glue) = tie.kind else {
            panic!("source glue became rigid");
        };
        assert_eq!(glue.normal_n_m3, 2e7);
        assert_eq!(glue.shear_n_m3, [3e7, 8e6]);
        assert_eq!(glue.viscosity_n_s_m3, Some([2000.0, 600.0, 3400.0]));
        assert_eq!(glue.width_m, 0.012);
        for other in [&dry, &stiff] {
            assert_eq!(
                wet.performance.bank().cards()[0].fundamental_hz(),
                other.performance.bank().cards()[0].fundamental_hz()
            );
            assert_eq!(
                wet.catalog["wire.fsmatpk"].pack.content_hash(),
                other.catalog["wire.fsmatpk"].pack.content_hash()
            );
        }
        let p = pressure(&mut wet, 7);
        assert!(p.iter().any(|p| p.abs() > 1e-14));
        assert_ne!(p, pressure(&mut dry, 31));
        assert_ne!(p, pressure(&mut stiff, 128));
        for r in [&dry, &wet, &stiff] {
            let frame = r.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        let mut replay = glue_ready(cx, gate, 2e7, 2000.0);
        let expected =
            fs_couple::pcm_wav::encode_pcm16_wav(&p, replay.spec.rate, replay.spec.full_scale)
                .unwrap()
                .0;
        let (wav, _) = replay
            .write(std::io::Cursor::new(Vec::<u8>::new()), gate)
            .unwrap();
        assert_eq!(wav.into_inner(), expected);
    });
}

#[test]
fn g0_glue_addresses_and_ordered_bulk_materials_are_not_inferred() {
    for (from, to) in [
        ("joint.fsintpk", "../joint.fsintpk"),
        ("left right 0.012", "missing right 0.012"),
        ("left right 0.012", "left right 0"),
        ("master-node 5", "master-node 4"),
    ] {
        assert!(
            Profile::parse(text().replace(from, to).as_bytes()).is_err(),
            "accepted {to}"
        );
    }
    scoped(|cx, gate| {
        let s = performance();
        let profile = loaded(&s, gate, &face_text(), source(2e7, 2000.0));
        let patch = s
            .compile_with_body(cx, gate, Some(profile), |name| Ok(pack(name, 0.0, None)))
            .unwrap();
        let body = patch
            .performance
            .bank()
            .source()
            .bank()
            .material_shell_assembly()
            .unwrap();
        assert!(body.specimen.edge_ties.is_empty());
        assert!(matches!(
            body.specimen.face_ties[0].kind,
            TieKind::Compliant(_)
        ));
        assert!(
            Profile::parse(
                face_text()
                    .replace("slave-triangle 7", "slave-triangle 8")
                    .as_bytes()
            )
            .is_err()
        );
        let distinct = text().replace("region right @body", "region right right-body.fsmatpk");
        for (roles, accepted) in [
            (["body.fsmatpk", "right-body.fsmatpk"], true),
            (["right-body.fsmatpk", "body.fsmatpk"], false),
            (["body.fsmatpk", "unrelated.fsmatpk"], false),
        ] {
            let s = performance();
            let profile = loaded(
                &s,
                gate,
                &distinct,
                glue_pack(2e7, 2000.0, roles, 200_000.0, None, false),
            );
            let result =
                s.compile_with_body(cx, gate, Some(profile), |name| Ok(pack(name, 0.0, None)));
            assert_eq!(result.is_ok(), accepted);
            if !accepted {
                assert!(
                    result
                        .err()
                        .unwrap()
                        .to_string()
                        .contains("ordered glue card surface mismatch")
                );
            }
        }
        let mut reads = 0;
        assert!(
            performance()
                .compile_with_body(
                    cx,
                    gate,
                    Some(Profile::parse(text().as_bytes()).unwrap()),
                    |name| {
                        reads += 1;
                        Ok(pack(name, 0.0, None))
                    }
                )
                .is_err()
        );
        assert_eq!(
            reads, 0,
            "an unloaded glue declaration must refuse before bulk I/O"
        );
    });
}

#[test]
fn g0_glue_requires_complete_time_domain_coefficients_and_replacement_is_atomic() {
    scoped(|cx, gate| {
        let s = performance();
        let mut profile = loaded(&s, gate, &text(), source(2e7, 2000.0));
        let bytes = profile.interface_pack_bytes().unwrap();
        for invalid in [
            glue_pack(
                2e7,
                2000.0,
                ["body.fsmatpk"; 2],
                200_000.0,
                Some(GLUE_SHEAR_VISCOSITY_2_PROPERTY),
                false,
            ),
            glue_pack(2e7, 2000.0, ["body.fsmatpk"; 2], 200_000.0, None, true),
        ] {
            assert!(
                profile
                    .load_interfaces(&s, gate, |_| Ok(invalid.clone()))
                    .is_err()
            );
            assert_eq!(profile.interface_pack_bytes().unwrap(), bytes);
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            profile
                .load_interfaces(&s, &cancelled, |_| Err("must not read".into()))
                .is_err()
        );
        let mut kept = s
            .compile_with_body(cx, gate, Some(profile), |name| Ok(pack(name, 0.0, None)))
            .unwrap();
        let mut fresh = glue_ready(cx, gate, 2e7, 2000.0);
        assert_eq!(pressure(&mut kept, 7), pressure(&mut fresh, 31));
    });
}

#[test]
fn g0_glue_checks_actual_modes_and_rejects_stale_solid_state() {
    scoped(|cx, gate| {
        let reference = glue_ready(cx, gate, 2e7, 2000.0);
        let image = reference
            .performance
            .bank()
            .source()
            .bank()
            .material_shell_assembly_image()
            .unwrap();
        let upper = 0.9 * fs_math::det::sqrt(image.eigenvalue_intervals().last().unwrap().0);
        assert!(upper > 1.0);
        let s = performance();
        let narrow = loaded(
            &s,
            gate,
            &text(),
            glue_pack(2e7, 2000.0, ["body.fsmatpk"; 2], upper, None, false),
        );
        let error = s
            .compile_with_body(cx, gate, Some(narrow), |name| Ok(pack(name, 0.0, None)))
            .err()
            .expect("an admitted query is not support for every assembled mode")
            .to_string();
        assert!(error.contains("original glue frequency support"), "{error}");
        let mut s = performance();
        let stale = loaded(&s, gate, &text(), source(2e7, 2000.0));
        s.temperature += 1.0;
        let error = s
            .compile_with_body(cx, gate, Some(stale), |name| Ok(pack(name, 0.0, None)))
            .err()
            .expect("a changed solid state requires original source resolution")
            .to_string();
        assert!(
            error.contains("reloaded at the current declared solid state"),
            "{error}"
        );
    });
}

fn held_environment_glue_pack() -> NormalizedInterfacePack {
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, MoistureBasis, QuantityKind, SemanticType, ValueForm}};
    let base = source(2e7, 2000.0);
    let original = base.claims_pack().claims();
    let mut claims = ClaimSet::new();
    for id in original.observation_ids() {
        assert_eq!(claims.register_observation(original.observation(id).unwrap().clone()).unwrap(), id);
    }
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for (_, original) in original.claims_ordered() {
        let mut claim = original.clone();
        let PropertyValue::Scalar { value, dims } = claim.value else { unreachable!() };
        let factor = if dims == stiffness_per_area_dims() { 4.0 } else { 2.0 };
        claim.value = PropertyValue::Curve { abscissa: "moisture_content".into(),
            abscissa_dims: Dims::NONE, knots: vec![(0.0, value), (0.2, factor*value)], dims };
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 280.0, 320.0)
            .with_quantity("omega", quantity(QuantityKind::Frequency(FrequencyConvention::Angular)), 0.0, 200_000.0)
            .with_quantity("moisture_content", quantity(QuantityKind::Moisture(MoistureBasis::DryMass)), 0.0, 0.2);
        claim.interpolation = InterpolationPolicy::LinearInside;
        claims.insert_claim(claim).unwrap();
    }
    let claims = NormalizedPack::new("synthetic-glue-environment", "synthetic-interface-v1",
        base.claims_pack().source_artifact(), "CC0-1.0", claims, vec![], vec![]).unwrap();
    NormalizedInterfacePack::new(base.card().surface_a().clone(), base.card().surface_b().clone(),
        base.card().context().clone(), claims).unwrap()
}

#[test]
fn g3_held_glue_environment_changes_the_assembled_body_and_rejects_stale_selection() {
    scoped(|cx, gate| {
        // Original synthetic artifact, with no inferred moisture transport.
        let source = held_environment_glue_pack();
        let state = |moisture| {
            let mut spec = performance();
            spec.material_states.insert("joint.fsintpk".into(),
                vec![("moisture_content".into(), moisture)]);
            spec
        };
        let prepare = |moisture| {
            let spec = state(moisture);
            let profile = loaded(&spec, gate, &text(), source.clone());
            assert!(profile.has_interface_source("joint.fsintpk"));
            spec.compile_with_body(cx, gate, Some(profile), |name| Ok(pack(name, 0.0, None))).unwrap()
        };
        let mut low = prepare(0.0);
        let mut high = prepare(0.2);
        let joint = |ready: &Ready| {
            let tie = &ready.performance.bank().source().bank().material_shell_assembly().unwrap()
                .specimen.edge_ties[0];
            let TieKind::Compliant(glue) = tie.kind else { panic!("original glue joint disappeared") };
            glue
        };
        let a = joint(&low);
        let b = joint(&high);
        assert_eq!(b.normal_n_m3, 4.0*a.normal_n_m3);
        assert_eq!(b.viscosity_n_s_m3.unwrap(), a.viscosity_n_s_m3.unwrap().map(|value| 2.0*value));
        assert_eq!(low.catalog["body.fsmatpk"].pack.content_hash(),
            high.catalog["body.fsmatpk"].pack.content_hash());
        let low_pressure = pressure(&mut low, 7);
        let high_pressure = pressure(&mut high, 31);
        assert_ne!(low_pressure, high_pressure);
        assert_eq!(high_pressure, pressure(&mut prepare(0.2), 1));
        let mut profile = loaded(&state(0.0), gate, &text(), source.clone());
        let original_bytes = profile.interface_pack_bytes().unwrap();
        assert!(profile.load_interfaces(&state(0.3), gate, |_| Ok(source.clone())).is_err());
        assert_eq!(profile.interface_pack_bytes().unwrap(), original_bytes);
        assert!(state(0.2).compile_with_body(cx, gate, Some(profile.clone()),
            |name| Ok(pack(name, 0.0, None))).is_err());
        state(0.0).compile_with_body(cx, gate, Some(profile),
            |name| Ok(pack(name, 0.0, None))).unwrap();
        let mut missing = Profile::parse(text().as_bytes()).unwrap();
        assert!(missing.load_interfaces(&performance(), gate, |_| Ok(source.clone())).is_err());
    });
}
