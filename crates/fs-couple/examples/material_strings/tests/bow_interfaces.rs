//! G3: independent source/interface changes reach the actual coupled renderer.
use super::*;
use fs_matdb::{SurfaceSpec, SystemContext};

fn interface(mu: f64, surface: &str, p_min: f64) -> NormalizedInterfacePack {
    let artifact = fs_blake3::ContentHash([11; 32]);
    let provenance = Provenance {
        source: "synthetic bow-interface test".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("mu_static", 0.8, Dims::NONE),
        ("mu_kinetic", mu, Dims::NONE),
        ("stribeck_speed", 0.1, fs_qty::Velocity::DIMS),
        (
            "viscous_per_speed",
            0.01,
            Dims::NONE.saturating_minus(fs_qty::Velocity::DIMS),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("p", p_min, 1e6)
                    .with("v", 0.0, 10.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let claims = NormalizedPack::new(
        "synthetic-bow-interface",
        "synthetic-v1",
        artifact,
        "CC0-1.0 synthetic redistribution permitted",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    let material = |chemistry| MaterialStateId {
        chemistry,
        phase: "declared solid model".into(),
        process: "frozen source-pack realization".into(),
        revision: 0,
    };
    let source = NormalizedInterfacePack::new(
        SurfaceSpec {
            material: material(format!("source-pack:{surface}")),
            texture_frame: "synthetic string surface".into(),
        },
        SurfaceSpec {
            material: material("synthetic rosined bow".into()),
            texture_frame: "synthetic hair direction".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: Some("synthetic rosin".into()),
            environment: "declared gas".into(),
            history: "synthetic fixed history".into(),
        },
        claims,
    )
    .unwrap();
    NormalizedInterfacePack::from_bytes(&source.to_bytes()).unwrap()
}
fn declaration(track: &str, path: &str) -> String {
    format!("bow-interface {track} {path} 293.15 0.0001 2 0 coulomb mu_static mu_kinetic")
}
fn source_spec() -> input::Spec {
    let text = String::from_utf8(spec().bytes)
        .unwrap()
        .replace(
            "friction stribeck 0.8 0.4 0.1 0",
            &format!(
                "friction interfaces 1\n{}",
                declaration("bow", "bow.fsintpk")
            ),
        )
        .replace(
            "initial\tbow\t0e0\t0e0\t2e-1",
            "initial\tbow\t3e-1\t1e-4\t2e-1",
        );
    input::Spec::parse(text.as_bytes()).unwrap()
}
fn sourced(cx: &Cx<'_>, gate: &CancelGate, mut spec: input::Spec, mu: f64, surface: &str) -> Ready {
    spec.friction
        .load(gate, |_| Ok(interface(mu, surface, 0.1)))
        .unwrap();
    compile(cx, gate, spec, 7800.0)
}

#[test]
fn canonical_material_card_binds_the_same_named_surface_to_the_bow() {
    scoped(|cx, gate| {
        let wire = super::material_cards::canonical("wire.fsmatpk", 7800.0);
        let legacy = interface(0.4, "wire.fsmatpk", 0.1);
        let mut surface = legacy.card().surface_a().clone();
        surface.material = wire.card().id().clone();
        let source = NormalizedInterfacePack::new(surface,
            legacy.card().surface_b().clone(), legacy.card().context().clone(),
            legacy.claims_pack().clone()).unwrap();
        let mut input = source_spec();
        input.friction.load(gate, |_| Ok(source.clone())).unwrap();
        let mut ready = input.compile(cx, gate, |name| Ok(if name == "wire.fsmatpk" {
            wire.clone()
        } else { super::material_cards::canonical(name, 7800.0) })).unwrap();
        assert_eq!(ready.catalog["wire.fsmatpk"].card.id(), wire.card().id());
        assert!(ready.performance.bank().bow_interface(0).is_some());
        assert!(pressure(&mut ready, 37).iter().any(|p| p.abs() > 1e-10));
        let mut mismatched = source_spec();
        mismatched.friction.load(gate, |_| Ok(legacy.clone())).unwrap();
        let error = mismatched.compile(cx, gate, |name|
            Ok(super::material_cards::canonical(name, 7800.0))).err().unwrap();
        assert!(error.to_string().contains("does not match"), "{error}");
    });
}

fn ribbon_spec(spec: input::Spec, record: &str) -> Result<input::Spec> {
    let text = String::from_utf8(spec.bytes)?.replace(
        "\nforces ", &format!("\nbow_ribbon {record}\nforces "));
    input::Spec::parse(text.as_bytes())
}

#[test]
fn ribbon_file_keeps_original_track_source_and_pressure_while_resolving_spatial_slip() {
    scoped(|cx, gate| {
        let mut point = sourced(cx, gate, source_spec(), 0.4, "wire.fsmatpk");
        let mut single = sourced(cx, gate, ribbon_spec(source_spec(), "0.04 1").unwrap(),
            0.4, "wire.fsmatpk");
        let mut wide = sourced(cx, gate, ribbon_spec(source_spec(), "0.04 5").unwrap(),
            0.4, "wire.fsmatpk");
        assert_eq!(wide.performance.schedule().to_canonical_bytes(),
            point.performance.schedule().to_canonical_bytes());
        assert_eq!(wide.performance.schedule().tracks().len(), 1);
        assert_eq!(wide.performance.bank().source().joint().max_ports, 5);
        let ribbon = wide.performance.bank().source().ribbon().unwrap();
        assert_eq!(ribbon.width_m(), 0.04);
        assert_eq!(ribbon.stations().len(), 5);
        let original = point.performance.bank().bow_interface(0).unwrap();
        for port in 0..5 {
            let binding = wide.performance.bank().bow_interface(port).unwrap();
            assert_eq!(binding.card(), original.card());
            assert_eq!(binding.temperature_kelvin(), original.temperature_kelvin());
            assert_eq!(binding.nominal_area_m2(), original.nominal_area_m2() / 5.0);
        }
        assert_ne!(wide.performance.bank().system().port(0).unwrap().left.shapes,
            wide.performance.bank().system().port(1).unwrap().left.shapes);
        assert_eq!(wide.performance.events().len(), point.performance.events().len() * 5);
        let p = pressure(&mut point, 37);
        assert_eq!(p, pressure(&mut single, 128));
        let q = pressure(&mut wide, 19);
        assert!(p.iter().zip(q).any(|(a, b)| (a - b).abs() > 1e-15));
        let frame = wide.performance.bank().system().last_frame().unwrap();
        let total: f64 = frame.ports.iter().map(|p| p.normal_force_n).sum();
        assert!((total - 0.0001).abs() < 1e-18);
        for (port, contact) in frame.ports.iter().enumerate() {
            let area = wide.performance.bank().bow_interface(port).unwrap().nominal_area_m2();
            assert!((contact.normal_force_n / area - 1.0).abs() < 1e-12);
        }
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn ribbon_source_rejects_invalid_width_resolution_or_ambiguous_normal_loads() {
    for record in ["0 3", "NaN 3", "0.01 0", "0.01 2", "0.01 65", "0.2 3"] {
        assert!(ribbon_spec(spec(), record).is_err(), "accepted {record}");
    }
    let source = spec();
    let mut tracks = source.schedule.tracks().to_vec();
    let mut second = tracks[0].clone(); second.id = "second".into(); tracks.push(second);
    let schedule = GestureSchedule::try_new(source.schedule.control_rate_hz, tracks).unwrap();
    let text = String::from_utf8(source.bytes).unwrap();
    let header = text.split_once("\nschedule\n").unwrap().0.replace("\nforces ",
        "\nbow_ribbon 0.01 3\nforces ");
    let mut bytes = format!("{header}\nschedule\n").into_bytes();
    bytes.extend(schedule.to_canonical_bytes());
    assert!(input::Spec::parse(&bytes).is_err());
    assert!(input::Spec::parse(include_bytes!("../source-bow.material-strings")).is_ok());
    scoped(|cx, gate| {
        let spec = input::Spec::parse(include_bytes!("../contact-bow.material-strings")).unwrap();
        let spec = ribbon_spec(spec, "0.01 3").unwrap();
        let stops = crate::stops::Profile::parse(include_bytes!("../contact-bow.stops")).unwrap();
        let mut reads = 0;
        assert!(spec.compile_with_profiles(cx, gate, None, Some(stops), |name| {
            reads += 1; Ok(pack(name, 7800.0, 200_000.0, None))
        }).is_err());
        assert_eq!(reads, 0);
    });
}

#[test]
fn interface_pack_changes_excitation_twist_and_pressure_with_identical_bulk_materials() {
    scoped(|cx, gate| {
        let mut a = sourced(cx, gate, source_spec(), 0.4, "wire.fsmatpk");
        let mut b = sourced(cx, gate, source_spec(), 0.2, "wire.fsmatpk");
        assert_eq!(a.spec.bytes, b.spec.bytes);
        for (x, y) in a
            .performance
            .bank()
            .system()
            .components()
            .iter()
            .zip(b.performance.bank().system().components())
        {
            assert_eq!(x.modes(), y.modes());
        }
        for name in ["wire.fsmatpk", "body.fsmatpk"] {
            assert_eq!(
                a.catalog[name].pack.content_hash(),
                b.catalog[name].pack.content_hash()
            );
        }
        let p = pressure(&mut a, 37);
        let q = pressure(&mut b, 128);
        assert!(p.iter().zip(q).any(|(p, q)| (p - q).abs() > 1e-15));
        assert_ne!(
            a.performance.bank().system().components()[0].states(),
            b.performance.bank().system().components()[0].states()
        );
        for r in [&a, &b] {
            let binding = r.performance.bank().bow_interface(0).unwrap();
            assert_eq!(
                binding.card().interface().history_id(),
                "synthetic fixed history"
            );
            assert_eq!(binding.nominal_area_m2(), 0.0001);
            let frame = r.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        let mut replay = sourced(cx, gate, source_spec(), 0.4, "wire.fsmatpk");
        let expected = fs_couple::pcm_wav::encode_pcm16_wav(&p, 48000, 0.01)
            .unwrap()
            .0;
        assert_eq!(
            replay
                .write(Cursor::new(Vec::new()), gate)
                .unwrap()
                .0
                .into_inner(),
            expected
        );
    });
}

#[test]
fn distinct_contact_sources_are_bound_by_track_name_not_declaration_order() {
    scoped(|cx, gate| {
        let source = source_spec();
        let schedule = GestureSchedule::try_new(
            700,
            vec![
                GestureTrack {
                    id: "first".into(),
                    target: GestureTarget::BowStroke { string: 0 },
                    initial: GestureValue::Bow {
                        velocity_m_per_s: 0.3,
                        normal_force_n: 0.0001,
                        station: 0.2,
                    },
                    events: vec![],
                },
                GestureTrack {
                    id: "second".into(),
                    target: GestureTarget::BowStroke { string: 0 },
                    initial: GestureValue::Bow {
                        velocity_m_per_s: -0.2,
                        normal_force_n: 0.0001,
                        station: 0.31,
                    },
                    events: vec![],
                },
            ],
        )
        .unwrap();
        let text = String::from_utf8(source.bytes).unwrap();
        let header = text.split_once("\nschedule\n").unwrap().0.replace(
            &format!(
                "friction interfaces 1\n{}",
                declaration("bow", "bow.fsintpk")
            ),
            &format!(
                "friction interfaces 2\n{}\n{}",
                declaration("second", "second.fsintpk"),
                declaration("first", "first.fsintpk")
            ),
        );
        let mut bytes = format!("{header}\nschedule\n").into_bytes();
        bytes.extend(schedule.to_canonical_bytes());
        let mut spec = input::Spec::parse(&bytes).unwrap();
        spec.friction
            .load(gate, |path| {
                Ok(interface(
                    if path == "first.fsintpk" { 0.4 } else { 0.2 },
                    "wire.fsmatpk",
                    0.1,
                ))
            })
            .unwrap();
        let mut ready = compile(cx, gate, spec, 7800.0);
        for (port, mu) in [(0, 0.4), (1, 0.2)] {
            assert_eq!(
                ready
                    .performance
                    .bank()
                    .bow_interface(port)
                    .unwrap()
                    .card()
                    .friction_law(),
                &FrictionLaw::Coulomb {
                    static_mu: 0.8,
                    kinetic_mu: mu
                }
            );
        }
        assert!(pressure(&mut ready, 19).iter().any(|p| p.abs() > 1e-15));
    });
}

#[test]
fn bonded_contact_binds_the_exposed_sleeve_and_refuses_an_old_core_interface() {
    scoped(|cx, gate| {
        let text=String::from_utf8(source_spec().bytes).unwrap()
            .replace("frankensim-material-strings-v1","frankensim-material-strings-v2")
            .replace("string wire.fsmatpk 0.5 0.0003 20 4 3 0.8 3000 0.1",
                "string wire.fsmatpk 0.5 0.0003 tension 20 4 3 0.8 3000 0.1\nsection bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003");
        let mut ready = sourced(
            cx,
            gate,
            input::Spec::parse(text.as_bytes()).unwrap(),
            0.4,
            "sleeve.fsmatpk",
        );
        assert!(pressure(&mut ready, 37).iter().any(|p| p.abs() > 1e-15));
        let mut wrong = input::Spec::parse(text.as_bytes()).unwrap();
        wrong
            .friction
            .load(gate, |_| Ok(interface(0.4, "wire.fsmatpk", 0.1)))
            .unwrap();
        assert!(
            wrong
                .compile(cx, gate, |name| Ok(pack(name, 7800.0, 200_000.0, None)))
                .is_err()
        );
    });
}

#[test]
fn absent_sources_and_runtime_pressure_outside_the_source_domain_refuse() {
    scoped(|cx, gate| {
        let mut reads = 0;
        assert!(
            source_spec()
                .compile(cx, gate, |name| {
                    reads += 1;
                    Ok(pack(name, 7800.0, 200_000.0, None))
                })
                .is_err()
        );
        assert_eq!(
            reads, 0,
            "missing interface must refuse before bulk construction"
        );
        let mut spec = source_spec();
        // Query at 2 Pa is supported, but actual N/area = 1 Pa is not.
        spec.friction
            .load(gate, |_| Ok(interface(0.4, "wire.fsmatpk", 2.0)))
            .unwrap();
        let mut ready = compile(cx, gate, spec, 7800.0);
        assert!(ready.performance.block(&mut [0.0; 1]).is_err());
        assert_eq!(ready.performance.samples_rendered(), 0);
        let text = String::from_utf8(source_spec().bytes).unwrap();
        for (from, to) in [
            ("bow-interface bow ", "bow-interface missing "),
            ("bow.fsintpk", "../bow.fsintpk"),
            ("mu_static mu_kinetic", "mu_static mu_static"),
            ("friction interfaces 1", "friction interfaces 17"),
        ] {
            assert!(input::Spec::parse(text.replace(from, to).as_bytes()).is_err());
        }
        let mut stribeck = input::Spec::parse(
            text.replace(
                "coulomb mu_static mu_kinetic",
                "stribeck mu_static mu_kinetic stribeck_speed viscous_per_speed",
            )
            .as_bytes(),
        )
        .unwrap();
        stribeck
            .friction
            .load(gate, |_| Ok(interface(0.4, "wire.fsmatpk", 0.1)))
            .unwrap();
        assert!(matches!(
            stribeck.bow_law(0).unwrap(),
            FrictionLaw::Stribeck { .. }
        ));
    });
}

mod hyperbolic {
    use super::*;

    fn source_spec_hyperbolic() -> input::Spec {
        // The existing fixture is synthetic. Its speed property is explicitly
        // selected for this manufactured law; no measured Stribeck fit is reused.
        let text = String::from_utf8(source_spec().bytes).unwrap().replace(
            "coulomb mu_static mu_kinetic",
            "hyperbolic mu_static mu_kinetic stribeck_speed",
        );
        input::Spec::parse(text.as_bytes()).unwrap()
    }
    fn numeric_text() -> String {
        let source = source_spec_hyperbolic();
        String::from_utf8(source.bytes).unwrap().replace(
            "friction interfaces 1\nbow-interface bow bow.fsintpk 293.15 0.0001 2 0 hyperbolic mu_static mu_kinetic stribeck_speed",
            "friction hyperbolic 0.8 0.4 0.1",
        )
    }

    #[test]
    fn hyperbolic_source_and_numeric_inputs_share_pressure_and_callback_history() {
        scoped(|cx, gate| {
            let mut source = sourced(cx,gate,source_spec_hyperbolic(),0.4,"wire.fsmatpk");
            let mut numeric = compile(cx,gate,input::Spec::parse(numeric_text().as_bytes()).unwrap(),7800.0);
            assert_eq!(source.spec.bow_law(0).unwrap(), &FrictionLaw::Hyperbolic {
                static_mu:0.8, kinetic_mu:0.4, characteristic_speed:0.1,
            });
            let binding = source.performance.bank().bow_interface(0).unwrap();
            assert_eq!(binding.card().interface().history_id(), "synthetic fixed history");
            assert_eq!(binding.nominal_area_m2(), 0.0001);
            assert!(numeric.performance.bank().bow_interface(0).is_none());
            let expected = pressure(&mut source,37);
            assert!(expected.iter().all(|p|p.is_finite()));
            assert!(expected.iter().any(|p|p.abs()>1e-15));
            assert_eq!(expected,pressure(&mut numeric,128));
            let mut replay = sourced(cx,gate,source_spec_hyperbolic(),0.4,"wire.fsmatpk");
            assert_eq!(expected,pressure(&mut replay,1));
        });
    }

    #[test]
    fn hyperbolic_interface_change_reaches_motion_without_replacing_bulk_material() {
        scoped(|cx, gate| {
            let mut a = sourced(cx,gate,source_spec_hyperbolic(),0.4,"wire.fsmatpk");
            let mut b = sourced(cx,gate,source_spec_hyperbolic(),0.2,"wire.fsmatpk");
            assert_eq!(a.spec.bytes,b.spec.bytes);
            for name in ["wire.fsmatpk","body.fsmatpk"] {
                assert_eq!(a.catalog[name].pack.content_hash(),b.catalog[name].pack.content_hash());
            }
            for (x,y) in a.performance.bank().system().components().iter()
                .zip(b.performance.bank().system().components()) {
                assert_eq!(x.modes(),y.modes());
            }
            let x=pressure(&mut a,19); let y=pressure(&mut b,128);
            assert!(x.iter().zip(&y).any(|(x,y)|(x-y).abs()>1e-15));
            assert_ne!(a.performance.bank().system().components()[0].states(),
                b.performance.bank().system().components()[0].states());
            for ready in [&a,&b] {
                let frame=ready.performance.bank().system().last_frame().unwrap();
                assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            }
        });
    }

    #[test]
    fn hyperbolic_selection_keeps_source_and_runtime_refusals() {
        let numeric=numeric_text();
        for record in ["hyperbolic 0.8 0.4 0", "hyperbolic 0.8 0.9 0.1",
            "hyperbolic -0.8 0.4 0.1", "hyperbolic 0.8 0.4 NaN",
            "hyperbolic 0.8 0.4 0.1 0"] {
            assert!(input::Spec::parse(numeric.replace("hyperbolic 0.8 0.4 0.1",record).as_bytes()).is_err());
        }
        scoped(|cx,gate| {
            let text=String::from_utf8(source_spec_hyperbolic().bytes).unwrap();
            let mut missing=input::Spec::parse(text.replace("stribeck_speed","missing_speed").as_bytes()).unwrap();
            assert!(missing.friction.load(gate,|_|Ok(interface(0.4,"wire.fsmatpk",0.1))).is_err());
            let mut invalid=source_spec_hyperbolic();
            invalid.friction.load(gate,|_|Ok(interface(0.4,"wire.fsmatpk",2.0))).unwrap();
            let mut ready=compile(cx,gate,invalid,7800.0);
            assert!(ready.performance.block(&mut [0.0;1]).is_err());
            assert_eq!(ready.performance.samples_rendered(),0);
            let mut reads=0;
            assert!(source_spec_hyperbolic().compile(cx,gate,|name| {
                reads+=1; Ok(pack(name,7800.0,200_000.0,None))
            }).is_err());
            assert_eq!(reads,0);
        });
    }
}


mod held_environment {
    use super::*;
    use fs_qty::{QuantitySpec, semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm}};

    // Manufactured source laws exercise transport and coupling; they are not
    // a measured rosin/humidity model. The same artifact serves every state.
    fn source(response: Option<&str>) -> NormalizedInterfacePack {
        let base = interface(0.4, "wire.fsmatpk", 0.1);
        let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let mut claims = ClaimSet::new();
        for id in base.claims_pack().claims().observation_ids() {
            claims.register_observation(base.claims_pack().claims().observation(id).unwrap().clone()).unwrap();
        }
        for (_, original) in base.claims_pack().claims().claims_ordered() {
            let mut claim = original.clone();
            let PropertyValue::Scalar { value, dims } = claim.value else { unreachable!() };
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 280.0, 320.0)
                .with_quantity("p", quantity(QuantityKind::Pressure), 0.1, 1e6)
                .with_quantity("v", QuantitySpec::dimensional(fs_qty::Velocity::DIMS), 0.0, 10.0)
                .with_quantity("moisture_content", quantity(QuantityKind::Moisture(MoistureBasis::DryMass)), 0.0, 0.2);
            claim.value = if let Some(axis) = response {
                let (lo, hi) = claim.validity.bound(axis).unwrap();
                PropertyValue::Surface2D {
                    axes: [axis.into(), "moisture_content".into()],
                    axis_dims: [claim.validity.axis_quantities()[axis].dims(), Dims::NONE],
                    samples: vec![([lo, 0.0], value), ([hi, 0.0], 0.9*value),
                        ([lo, 0.2], 0.5*value), ([hi, 0.2], 0.4*value)],
                    triangles: vec![[0, 1, 2], [1, 3, 2]], dims,
                }
            } else {
                PropertyValue::Curve { abscissa: "moisture_content".into(), abscissa_dims: Dims::NONE,
                    knots: vec![(0.0, value), (0.2, 0.5*value)], dims }
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
            claims.insert_claim(claim).unwrap();
        }
        let claims = NormalizedPack::new("synthetic-bow-environment", "synthetic-v1",
            base.claims_pack().source_artifact(), "CC0-1.0", claims, vec![], vec![]).unwrap();
        NormalizedInterfacePack::new(base.card().surface_a().clone(), base.card().surface_b().clone(),
            base.card().context().clone(), claims).unwrap()
    }

    fn input(moisture: f64) -> input::Spec {
        let text = String::from_utf8(source_spec().bytes).unwrap().replace("state 293.15 1000",
            &format!("state 293.15 1000\nmaterial-states 1\nmaterial-state bow.fsintpk 1 moisture_content {moisture}"));
        input::Spec::parse(text.as_bytes()).unwrap()
    }
    fn load(s: &mut input::Spec, gate: &CancelGate, source: &NormalizedInterfacePack) -> Result<()> {
        s.friction.load_with_environment(&s.material_states, gate, |_| Ok(source.clone()))
    }

    #[test]
    fn one_bow_source_changes_actual_pressure_with_held_moisture_and_preserves_callback_history() {
        scoped(|cx, gate| {
            let source = source(None);
            let make = |moisture| {
                let mut s = input(moisture);
                load(&mut s, gate, &source).unwrap();
                compile(cx, gate, s, 7800.0)
            };
            let mut dry = make(0.0);
            let mut moist = make(0.2);
            for name in ["wire.fsmatpk", "body.fsmatpk"] {
                assert_eq!(dry.catalog[name].pack.content_hash(), moist.catalog[name].pack.content_hash());
            }
            for (a, b) in dry.performance.bank().system().components().iter()
                .zip(moist.performance.bank().system().components()) {
                assert_eq!(a.modes(), b.modes());
            }
            assert_eq!(moist.spec.bow_law(0).unwrap(),
                &FrictionLaw::Coulomb { static_mu: 0.4, kinetic_mu: 0.2 });
            let expected = pressure(&mut dry, 37);
            let changed = pressure(&mut moist, 128);
            assert!(expected.iter().zip(&changed).any(|(a,b)| (a-b).abs() > 1e-15));
            assert_eq!(changed, pressure(&mut make(0.2), 1));
            for ready in [&dry, &moist] {
                assert!(ready.performance.bank().bow_interface(0).unwrap().temperature_source().is_none());
                let frame = ready.performance.bank().system().last_frame().unwrap();
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
            for model in ["stribeck mu_static mu_kinetic stribeck_speed viscous_per_speed",
                "hyperbolic mu_static mu_kinetic stribeck_speed"] {
                let text = String::from_utf8(input(0.2).bytes).unwrap()
                    .replace("coulomb mu_static mu_kinetic", model);
                let mut s = input::Spec::parse(text.as_bytes()).unwrap();
                load(&mut s, gate, &source).unwrap();
                s.friction.validate_thermal_mode(false).unwrap();
            }
        });
    }

    #[test]
    fn thermal_bow_refresh_retains_source_context_and_refuses_dynamic_coefficients_without_their_owner() {
        scoped(|_, gate| {
            let mut s = input(0.1);
            load(&mut s, gate, &source(Some("T"))).unwrap();
            assert!(s.friction.validate_thermal_mode(false).is_err());
            s.friction.validate_thermal_mode(true).unwrap();
            let binding = s.friction.binding("bow").unwrap().unwrap();
            let owner = binding.temperature_source().unwrap();
            let point = owner.initial_source().query().with("T", 310.0).unwrap()
                .with("p", 3.0).unwrap().with("v", 0.2).unwrap();
            let (resolved, law) = owner.coefficients_at_query(&point).unwrap();
            assert_eq!(resolved.query().axes()["moisture_content"], 0.1);
            assert_eq!(resolved.card_identity(), owner.initial_source().card_identity());
            assert_ne!(law, *binding.card().friction_law());
            for property in resolved.properties() {
                assert_eq!(property.answer().receipt.selected,
                    owner.initial_source().property(property.requirement().name()).unwrap().answer().receipt.selected);
            }
            assert!(owner.coefficients_at_query(&point.with("moisture_content", 0.2).unwrap()).is_err());
            for dynamic in ["p", "v"] {
                assert!(load(&mut input(0.1), gate, &source(Some(dynamic))).is_err());
            }
        });
    }

    #[test]
    fn changed_or_missing_held_bow_state_refuses_before_bulk_loading() {
        scoped(|cx, gate| {
            let source = source(None);
            assert!(source_spec().friction.load(gate, |_| Ok(source.clone())).is_err());
            let mut s = input(0.1);
            load(&mut s, gate, &source).unwrap();
            s.material_states.get_mut("bow.fsintpk").unwrap()[0].1 = 0.2;
            let mut reads = 0;
            assert!(s.compile(cx, gate, |name| {
                reads += 1; Ok(pack(name, 7800.0, 200_000.0, None))
            }).is_err());
            assert_eq!(reads, 0);
            assert!(load(&mut input(0.3), gate, &source).is_err());
            for axis in ["p", "v"] {
                let mut s = input(0.1);
                s.material_states.get_mut("bow.fsintpk").unwrap().push((axis.into(), 0.2));
                assert!(load(&mut s, gate, &source).is_err());
            }
        });
    }
}
