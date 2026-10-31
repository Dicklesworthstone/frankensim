//! Source-bound holding footprints through the existing file-to-pressure path.
use super::*;
use fs_couple::shear_shell::connections::{
    stiffness_per_area_dims, viscosity_per_area_dims, ShellInterfaceSide,
};
use fs_matdb::{NormalizedInterfacePack, SurfaceSpec, SystemContext};

fn held_text() -> &'static str { include_str!("../../held-arched.body") }
fn held_profile() -> Profile { Profile::parse(held_text().as_bytes()).unwrap() }

fn holding_pack(k: f64, c: f64, body: &str, upper: f64) -> NormalizedInterfacePack {
    let artifact = fs_blake3::ContentHash([9; 32]);
    let provenance = Provenance { source: "manufactured compliant holding fixture".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact) };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic shell surface against immobile support".into(),
        method: "manufactured bilateral linear interface".into(), artifact,
        caveats: "not measured holding or glue data".into(), provenance: provenance.clone(),
    }).unwrap();
    for (name, value, dims) in [
        ("normal_stiffness_per_area", k, stiffness_per_area_dims()),
        ("normal_viscosity_per_area", c, viscosity_per_area_dims()),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, upper),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation.clone()], provenance: provenance.clone(),
        }).unwrap();
    }
    let claims = NormalizedPack::new("synthetic-holding", "synthetic-interface-v1",
        artifact, "CC0-1.0 synthetic redistribution permitted", claims, vec![], vec![]).unwrap();
    let pack = NormalizedInterfacePack::new(
        SurfaceSpec {
            material: MaterialStateId { chemistry: format!("source-pack:{body}"),
                phase: "declared solid model".into(), process: "frozen source-pack realization".into(),
                revision: 0 },
            texture_frame: "authored global shell footprint".into(),
        },
        SurfaceSpec {
            material: MaterialStateId { chemistry: "synthetic support".into(),
                phase: "solid".into(), process: "immobile fixture".into(), revision: 0 },
            texture_frame: "authored fixed counterpart".into(),
        },
        SystemContext { medium: "manufactured bilateral layer".into(), third_body: None,
            environment: "fixed source-state fixture".into(), history: "stress-free reference".into() },
        claims,
    ).unwrap();
    NormalizedInterfacePack::from_bytes(&pack.to_bytes()).unwrap()
}

fn held_ready(cx: &Cx<'_>, gate: &CancelGate, k: f64, c: f64) -> Ready {
    let spec = shell_spec();
    let mut profile = held_profile();
    let mut reads = 0;
    profile.load_interfaces(&spec, gate, |name| {
        assert_eq!(name, "holding.fsintpk");
        reads += 1;
        Ok(holding_pack(k, c, "body.fsmatpk", 200_000.0))
    }).unwrap();
    assert_eq!(reads, 1);
    assert_eq!(profile.additional_body_connections(), 1);
    spec.compile_with_body(cx, gate, Some(profile), |name| Ok(pack(name, 1.0, None))).unwrap()
}

#[test]
fn authored_interface_source_changes_the_complete_body_response_and_streamed_pressure() {
    scoped(|cx, gate| {
        let mut soft = held_ready(cx, gate, 1e8, 2e4);
        let mut stiff = held_ready(cx, gate, 4e8, 2e4);
        let mut viscous = held_ready(cx, gate, 1e8, 8e4);
        let source = soft.performance.bank().source().bank();
        let shell = source.material_shell().unwrap();
        let image = source.material_shell_image().unwrap();
        let link = &shell.specimen.interfaces()[0];
        assert_eq!(link.input().name, "holding");
        assert_eq!(link.law().stiffness_n_m, 1e4);
        assert_eq!(link.law().damping_n_s_m, 2.0);
        assert_eq!(link.omega_band_rad_s(), (0.0, 200_000.0));
        assert_eq!(link.input().interface.surface_materials()[0],
            shell.specimen.regions()[0].input().region.material.material().clone());
        let ShellInterfaceSide::Fixed { material } = &link.input().right else {
            panic!("fixed ordered surface role disappeared");
        };
        assert_eq!(material, &link.input().interface.surface_materials()[1]);
        assert_eq!(material.chemistry, "synthetic support");
        let row = &image.interface_columns()[0];
        for (i, mode) in image.physical_modes().iter().enumerate() {
            let coordinate = 0.5 * (mode[6 * 11 + 2] - 0.01 * mode[6 * 11 + 4])
                + 0.5 * (mode[6 * 12 + 2] - 0.01 * mode[6 * 12 + 4]);
            assert!((row[i] - coordinate).abs() < 1e-12);
        }
        let channel = source.body_connections().last().unwrap();
        assert_eq!(&channel.left.shapes, row);
        assert_eq!(channel.stiffness_n_m, 0.0, "spring belongs in the shared shell eigenproblem");
        assert_eq!(channel.damping_n_s_m, 2.0);
        let modes = |r: &Ready| r.performance.bank().system().components()[1].modes().to_vec();
        assert_ne!(modes(&soft), modes(&stiff));
        assert_eq!(modes(&soft), modes(&viscous));
        for other in [&stiff, &viscous] {
            assert_eq!(soft.performance.bank().cards()[0].fundamental_hz(),
                other.performance.bank().cards()[0].fundamental_hz());
            assert_eq!(source.gas(), other.performance.bank().source().bank().gas());
            for name in ["wire.fsmatpk", "sleeve.fsmatpk", "body.fsmatpk"] {
                assert_eq!(soft.catalog[name].pack.content_hash(), other.catalog[name].pack.content_hash());
            }
        }
        let p = pressure(&mut soft, 7);
        assert!(p.iter().any(|v| v.abs() > 1e-14));
        assert!(p.iter().zip(pressure(&mut stiff, 31)).any(|(a, b)| *a != b));
        assert!(p.iter().zip(pressure(&mut viscous, 96)).any(|(a, b)| *a != b));
        for r in [&soft, &stiff, &viscous] {
            let frame = r.performance.bank().system().last_frame().unwrap();
            assert_eq!(frame.external_work_j, 0.0);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        let mut replay = held_ready(cx, gate, 1e8, 2e4);
        let expected = fs_couple::pcm_wav::encode_pcm16_wav(&p, replay.spec.rate, replay.spec.full_scale).unwrap().0;
        let (wav, report) = replay.write(std::io::Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), expected);
        assert_eq!(report.samples, 96);
        let profile = replay.body_profile.as_ref().unwrap();
        assert_eq!(profile.bytes, held_text().as_bytes());
        assert_eq!(profile.interface_pack_bytes().unwrap(),
            holding_pack(1e8, 2e4, "body.fsmatpk", 200_000.0).to_bytes().len());
    });
}

#[test]
fn pure_dashpot_sources_admit_while_wrong_ordered_materials_and_bands_refuse() {
    scoped(|cx, gate| {
        let mut dashpot = held_ready(cx, gate, 0.0, 2e4);
        let source = dashpot.performance.bank().source().bank();
        let link = &source.material_shell().unwrap().specimen.interfaces()[0];
        assert_eq!(link.law().stiffness_n_m, 0.0);
        assert_eq!(source.body_connections().last().unwrap().damping_n_s_m, 2.0);
        assert!(pressure(&mut dashpot, 13).iter().any(|p| p.abs() > 1e-14));

        for (body, upper) in [("other-original-body-pack", 200_000.0), ("body.fsmatpk", 1.0)] {
            let spec = shell_spec();
            let mut profile = held_profile();
            profile.load_interfaces(&spec, gate, |_| Ok(holding_pack(1e8, 2e4, body, upper))).unwrap();
            assert!(spec.compile_with_body(cx, gate, Some(profile),
                |name| Ok(pack(name, 1.0, None))).is_err());
        }
        // No interface cards are guessed from the bulk material or holding label.
        let mut reads = 0;
        assert!(shell_spec().compile_with_body(cx, gate, Some(held_profile()),
            |name| { reads += 1; Ok(pack(name, 1.0, None)) }).is_err());
        assert_eq!(reads, 0);
    });
}

#[test]
fn interface_geometry_and_complete_allowances_precede_source_io_and_paths_cache_once() {
    let text = held_text();
    for (from, to) in [
        ("interfaces 1 2 1048576", "interfaces 1 1 1048576"),
        ("interfaces 1 2 1048576", "interfaces 65 2 1048576"),
        ("holding.fsintpk", "../holding.fsintpk"),
        ("side region body 2", "side region absent 2"),
        ("point 11 0.5 0.01 0 0", "point 11 0.4 0.01 0 0"),
        ("point 12 0.5 0.01 0 0", "point 11 0.5 0.01 0 0"),
    ] {
        assert!(Profile::parse(text.replace(from, to).as_bytes()).is_err(), "accepted {to}");
    }
    // A source-authored compliant support may replace all homogeneous supports.
    // Rank and rigid-body admission remain the shared eigensolver's decision.
    let supports = text.split_once("supports 16\n").unwrap().1.split_once("interfaces ").unwrap().0;
    assert!(Profile::parse(text.replace(&format!("supports 16\n{supports}"), "supports 0\n").as_bytes()).is_ok());
    scoped(|_, gate| {
        let mut over = shell_spec();
        // More modal coordinates than the whole instrument admits, so the reservation refuses before any pack I/O.
        over.strings[0].bending_modes = input::MAX_TOTAL_MODES;
        let mut reads = 0;
        assert!(held_profile().load_interfaces(&over, gate, |_| {
            reads += 1; Ok(holding_pack(1e8, 2e4, "body.fsmatpk", 200_000.0))
        }).is_err());
        assert_eq!(reads, 0);
        let spec = shell_spec();
        let normal = spec.network_dimensions(true).unwrap();
        let connected = spec.network_dimensions_with_body_connections(true, 1).unwrap();
        assert_eq!((connected.0, connected.1), (normal.0, normal.1 + 1));
        assert_eq!(connected.2, connected.0 * (connected.1 + 1).pow(2));
        assert!(spec.network_dimensions_with_body_connections(true, MAX_MODAL_CONNECTIONS).is_err());

        let block = text.split_once("interfaces 1 2 1048576\n").unwrap().1.split_once("bridges 1\n").unwrap().0;
        let repeated = text.replace("interfaces 1 2 1048576", "interfaces 2 4 1048576")
            .replace("bridges 1\n", &format!("{}bridges 1\n",
                block.replace("interface holding ", "interface second ")));
        let mut profile = Profile::parse(repeated.as_bytes()).unwrap();
        let mut reads = 0;
        profile.load_interfaces(&spec, gate, |_| {
            reads += 1; Ok(holding_pack(1e8, 2e4, "body.fsmatpk", 200_000.0))
        }).unwrap();
        assert_eq!(reads, 1);
        assert_eq!(profile.additional_body_connections(), 2);
        let bytes = profile.interface_pack_bytes().unwrap();
        assert_eq!(bytes, holding_pack(1e8, 2e4, "body.fsmatpk", 200_000.0).to_bytes().len());
        assert!(profile.load_interfaces(&spec, gate, |_| Err("source unavailable".into())).is_err());
        assert_eq!(profile.interface_pack_bytes().unwrap(), bytes);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        let mut reads = 0;
        assert!(profile.load_interfaces(&spec, &cancelled, |_| {
            reads += 1; Ok(holding_pack(1e8, 2e4, "body.fsmatpk", 200_000.0))
        }).is_err());
        assert_eq!(reads, 0);
    });
}

fn held_environment_pack() -> NormalizedInterfacePack {
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, MoistureBasis, QuantityKind, SemanticType, ValueForm}};
    let base = holding_pack(1e8, 2e4, "body.fsmatpk", 200_000.0);
    let original = base.claims_pack().claims();
    let mut claims = ClaimSet::new();
    for id in original.observation_ids() {
        assert_eq!(claims.register_observation(original.observation(id).unwrap().clone()).unwrap(), id);
    }
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for (_, original) in original.claims_ordered() {
        let mut claim = original.clone();
        let PropertyValue::Scalar { value, dims } = claim.value else { unreachable!() };
        claim.value = PropertyValue::Curve { abscissa: "moisture_content".into(),
            abscissa_dims: Dims::NONE, knots: vec![(0.0, value), (0.2, 4.0*value)], dims };
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 280.0, 320.0)
            .with_quantity("omega", quantity(QuantityKind::Frequency(FrequencyConvention::Angular)), 0.0, 200_000.0)
            .with_quantity("moisture_content", quantity(QuantityKind::Moisture(MoistureBasis::DryMass)), 0.0, 0.2);
        claim.interpolation = InterpolationPolicy::LinearInside;
        claims.insert_claim(claim).unwrap();
    }
    let claims = NormalizedPack::new("synthetic-holding-environment", "synthetic-interface-v1",
        base.claims_pack().source_artifact(), "CC0-1.0", claims, vec![], vec![]).unwrap();
    NormalizedInterfacePack::new(base.card().surface_a().clone(), base.card().surface_b().clone(),
        base.card().context().clone(), claims).unwrap()
}

#[test]
fn g3_held_body_interface_state_changes_support_and_pressure_through_the_original_source() {
    scoped(|cx, gate| {
        // Manufactured state response, not measured tissue or sorption data.
        let source = held_environment_pack();
        let prepare = |moisture: Option<f64>| -> Result<Ready> {
            let mut spec = shell_spec();
            if let Some(value) = moisture {
                spec.material_states.insert("holding.fsintpk".into(),
                    vec![("moisture_content".into(), value)]);
            }
            let mut profile = held_profile();
            profile.load_interfaces(&spec, gate, |_| Ok(source.clone()))?;
            assert!(profile.has_interface_source("holding.fsintpk"));
            spec.compile_with_body(cx, gate, Some(profile), |name| Ok(pack(name, 1.0, None)))
        };
        let mut low = prepare(Some(0.0)).unwrap();
        let mut high = prepare(Some(0.2)).unwrap();
        let link = |ready: &Ready| ready.performance.bank().source().bank()
            .material_shell().unwrap().specimen.interfaces()[0].clone();
        let a = link(&low);
        let b = link(&high);
        assert_eq!(a.input().interface.card_identity(), b.input().interface.card_identity());
        assert_eq!(b.input().interface.query().axes()["moisture_content"], 0.2);
        assert_eq!(b.law().stiffness_n_m, 4.0*a.law().stiffness_n_m);
        assert_eq!(b.law().damping_n_s_m, 4.0*a.law().damping_n_s_m);
        for name in ["wire.fsmatpk", "sleeve.fsmatpk", "body.fsmatpk"] {
            assert_eq!(low.catalog[name].pack.content_hash(), high.catalog[name].pack.content_hash());
        }
        let low_pressure = pressure(&mut low, 7);
        let high_pressure = pressure(&mut high, 31);
        assert_ne!(low_pressure, high_pressure);
        assert_eq!(high_pressure, pressure(&mut prepare(Some(0.2)).unwrap(), 1));
        for ready in [&low, &high] {
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(prepare(None).is_err());
        assert!(prepare(Some(0.3)).is_err());
    });
}
