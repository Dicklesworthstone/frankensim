//! G3/G4: complete source data -> curved shell/string/gas mechanics -> pressure.
use super::*;
use super::super::shell::{MAX_MODES, MAX_NODES, MAX_SHELL_PROJECTION_WORK, MAX_TRIANGLES};

fn shell_text() -> &'static str { include_str!("../../arched.body") }
fn shell_profile() -> Profile { Profile::parse(shell_text().as_bytes()).unwrap() }
fn shell_spec() -> input::Spec {
    let bytes = include_str!("../../cylinder-flow.material-strings")
        .replace("48000 4096 0.01 128", "48000 96 0.01 128")
        .replace("cylinder-gas-loading 0.001 30 12", "cylinder-gas-loading 0.001 30 4");
    input::Spec::parse(bytes.as_bytes()).unwrap()
}
fn ready(cx: &Cx<'_>, gate: &CancelGate, loss: f64, shear: f64) -> Ready {
    shell_spec().compile_with_body(cx, gate, Some(shell_profile()),
        |name| Ok(if name == "body.fsmatpk" { pack_with_shear(name, loss, shear, None) }
            else { pack(name, 1.0, None) })).unwrap()
}

const SHELL_ROTATION: &str = "shell-rotational-bridge 0 0.8 0 10 0.2 0.3 0.5 0 1 0 0.001 1e-7 0";
fn shell_rotation_text() -> String {
    let spec = shell_spec();
    std::str::from_utf8(&spec.bytes).unwrap().replace("cylinder-gas-loading",
        &format!("rotational-bridges 1\n{SHELL_ROTATION}\ncylinder-gas-loading"))
}

#[test]
fn original_shell_rotation_drives_the_existing_torque_string_body_gas_and_pcm_solve() {
    scoped(|cx,gate| {
        let text = shell_rotation_text();
        let make = || input::Spec::parse(text.as_bytes()).unwrap().compile_with_body(cx,gate,
            Some(shell_profile()),|name|Ok(pack_with_shear(name,1.0,1.0,None))).unwrap();
        let mut connected = make();
        assert_eq!(connected.rotational_bridge_baselines_m,vec![None]);
        let source = connected.performance.bank().source().bank();
        let input::RotationalBody::Shell(geometry) = &connected.spec.rotational_bridges[0].body
            else { panic!("physical rotation input lost its original geometry"); };
        let actual = &source.rotational_body_connections()[0];
        assert_eq!(actual.body,source.shell_rotation(geometry,gate).unwrap());
        assert_eq!(actual.stiffness_n_m_per_rad,0.001);
        assert_eq!(actual.damping_n_m_s_per_rad,1e-7);
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.material_cavity().is_some());
        let mut unconnected = ready(cx,gate,1.0,1.0);
        let original = pressure(&mut unconnected,17);
        let samples = pressure(&mut connected,17);
        assert!(samples.iter().all(|v|v.is_finite()));
        assert!(samples.iter().zip(&original).any(|(a,b)|(a-b).abs()>1e-14));
        let frame = connected.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let mut partitioned = make();
        assert_eq!(pressure(&mut partitioned,37),samples);
    });
}

#[test]
fn physical_rotation_records_and_missing_or_fixed_geometry_refuse_before_pack_io() {
    let text = shell_rotation_text();
    for bad in [
        "shell-rotational-bridge 0 0 0 10 0.2 0.3 0.5 0 1 0 0.001 1e-7 0",
        "shell-rotational-bridge 0 0.8 0 10 -0.2 0.3 0.9 0 1 0 0.001 1e-7 0",
        "shell-rotational-bridge 0 0.8 0 10 0.2 0.3 0.4 0 1 0 0.001 1e-7 0",
        "shell-rotational-bridge 0 0.8 0 10 0.2 0.3 0.5 0 2 0 0.001 1e-7 0",
        "shell-rotational-bridge 0 0.8 0 10 0.2 0.3 0.5 0 1 0 0 0 0",
    ] {
        assert!(input::Spec::parse(text.replace(SHELL_ROTATION,bad).as_bytes()).is_err());
    }
    let legacy = rotational_input().replace(ROTATIONAL_BRIDGE_BLOCK,
        &format!("rotational-bridges 2\n{ROTATIONAL_BRIDGE_ROW}\n{SHELL_ROTATION}\n"));
    assert_eq!(input::Spec::parse(legacy.as_bytes()).unwrap().rotational_bridges.len(),2);
    scoped(|cx,gate| {
        for (row,profile) in [
            (SHELL_ROTATION,None),
            (SHELL_ROTATION,Some(super::profile())),
            ("shell-rotational-bridge 0 0.8 1 10 0.2 0.3 0.5 0 1 0 0.001 1e-7 0",Some(shell_profile())),
            ("shell-rotational-bridge 0 0.8 0 32 0.2 0.3 0.5 0 1 0 0.001 1e-7 0",Some(shell_profile())),
            ("shell-rotational-bridge 0 0.8 0 0 1 0 0 0 1 0 0.001 1e-7 0",Some(shell_profile())),
        ] {
            let spec = input::Spec::parse(text.replace(SHELL_ROTATION,row).as_bytes()).unwrap();
            let mut loads = 0;
            assert!(spec.compile_with_body(cx,gate,profile,|name| {
                loads += 1; Ok(pack(name,1.0,None))
            }).is_err());
            assert_eq!(loads,0);
        }
    });
}

#[test]
fn triangle_bridge_loads_the_same_material_modes_with_its_offset_moment() {
    scoped(|cx, gate| {
        let text = shell_text().replace("bridge 12 0 0 1",
            "bridge-triangle 10 0.2 0.3 0.5 0 0 1 0.02 0 0");
        let mut ready = shell_spec().compile_with_body(cx, gate,
            Some(Profile::parse(text.as_bytes()).unwrap()),
            |name| Ok(pack_with_shear(name, 1.0, 1.0, None))).unwrap();
        let source = ready.performance.bank().source().bank();
        let body = source.material_shell().unwrap();
        let image = source.material_shell_image().unwrap();
        let full = &body.force_ports[0];
        let mut expected = vec![0.0; 150];
        for (node, weight) in [(6,0.2),(7,0.3),(12,0.5)] {
            expected[6*node+2] = weight;
            expected[6*node+4] = -0.02*weight;
        }
        assert_eq!(full, &expected);
        for (mode, &projected) in image.physical_modes().iter().zip(&image.force_ports()[0]) {
            let direct: f64 = mode.iter().zip(full).map(|(u,f)|u*f).sum();
            assert!((direct-projected).abs() < 1e-12*direct.abs().max(1.0));
        }
        assert_eq!(ready.body_profile.as_ref().unwrap().bytes, text.as_bytes());
        let samples = pressure(&mut ready, 17);
        assert!(samples.iter().all(|v|v.is_finite()));
        assert!(samples.iter().any(|v|v.abs()>1e-14));
        let frame = ready.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn response_projection_reuses_playback_body_without_loading_string_sources() {
    scoped(|cx, gate| {
        let profile = shell_profile();
        let spec = shell_spec();
        let mut catalog = Catalog::new();
        let mut loaded = Vec::new();
        let projected = profile.project_shell(&spec, &mut catalog, &mut |name| {
            assert_eq!(name, "body.fsmatpk", "structural response requested a string source");
            loaded.push(name.to_owned());
            Ok(pack_with_shear(name, 1.0, 1.0, None).into())
        }, gate).unwrap();
        assert_eq!(loaded, ["body.fsmatpk"]);
        assert!(projected.assembly().is_none());
        let projected = projected.shell().unwrap();
        let performer = ready(cx, gate, 1.0, 1.0);
        let image = performer.performance.bank().source().bank().material_shell_image().unwrap();
        assert_eq!(projected.eigenvalue_intervals(), image.eigenvalue_intervals());
        assert_eq!(projected.physical_modes(), image.physical_modes());
        assert_eq!(projected.force_ports(), image.force_ports());
        assert_eq!(projected.volume_shapes_m2_per_sqrt_kg(), image.volume_shapes_m2_per_sqrt_kg());
        assert_eq!(projected.cavity_volume_shapes_m2_per_sqrt_kg(), image.cavity_volume_shapes_m2_per_sqrt_kg());
        assert_eq!(projected.damping().matrix_per_s(), image.damping().matrix_per_s());
    });
}

#[test]
fn loaded_analysis_keeps_the_playback_bank_without_a_time_propagator() {
    use fs_math::c64::C64;
    scoped(|cx, gate| {
        let mut spec = shell_spec();
        let budget = PreparedModalBudget {
            max_modes: 64, max_dense_bytes: 268_435_456,
            max_multiply_terms: 1_000_000_000_000,
        };
        spec.prepared_linear = Some(budget);
        let load = |name: &str| Ok(if name == "body.fsmatpk" {
            pack_with_shear(name, 1.0, 1.0, None)
        } else { pack(name, 1.0, None) });
        let bound = spec.clone().compile_bank_with_profiles(cx, gate,
            Some(shell_profile()), None, true, load).unwrap();
        let playback = spec.compile_with_body(cx, gate, Some(shell_profile()), load).unwrap();
        let source = bound.bank.source().bank();
        assert_eq!(source.preparation(), ModalPreparation::HarmonicOnly(budget));
        assert_eq!(source.prepared_linear_budget(), Some(budget));
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.material_cavity().is_some());
        assert_eq!(bound.bank.cards().len(), playback.performance.bank().cards().len());
        for (a, b) in bound.bank.cards().iter().zip(playback.performance.bank().cards()) {
            assert_eq!(a.fundamental_hz(), b.fundamental_hz());
            assert_eq!(a.linear_density_kg_m, b.linear_density_kg_m);
            assert_eq!(a.mode_count, b.mode_count);
        }
        let body = source.material_shell_image().unwrap();
        let other = playback.performance.bank().source().bank().material_shell_image().unwrap();
        assert_eq!(body.physical_modes(), other.physical_modes());
        assert_eq!(body.force_ports(), other.force_ports());
        let before: Vec<_> = bound.bank.system().components().iter()
            .map(|model| model.states().to_vec()).collect();
        for omega in [700.0, 2000.0, 6000.0] {
            let a = bound.bank.material_shell_bridge_harmonic_response(
                omega, 0, C64::from_re(0.001), 100_000_000, gate).unwrap();
            let b = playback.performance.bank().material_shell_bridge_harmonic_response(
                omega, 0, C64::from_re(0.001), 100_000_000, gate).unwrap();
            for (x, y) in a.released.dynamics.displacement_m_sqrt_kg.iter()
                .zip(&b.released.dynamics.displacement_m_sqrt_kg) {
                assert!((*x - *y).abs() <= 1e-10 * x.abs().max(y.abs()).max(1e-20));
            }
            let x = a.released.observer_pressure_pa;
            let y = b.released.observer_pressure_pa;
            assert!((x-y).abs() <= 1e-10 * x.abs().max(y.abs()).max(1e-20));
            assert_eq!(a.source_band_rad_s, b.source_band_rad_s);
            let body_q = a.body_displacement_m_sqrt_kg().unwrap();
            assert_eq!(body_q.len(), body.physical_modes().len());
            for (port, velocity) in body.force_ports().iter().zip(&a.bridge_velocity_m_s) {
                let displacement = port.iter().zip(body_q).fold(C64::ZERO,
                    |sum, (&shape, &q)| sum + q.scale(shape));
                let expected = C64::new(displacement.im * omega, -displacement.re * omega);
                assert!((expected - *velocity).abs() <= 1e-10 * velocity.abs().max(1e-20),
                    "body observation must borrow the original coupled shell coordinates");
            }
            let mut truncated = a.clone();
            truncated.released.dynamics.displacement_m_sqrt_kg.clear();
            assert!(truncated.body_displacement_m_sqrt_kg().is_err());
            let matrix = bound.bank.material_shell_bridge_harmonic_matrix(omega, 100_000_000, gate).unwrap();
            assert_eq!(matrix.len(), body.force_ports().len());
            let drives: Vec<_> = (0..matrix.len()).map(|index| C64::new(0.2, -0.1 * (index + 1) as f64)).collect();
            let coherent = bound.bank.material_shell_bridge_forces_harmonic_response(
                omega, &drives, 100_000_000, gate).unwrap();
            let mut expected_pressure = C64::ZERO;
            for (bridge, column) in matrix.iter().enumerate() {
                let single = bound.bank.material_shell_bridge_harmonic_response(
                    omega, bridge, C64::from_re(1.0), 100_000_000, gate).unwrap();
                assert_eq!(*column, single, "matrix columns must retain the single-drive physical path");
                expected_pressure = expected_pressure + column.released.observer_pressure_pa * drives[bridge];
            }
            assert!((coherent.released.observer_pressure_pa - expected_pressure).abs()
                <= 1e-10 * expected_pressure.abs().max(1e-20));
        }
        assert_eq!(before, bound.bank.system().components().iter()
            .map(|model| model.states().to_vec()).collect::<Vec<_>>());
        assert!(bound.into_performance(gate).is_err());
    });
}

#[test]
fn independent_bridge_axis_projects_the_original_shell_and_reaches_pressure() {
    scoped(|cx, gate| {
        let text = include_str!("../../cylinder-flow.material-strings")
            .replace("48000 4096 0.01 128", "48000 96 0.01 128")
            .replace("cylinder-gas-loading 0.001 30 12",
                "translational-bridges 1\ntranslational-bridge 0 0.8 0 1 0 11 1 0 0 1500 0.02 0\ncylinder-gas-loading 0.001 30 4");
        let build = |body_shear| input::Spec::parse(text.as_bytes()).unwrap()
            .compile_with_body(cx, gate, Some(shell_profile()), |name|
                Ok(if name == "body.fsmatpk" { pack_with_shear(name, 1.0, body_shear, None) }
                    else { pack(name, 1.0, None) })).unwrap();
        let mut directed = build(1.0);
        let source = directed.performance.bank().source().bank();
        let joints = source.translational_body_connections();
        assert_eq!(joints.len(), 1);
        assert_eq!(joints[0].direction.components(), [0.0, 1.0]);
        let image = source.material_shell_image().unwrap();
        let expected: Vec<_> = image.physical_modes().iter().map(|mode| mode[6*11]).collect();
        assert_eq!(joints[0].body.shapes, expected);
        assert!(expected.iter().any(|value| *value != 0.0));
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.material_cavity().is_some());
        let actual = pressure(&mut directed, 7);
        let mut same = build(1.0);
        assert_eq!(actual, pressure(&mut same, 31));
        let mut original = ready(cx, gate, 1.0, 1.0);
        assert_ne!(actual, pressure(&mut original, 7));
        let mut swapped = build(0.1);
        assert_ne!(joints_for(&same), joints_for(&swapped));
        assert_ne!(actual, pressure(&mut swapped, 13));
        let frame = directed.performance.bank().system().last_frame().unwrap();
        assert_eq!(frame.external_work_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

fn joints_for(ready: &Ready) -> Vec<f64> {
    ready.performance.bank().source().bank().translational_body_connections()[0].body.shapes.clone()
}

#[test]
fn independent_bridge_axes_require_real_shell_geometry_and_physical_directions() {
    let base = include_str!("../../directional-pluck.material-strings");
    let valid = "translational-bridges 1\ntranslational-bridge 0 0.8 0 1 0 11 1 0 0 1500 0.02 0\n";
    let text = base.replace("cavity neck", &format!("{valid}cavity neck"));
    for (old, new) in [
        ("0 0.8 0 1", "0 0.8 0 2"),
        ("0 11 1 0 0", "0 11 2 0 0"),
        ("1500 0.02 0", "0 0 0"),
        ("0 0.8 0 1", "0 1 0 1"),
        ("1500 0.02 0", "-1 0.02 0"),
    ] {
        assert!(input::Spec::parse(text.replace(old, new).as_bytes()).is_err(), "accepted {new}");
    }
    scoped(|cx, gate| {
        let result = input::Spec::parse(text.as_bytes()).unwrap().compile(cx, gate,
            |name| Ok(pack(name, 1.0, None)));
        assert!(result.is_err(), "a plate profile cannot supply an unmodeled lateral shell displacement");
    });
}

#[test]
fn curved_material_body_keeps_layers_planes_torsion_cavity_and_gas_in_actual_pressure() {
    scoped(|cx, gate| {
        let mut a = ready(cx, gate, 1.0, 1.0);
        let source = a.performance.bank().source().bank();
        let shell = source.material_shell().unwrap();
        let image = source.material_shell_image().unwrap();
        assert_eq!(a.body_profile.as_ref().unwrap().bytes, shell_text().as_bytes());
        assert_eq!(shell.specimen.mesh().nodes.len(), 25);
        assert_eq!(shell.specimen.mesh().tris.len(), 32);
        assert_eq!(shell.specimen.mesh().nodes[12], [0.2, 0.15, 0.006]);
        assert_eq!(shell.specimen.regions()[0].input().axis1_global, [1.0, 0.0, 0.0]);
        assert!(shell.specimen.mass_kg() > 0.0);
        assert_eq!(shell.force_ports.len(), 1);
        assert_eq!(shell.force_ports[0].len(), 150);
        assert_eq!(shell.force_ports[0][6*12+2], 1.0);
        assert_eq!(image.physical_modes().len(), 2);
        assert!(image.physical_modes().iter().all(|v| v.len() == 150));
        assert!(image.volume_shapes_m2_per_sqrt_kg().iter().any(|v| *v != 0.0));
        assert_eq!(image.damping().mode_count(), 2);
        let MaterialStringSource::WithPolarizations(planes) = &source.strings()[0].source else {
            panic!("two physical string planes disappeared");
        };
        assert!(matches!(planes.source(), MaterialStringSource::WithTorsion(_)));
        assert_eq!(source.cylinder_gas_loading().unwrap().strings()[0].flow().branches().len(), 4);
        assert_eq!(source.cylinder_gas_loading().unwrap().strings()[0].flow().gas(), source.gas());
        assert_eq!(source.material_cavity().unwrap().cavity.gas(), source.gas());
        assert!(a.performance.preload_energy_j().unwrap() > 0.0);

        let samples = pressure(&mut a, 7);
        assert!(samples.iter().any(|p| p.abs() > 1e-14));
        assert!(a.performance.bank().system().components()[1].states().iter()
            .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        let frame = a.performance.bank().system().last_frame().unwrap();
        assert_eq!(frame.external_work_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);

        let mut b = ready(cx, gate, 1.0, 1.0);
        let encoded = fs_couple::pcm_wav::encode_pcm16_wav(&samples, b.spec.rate, b.spec.full_scale).unwrap().0;
        let (wav, report) = b.write(std::io::Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), encoded);
        assert_eq!(report.samples, 96);

        // Change only sourced transverse shear, then only sourced viscosity.
        // The original material string, its gas and its installed load persist.
        let mut softer = ready(cx, gate, 1.0, 0.1);
        let mut lossy = ready(cx, gate, 4.0, 1.0);
        assert_eq!(a.performance.bank().cards()[0].fundamental_hz(),
            softer.performance.bank().cards()[0].fundamental_hz());
        for name in ["wire.fsmatpk", "sleeve.fsmatpk"] {
            assert_eq!(a.catalog[name].pack.content_hash(), softer.catalog[name].pack.content_hash());
            assert_eq!(a.catalog[name].pack.content_hash(), lossy.catalog[name].pack.content_hash());
        }
        let modes = |r: &Ready| r.performance.bank().system().components()[1].modes().to_vec();
        assert_ne!(modes(&a), modes(&softer));
        assert_eq!(modes(&a), modes(&lossy));
        assert_ne!(a.performance.bank().source().bank().material_shell_image().unwrap()
            .damping().matrix_per_s(), lossy.performance.bank().source().bank()
            .material_shell_image().unwrap().damping().matrix_per_s());
        assert!(samples.iter().zip(pressure(&mut softer, 31)).any(|(a, b)| *a != b));
        assert!(samples.iter().zip(pressure(&mut lossy, 31)).any(|(a, b)| *a != b));
    });
}

fn held_shell_pack(band: (f64, f64)) -> NormalizedPack {
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, MoistureBasis, QuantityKind, SemanticType, ValueForm}};
    // Manufactured conditioning response, not a wood measurement or sorption law.
    let original = pack_with_shear("body.fsmatpk", 1.0, 1.0, None);
    let mut claims = ClaimSet::new();
    for id in original.claims().observation_ids() {
        assert_eq!(claims.register_observation(original.claims().observation(id).unwrap().clone()).unwrap(), id);
    }
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for (_, original) in original.claims().claims_ordered() {
        let mut claim = original.clone();
        let PropertyValue::Scalar { value, dims } = claim.value else { unreachable!() };
        let scale = if dims == Density::DIMS { 1.2 }
            else if dims == Pressure::DIMS { 0.7 }
            else if dims == DynViscosity::DIMS { 2.0 }
            else { 1.0 };
        claim.value = PropertyValue::Curve { abscissa: "moisture_content".into(),
            abscissa_dims: Dims::NONE, knots: vec![(0.0, value), (0.2, scale*value)], dims };
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 280.0, 320.0)
            .with_quantity("omega", quantity(QuantityKind::Frequency(FrequencyConvention::Angular)), band.0, band.1)
            .with_quantity("moisture_content", quantity(QuantityKind::Moisture(MoistureBasis::DryMass)), 0.0, 0.2);
        claim.interpolation = InterpolationPolicy::LinearInside;
        claims.insert_claim(claim).unwrap();
    }
    NormalizedPack::new("synthetic-shell-environment", "synthetic-tensor-v1",
        original.source_artifact(), "CC0-1.0", claims, vec![], vec![]).unwrap()
}

#[test]
fn held_body_material_state_reaches_mass_stiffness_full_loss_and_pressure() {
    scoped(|cx, gate| {
        let source = held_shell_pack((0.0, 200_000.0));
        let prepare = |source: &NormalizedPack, moisture: Option<f64>| -> Result<Ready> {
            let mut spec = shell_spec();
            if let Some(value) = moisture {
                spec.material_states.insert("body.fsmatpk".into(),
                    vec![("moisture_content".into(), value)]);
            }
            spec.compile_with_body(cx, gate, Some(shell_profile()), |name|
                Ok(if name == "body.fsmatpk" { source.clone() } else { pack(name, 1.0, None) }))
        };
        let mut low = prepare(&source, Some(0.0)).unwrap();
        let mut high = prepare(&source, Some(0.2)).unwrap();
        let a = low.performance.bank().source().bank().material_shell().unwrap().specimen.regions()[0].clone();
        let b = high.performance.bank().source().bank().material_shell().unwrap().specimen.regions()[0].clone();
        assert_eq!(a.area_m2(), b.area_m2());
        assert_eq!(a.thickness_m(), b.thickness_m());
        assert!((b.mass_kg()/a.mass_kg() - 1.2).abs() < 1e-12);
        assert_eq!(a.input().region.material.card_identity(), b.input().region.material.card_identity());
        assert_eq!(b.input().region.material.state_coordinate("moisture_content"), Some(0.2));
        for (key, scale) in [("density", 1.2), ("young_modulus_1", 0.7),
            ("shear_modulus_23", 0.7), (VISCOSITY[0], 2.0),
            (fs_couple::shear_shell::TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES[1], 2.0)] {
            let scalar = |region: &fs_couple::shear_shell::ResolvedShellRegion|
                region.input().region.material.property(key).unwrap().value_si();
            assert!((scalar(&b)/scalar(&a) - scale).abs() < 1e-12);
        }
        let modes = |r: &Ready| r.performance.bank().system().components()[1].modes().to_vec();
        assert_ne!(modes(&low), modes(&high));
        assert_ne!(low.performance.bank().source().bank().material_shell_image().unwrap().damping().matrix_per_s(),
            high.performance.bank().source().bank().material_shell_image().unwrap().damping().matrix_per_s());
        for name in ["wire.fsmatpk", "sleeve.fsmatpk", "body.fsmatpk"] {
            assert_eq!(low.catalog[name].pack.content_hash(), high.catalog[name].pack.content_hash());
        }
        assert_eq!(low.performance.bank().cards()[0].fundamental_hz(), high.performance.bank().cards()[0].fundamental_hz());
        let low_pressure = pressure(&mut low, 7);
        let high_pressure = pressure(&mut high, 31);
        assert_ne!(low_pressure, high_pressure);
        assert_eq!(high_pressure, pressure(&mut prepare(&source, Some(0.2)).unwrap(), 1));
        for ready in [&low, &high] {
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(prepare(&source, None).is_err(), "no default moisture state");
        assert!(prepare(&source, Some(0.3)).is_err(), "no source extrapolation");
        // The reference query is in-band; actual shell modes still need support.
        assert!(prepare(&held_shell_pack((999.999, 1000.001)), Some(0.1)).is_err());
    });
}

#[test]
fn shell_geometry_and_resource_records_refuse_without_borrowing_missing_material_data() {
    let over_nodes = format!("nodes {}", MAX_NODES + 1);
    let over_triangles = format!("triangles {}", MAX_TRIANGLES + 1);
    let over_work = format!("projection-work {}", MAX_SHELL_PROJECTION_WORK + 1);
    for (from, to) in [
        ("nodes 25", over_nodes.as_str()),
        ("triangles 32", over_triangles.as_str()),
        ("mass consistent", "mass automatic"),
        ("director-min-dot 0.9", "director-min-dot 0"),
        ("drilling 0.001", "drilling 0"),
        ("projection-work 67108864", over_work.as_str()),
        ("orthotropic12 1 0.17", "orthotropic12 0 0.17"),
        ("0.8333333333333334 1 0 0", "0.8333333333333334 2 0 0"),
        ("region body @body", "region body ../body.fsmatpk"),
        ("triangle 0 1 6 body", "triangle 0 1 6 unknown"),
        ("support 0 111111", "support 0 11111"),
        ("bridge 12 0 0 1", "bridge 0 0 0 1"),
        ("bridge 12 0 0 1", "bridge 12 0 0 2"),
    ] {
        assert!(Profile::parse(shell_text().replace(from, to).as_bytes()).is_err(), "accepted {to}");
    }
    assert!(Profile::parse(format!("{}extra\n", shell_text()).as_bytes()).is_err());
    scoped(|cx, gate| {
        // Topology, actual physical directors and bridge membership are checked
        // before the first source pack is requested by the executable.
        for text in [
            shell_text().replace("triangle 0 1 6 body", "triangle 0 0 6 body"),
            shell_text().replace("triangle 0 1 6 body", "triangle 0 6 1 body"),
            shell_text().replace("bridges 1\nbridge 12 0 0 1", "bridges 2\nbridge 12 0 0 1\nbridge 11 0 0 1"),
        ] {
            let mut loads = 0;
            assert!(shell_spec().compile_with_body(cx, gate, Some(Profile::parse(text.as_bytes()).unwrap()),
                |name| { loads += 1; Ok(pack(name, 1.0, None)) }).is_err());
            assert_eq!(loads, 0);
        }
        for key in ["shear_modulus_13", "transverse_shear_viscosity_23"] {
            assert!(shell_spec().compile_with_body(cx, gate, Some(shell_profile()),
                |name| Ok(pack(name, 1.0, Some(key)))).is_err(), "borrowed missing {key}");
        }
        let too_small = shell_text().replace("projection-work 67108864", "projection-work 1");
        assert!(shell_spec().compile_with_body(cx, gate,
            Some(Profile::parse(too_small.as_bytes()).unwrap()),
            |name| Ok(pack(name, 1.0, None))).is_err());
    });
}

#[test]
fn shell_rotational_chord_uses_actual_3d_nodes_and_common_force_direction() {
    let mut spec = shell_spec();
    spec.strings.push(spec.strings[0].clone());
    let text = shell_text().replace("bridges 1\nbridge 12 0 0 1",
        "bridges 2\nbridge 12 0 0 1\nbridge 11 0 0 1");
    let profile = Profile::parse(text.as_bytes()).unwrap();
    let lever = profile.bridge_separation_m(&spec, 0, 1).unwrap();
    assert!((lever - 0.1).abs() < 1e-15);
    assert!(lever < fs_math::det::hypot(0.1, 0.0015));
    // Both points at equal height with a force direction parallel to their
    // baseline measure extension, not a work-conjugate angular coordinate.
    let parallel = text.replace("bridge 12 0 0 1", "bridge 13 1 0 0")
        .replace("bridge 11 0 0 1", "bridge 11 1 0 0");
    let parallel = Profile::parse(parallel.as_bytes()).unwrap();
    assert!(parallel.bridge_separation_m(&spec, 0, 1).is_err());
    let nearly_parallel = text.replace("bridge 12 0 0 1", "bridge 13 1 0 1e-16")
        .replace("bridge 11 0 0 1", "bridge 11 1 0 1e-16");
    let nearly_parallel = Profile::parse(nearly_parallel.as_bytes()).unwrap();
    assert!(nearly_parallel.bridge_separation_m(&spec, 0, 1).is_err());
    let different = Profile::parse(text.replace("bridge 11 0 0 1", "bridge 11 1 0 0").as_bytes()).unwrap();
    assert!(different.bridge_separation_m(&spec, 0, 1).is_err());
    let coincident = Profile::parse(text.replace("bridge 11", "bridge 12").as_bytes()).unwrap();
    assert!(coincident.bridge_separation_m(&spec, 0, 1).is_err());
}

#[path = "interface_tests.rs"]
mod interfaces;

#[path = "pressure_tests.rs"]
mod pressure_surfaces;

/// A flat `nx * ny` grid body in the layout of `arched.body`: `nodes` rows, `2 (nx-1)(ny-1) + extra`
/// triangle rows (the extras repeat triangle 0), one support and one bridge. Only counts and indices are
/// meaningful: this feeds the parser, not the geometry validator.
fn grid_body(nx: usize, ny: usize, declared_nodes: usize, extra_triangles: usize) -> String {
    let mut text: String = shell_text().lines().take(8).map(|l| format!("{l}\n")).collect();
    text.push_str(&format!("nodes {declared_nodes}\n"));
    for j in 0..ny {
        for i in 0..nx {
            text.push_str(&format!("node {} {} 0 0 0 1\n",
                i as f64 / (nx - 1) as f64, j as f64 / (ny - 1) as f64));
        }
    }
    let mut triangles = Vec::new();
    for j in 0..ny - 1 {
        for i in 0..nx - 1 {
            let n = j * nx + i;
            triangles.push([n, n + 1, n + nx + 1]);
            triangles.push([n, n + nx + 1, n + nx]);
        }
    }
    for _ in 0..extra_triangles { triangles.push(triangles[0]); }
    text.push_str(&format!("triangles {}\n", triangles.len()));
    for t in &triangles { text.push_str(&format!("triangle {} {} {} body\n", t[0], t[1], t[2])); }
    text.push_str(&format!("supports 1\nsupport 0 111111\nbridges 1\nbridge {} 0 0 1\n", nx * ny - 1));
    text
}

#[test]
fn shell_profile_admits_exactly_the_node_ceiling_and_refuses_one_more() {
    assert_eq!(MAX_NODES % 32, 0, "the grid below is 32 rows high");
    let nx = MAX_NODES / 32;
    assert!(Profile::parse(grid_body(nx, 32, MAX_NODES, 0).as_bytes()).is_ok());
    // One more node is refused when it is only declared and when it is really present.
    assert!(Profile::parse(grid_body(nx, 32, MAX_NODES + 1, 0).as_bytes()).is_err());
    // The real over-cap grid is `columns` wide, the smallest factor of MAX_NODES + 1.
    let columns = (2..MAX_NODES + 1).find(|f| (MAX_NODES + 1) % f == 0).expect("MAX_NODES + 1 must not be prime");
    assert!(Profile::parse(grid_body(columns, (MAX_NODES + 1) / columns, MAX_NODES + 1, 0).as_bytes()).is_err());
}

#[test]
fn shell_profile_admits_exactly_the_triangle_ceiling_and_refuses_one_more() {
    let nx = MAX_NODES / 32;
    let base = 2 * (nx - 1) * 31;
    assert!(base <= MAX_TRIANGLES, "the node-cap grid must leave room for the triangle cap");
    assert!(Profile::parse(grid_body(nx, 32, MAX_NODES, MAX_TRIANGLES - base).as_bytes()).is_ok());
    assert!(Profile::parse(grid_body(nx, 32, MAX_NODES, MAX_TRIANGLES - base + 1).as_bytes()).is_err());
}

#[test]
fn shell_profile_admits_exactly_the_mode_ceiling_and_refuses_one_more() {
    let profile = shell_profile();
    let mut spec = shell_spec();
    spec.body.modes = MAX_MODES;
    assert!(profile.preflight(&spec).is_ok());
    spec.body.modes = MAX_MODES + 1;
    assert!(profile.preflight(&spec).is_err());
    spec.body.modes = 0;
    assert!(profile.preflight(&spec).is_err());
    // A flat plate still admits no more than its own 16.
    let plate = Profile::parse(include_str!("../../regional.body").as_bytes()).unwrap();
    spec.body.modes = input::MAX_PLATE_BODY_MODES;
    assert!(plate.preflight(&spec).is_ok());
    spec.body.modes = input::MAX_PLATE_BODY_MODES + 1;
    assert!(plate.preflight(&spec).is_err());
}

/// A complete instrument admits exactly `MAX_TOTAL_MODES` modal coordinates (1024, the prepared owner's ceiling) and
/// refuses one more; at that ceiling it admits exactly `MAX_WORK` setup terms, one link per body mode at the body
/// ceiling, and refuses one more link. A `prepared-linear` record admits the library ceiling and refuses one more mode.
#[test]
fn the_complete_instrument_admits_exactly_the_mode_and_work_ceilings_and_refuses_one_more() {
    let spec = shell_spec();
    let dimensions = |links, modes| spec.network_dimensions_with_body_modes(false, links, modes);
    let (base, base_links, _) = dimensions(0, 0).unwrap();
    let (at, _, _) = dimensions(0, input::MAX_TOTAL_MODES - base).unwrap();
    assert_eq!(at, input::MAX_TOTAL_MODES);
    let over = dimensions(0, input::MAX_TOTAL_MODES - base + 1).unwrap_err();
    assert!(over.to_string().contains("exceeds 1024 modes"), "{over}");
    let extra = input::MAX_BODY_MODES - base_links;
    let (_, links, terms) = dimensions(extra, input::MAX_TOTAL_MODES - base).unwrap();
    assert_eq!(links, input::MAX_BODY_MODES);
    assert_eq!(terms as u64, input::MAX_WORK);
    let over = dimensions(extra + 1, input::MAX_TOTAL_MODES - base).unwrap_err();
    assert!(over.to_string().contains("setup work budget"), "{over}");
    let body = "body body.fsmatpk 0.4 0.3 0.001 2 1";
    let cavity = "cavity neck 0.003 0.01 0.02 100 1";
    let text = violin_performance(96, body, cavity, CANNONE_DIRECTIONS);
    let ceiling = fs_couple::render::schedule::force::coupled::MAX_PREPARED_MODAL_MODES;
    let record = |modes: usize| format!("prepared-linear {modes} ");
    let with = |modes: usize| text.replacen(&record(512), &record(modes), 1);
    assert!(input::Spec::parse(with(ceiling).as_bytes()).is_ok());
    assert!(input::Spec::parse(with(ceiling + 1).as_bytes()).is_err());
}

#[test]
fn shell_profile_admits_exactly_the_byte_ceiling_and_refuses_one_more() {
    let text = shell_text();
    // Pad one numeric token so the whole profile is exactly the ceiling, then one byte over.
    let padded = |extra: usize| text.replace("drilling 0.001", &format!("drilling 0.001{}", "0".repeat(extra)));
    let to_cap = super::super::MAX_SHELL_BYTES - text.len();
    assert_eq!(padded(to_cap).len(), super::super::MAX_SHELL_BYTES);
    assert!(Profile::parse(padded(to_cap).as_bytes()).is_ok());
    assert!(Profile::parse(padded(to_cap + 1).as_bytes()).is_err());
}

#[test]
fn shell_profile_admits_exactly_the_projection_work_ceiling_and_refuses_one_more() {
    let with = |work: usize| shell_text().replace("projection-work 67108864", &format!("projection-work {work}"));
    assert!(Profile::parse(with(MAX_SHELL_PROJECTION_WORK).as_bytes()).is_ok());
    assert!(Profile::parse(with(MAX_SHELL_PROJECTION_WORK + 1).as_bytes()).is_err());
    // A cap-sized component's own projection work fits the ceiling.
    let cost = fs_plate::shell::mitc::viscosity::projection_terms(MAX_NODES, MAX_TRIANGLES, MAX_MODES).unwrap();
    assert!(cost < MAX_SHELL_PROJECTION_WORK, "{cost} terms");
}

/// A flat clamped `nx * ny` plate in the layout of `arched.body`, with one bridge at its centre node.
/// Unlike `grid_body` this is a complete, geometrically valid shell: every boundary node is clamped.
fn clamped_plate_body(nx: usize, ny: usize, window: &str, work: usize) -> String {
    let mut text: String = shell_text().lines().take(8).map(|l| format!("{l}\n")).collect();
    text = text.replace("window 1 100000000", window)
        .replace("projection-work 67108864", &format!("projection-work {work}"));
    text.push_str(&format!("nodes {}\n", nx * ny));
    for j in 0..ny {
        for i in 0..nx {
            text.push_str(&format!("node {} {} 0 0 0 1\n",
                i as f64 / (nx - 1) as f64, j as f64 / (ny - 1) as f64));
        }
    }
    let mut triangles = Vec::new();
    for j in 0..ny - 1 {
        for i in 0..nx - 1 {
            let n = j * nx + i;
            triangles.push([n, n + 1, n + nx + 1]);
            triangles.push([n, n + nx + 1, n + nx]);
        }
    }
    text.push_str(&format!("triangles {}\n", triangles.len()));
    for t in &triangles { text.push_str(&format!("triangle {} {} {} body\n", t[0], t[1], t[2])); }
    let clamped: Vec<usize> = (0..nx * ny)
        .filter(|n| { let (i, j) = (n % nx, n / nx); i == 0 || j == 0 || i == nx - 1 || j == ny - 1 }).collect();
    text.push_str(&format!("supports {}\n", clamped.len()));
    for n in &clamped { text.push_str(&format!("support {n} 111111\n")); }
    text.push_str(&format!("bridges 1\nbridge {} 0 0 1\n", (ny / 2) * nx + nx / 2));
    text
}

/// Rendered samples of the end-to-end ceiling test. A COST choice, not a limit: nothing caps the sample count. Measured on ts1
/// (release, load average about 60, 98-coordinate prepared network, 2026-10-23, bead .358): construction takes 16 s and rendering
/// about 1.0 ms per sample; before the pre-summed loss (A3) and the one-integral construction they were 115 s and 19.5 ms. 96 is
/// the length the curved-body tests in this file already render, and the length the route's before/after pressure-bit checks
/// compare. Every dimension the ceilings bound (nodes, triangles, modes, gas memory, projection work) is at full value.
const CEILING_SAMPLES: u64 = 96;

/// The ceiling fixture: the clamped 15 x 11 plate with `modes` retained modes under the cylinder-flow performance (its 12 gas-memory
/// intervals) for `samples` samples, with the prepared-linear record that any network above 42 modal coordinates needs.
fn ceiling_inputs(modes: usize, samples: u64) -> (Profile, input::Spec) {
    let text = clamped_plate_body(15, 11, "window 1 400000000", 2_000_000_000);
    let profile = Profile::parse(text.as_bytes()).unwrap();
    // 70 body modes, the cavity, two string planes, twist and gas memory (3 modes x 12 intervals) make 110: declare room for 128.
    let performance = include_str!("../../cylinder-flow.material-strings")
        .replace("48000 4096 0.01 128", &format!("48000 {samples} 0.01 128"))
        .replace("cylinder-gas-loading 0.001 30 12",
            "cylinder-gas-loading 0.001 30 12\nprepared-linear 128 1073741824 10000000000000");
    let mut spec = input::Spec::parse(performance.as_bytes()).unwrap();
    spec.body.modes = modes;
    (profile, spec)
}

#[test]
fn a_mesh_and_mode_count_beyond_every_former_ceiling_run_through_the_whole_pipeline() {
    // 15 x 11 = 165 nodes (former ceiling 128), 280 triangles (256), 70 retained modes (16), and a projection
    // allowance of 2e9 terms (268,435,456). The window holds well over 70 bending modes of the fixture plate.
    let (nx, ny) = (15, 11);
    let (profile, spec) = ceiling_inputs(70, CEILING_SAMPLES);
    scoped(|cx, gate| {
        let constructing = std::time::Instant::now();
        let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(if name == "body.fsmatpk" {
            pack_with_shear(name, 1.0, 1.0, None) } else { pack(name, 1.0, None) })).unwrap();
        // Timing lines go to stderr (visible with --nocapture): this test is slow by the cost measured above.
        eprintln!("ceiling e2e: {}", cost_context(&r.spec));
        eprintln!("ceiling e2e: constructed in {:.1?}", constructing.elapsed());
        let source = r.performance.bank().source().bank();
        let image = source.material_shell_image().unwrap();
        assert_eq!(source.material_shell().unwrap().specimen.mesh().nodes.len(), nx * ny);
        assert_eq!(source.material_shell().unwrap().specimen.mesh().tris.len(), 2 * (nx - 1) * (ny - 1));
        assert_eq!(image.physical_modes().len(), 70);
        assert!(image.physical_modes().iter().all(|v| v.len() == 6 * nx * ny));
        assert_eq!(image.damping().mode_count(), 70);
        let rendering = std::time::Instant::now();
        let samples = pressure(&mut r, 7);
        eprintln!("ceiling e2e: {} samples in {:.1?}", samples.len(), rendering.elapsed());
        assert!(samples.iter().all(|p| p.is_finite()));
        assert!(samples.iter().any(|p| p.abs() > 1e-14), "the body radiated nothing");
        eprintln!("ceiling e2e: pressure bits fnv1a {:#018x}", pressure_bits(&samples));
    });
}

/// FNV-1a over the exact bits of every pressure sample: an identity of the rendered sound, not a tolerance. It is printed,
/// not pinned: an exact restructuring of the route is checked by comparing it before and after the change.
fn pressure_bits(samples: &[f64]) -> u64 {
    samples.iter().flat_map(|p| p.to_bits().to_le_bytes())
        .fold(0xcbf2_9ce4_8422_2325_u64, |hash, byte| (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3))
}

/// What a wall time was measured on: host, logical CPUs, the load average when the line is printed, the render's thread,
/// and the complete network the prepared route propagates (coordinates, connections).
fn cost_context(spec: &input::Spec) -> String {
    let host = std::fs::read_to_string("/etc/hostname").map_or_else(|_| "unknown".into(), |s| s.trim().to_string());
    let cpus = std::thread::available_parallelism().map_or(0, |n| n.get());
    let load = std::fs::read_to_string("/proc/loadavg")
        .map_or_else(|_| "unavailable".into(), |s| s.split_whitespace().take(3).collect::<Vec<_>>().join(" "));
    let (coordinates, links, _) = spec.network_dimensions(true).unwrap();
    format!("host {host}, {cpus} logical CPUs, load average {load}, render on the test thread, \
        {coordinates} coordinates, {links} connections")
}

/// The ceiling fixture rendered for half a second (24000 samples at 48 kHz) through the prepared exact route: the wall time
/// per audio second with what it was measured on, and the pressure bits of the whole half second.
#[test]
#[ignore = "half-second ceiling render: run with --ignored --nocapture"]
fn the_half_second_ceiling_render_reports_its_cost_and_pressure_bits() {
    let (profile, spec) = ceiling_inputs(70, 24_000);
    scoped(|cx, gate| {
        eprintln!("half-second ceiling: {}", cost_context(&spec));
        let constructing = std::time::Instant::now();
        let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(if name == "body.fsmatpk" {
            pack_with_shear(name, 1.0, 1.0, None) } else { pack(name, 1.0, None) })).unwrap();
        let built = constructing.elapsed();
        let rendering = std::time::Instant::now();
        let samples = pressure(&mut r, 7);
        let rendered = rendering.elapsed();
        assert_eq!(samples.len(), 24_000);
        assert!(samples.iter().all(|p| p.is_finite()));
        let audio_s = samples.len() as f64 / 48_000.0;
        eprintln!("half-second ceiling: constructed in {built:.1?}, {} samples rendered in {rendered:.1?} \
            ({:.1} s render per audio second, {:.1} s with construction)", samples.len(),
            rendered.as_secs_f64() / audio_s, (built + rendered).as_secs_f64() / audio_s);
        eprintln!("half-second ceiling: {}", cost_context(&r.spec));
        eprintln!("half-second ceiling: pressure bits fnv1a {:#018x}", pressure_bits(&samples));
    });
}

/// The four strings of the violin-sized rung, G3 D4 A4 E5: pitch [Hz] and radius [m] of plain strings in the synthetic steel
/// of these tests (7800 kg/m^3, E 200 GPa, G 80 GPa), NOT violin strings (their packs are bead .192; a wound string here
/// would refuse torsion). The tension that puts each speaking length on its pitch is derived below, not chosen: 36.5, 56.8,
/// 81.7 and 77.5 N, 252 N in all.
const VIOLIN_STRINGS: [(f64, f64); 4] = [(196.0, 0.000_3), (293.66, 0.000_25), (440.0, 0.000_2), (659.26, 0.000_13)];

/// Every string keeps every bending and twist mode of its full span below this frequency [Hz].
const VIOLIN_STRING_BAND_HZ: f64 = 10_000.0;

/// The violin-sized cost rung of bead .358, on synthetic data throughout (NOT a measured violin; the il Cannone body is
/// .65/.153). It uses the clamped 15 x 11 plate with 128 retained modes and one bridge node per string, a Helmholtz cavity
/// (2 L; the two f-holes as one 18.7 mm neck through 3 mm), and the four strings above. Each string has a 0.328 m speaking
/// length from the nut to its bridge foot and a 55 mm afterlength to the tailpiece (both ends fixed), two planes and its
/// twist, and every mode of that full span below `VIOLIN_STRING_BAND_HZ`. String 0 (G3) is bowed at 0.2 m/s and 0.5 N, a
/// tenth of the speaking length from the bridge. There is no cylinder gas loading: its memory (bending modes x planes x
/// intervals) would exceed the 512-coordinate prepared ceiling.
fn violin_inputs(samples: u64) -> (Profile, input::Spec) {
    let (nx, ny) = (15, 11);
    let centre = (ny / 2) * nx + nx / 2;
    let feet: String = [centre - 2, centre - 1, centre + 1, centre + 2].iter()
        .map(|node| format!("bridge {node} 0 0 1\n")).collect();
    let text = clamped_plate_body(nx, ny, "window 1 4000000000", 2_000_000_000)
        .replace(&format!("bridges 1\nbridge {centre} 0 0 1\n"), &format!("bridges 4\n{feet}"));
    let profile = Profile::parse(text.as_bytes()).unwrap();
    let performance = violin_performance(samples, "body body.fsmatpk 0.4 0.3 0.001 128 1",
        "cavity neck 0.002 0.0187 0.003 100 1", "0.6 0.8 1 0");
    (profile, input::Spec::parse(performance.as_bytes()).unwrap())
}

/// The violin-sized rung's performance (see `violin_inputs`: its strings, bow and schedule) around a body and a cavity
/// record, with every string's `polarizations two` directions (BOW_X BOW_Y BRIDGE_X BRIDGE_Y).
fn violin_performance(samples: u64, body: &str, cavity: &str, directions: &str) -> String {
    use std::f64::consts::PI;
    let (speaking, span) = (0.328_f64, 0.383_f64);
    let bridge = speaking / span;
    let (density, young, shear) = (7800.0, 200e9, 80e9);
    let mut strings = String::new();
    for &(pitch, radius) in &VIOLIN_STRINGS {
        // The bridge foot is the speaking length's end: f0 = sqrt(T / mu) / (2 L_speaking), mu = rho pi r^2 (bonded steel
        // core and steel sleeve).
        let tension = density * PI * radius * radius * (2.0 * speaking * pitch).powi(2);
        // The material string's modes are those of its full span: f_n = n c / (2 L) sqrt(1 + B n^2) with c = 2 L_speaking f0
        // and B = pi^2 E I / (T L^2); torsion f_m = m sqrt(G / rho) / (2 L).
        let fundamental = 2.0 * speaking * pitch / (2.0 * span);
        let inharmonicity = PI * PI * young * PI * radius.powi(4) / 4.0 / (tension * span * span);
        let bending = (1..=input::MAX_MODES)
            .take_while(|&n| n as f64 * fundamental * (1.0 + inharmonicity * (n * n) as f64).sqrt() < VIOLIN_STRING_BAND_HZ)
            .count();
        let twist = (VIOLIN_STRING_BAND_HZ * 2.0 * span / (shear / density).sqrt()) as usize;
        strings.push_str(&format!("string wire.fsmatpk {span} {radius} tension {tension} {bending} {twist} {bridge} 3000 0.1\n\
            section bonded {} 1\nlayer sleeve.fsmatpk {radius}\npolarizations two {directions}\n", radius * 2.0 / 3.0));
    }
    let bow = bridge * 0.9;
    format!("frankensim-material-strings-v3\naudio 48000 {samples} 0.01 128\nstate 293.15 1000\n\
        gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2\n{body}\nstrings 4\n{strings}\
        prepared-linear 512 4294967296 1000000000000000\n{cavity}\n\
        friction stribeck 0.8 0.4 0.1 0\nforces retain 0 0\nschedule\nfrankensim-gesture-schedule-v1\n\
        control_rate_hz\t700\ntracks\t1\ntrack\tbow\tbow\t0\ninitial\tbow\t2e-1\t5e-1\t{bow:e}\nevents\t0\n")
}

/// Every il Cannone string's BOW_X BOW_Y BRIDGE_X BRIDGE_Y. The bridge nodes push along x, across the body, which is also
/// the bow's tangent, so the bow drives the plane the bridge holds. The bow-normal plane has no bridge (the string record
/// carries one bridge direction) and nothing drives it. With the violin-sized rung's `0.6 0.8 1 0`, 0.8 of the bow drove
/// that free plane, and the G string sounded its full 0.383 m span (167 Hz) instead of its 0.328 m speaking length.
const CANNONE_DIRECTIONS: &str = "1 0 1 0";

/// The il Cannone body's assembly choices (`cannone-h8.assembly`, the file `material_strings assemble-shell` reads): the
/// window from 15 Hz to 2 kHz (the packs declare 100 Hz to 5 kHz and 150 Hz to 10 kHz), spruce as the performance body
/// pack and maple the other, the 3 mm reference thickness, and the four string contacts on the bridge, each driving its
/// node along x.
fn cannone_options() -> super::super::shell::assembly::table::OptionsFile {
    super::super::shell::assembly::table::OptionsFile::parse(include_str!("../../cannone-h8.assembly")).unwrap()
}

/// The il Cannone violin body (data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/shell-assembly-h8.tsv: 48
/// parts, 7042 nodes, 10192 triangles, 119 ties) written as a v4 assembled profile with `cannone_options`. The table's
/// frame has x across the body, y vertical and z along the strings.
fn cannone_profile() -> (Profile, super::super::shell::assembly::table::ConvertedTable) {
    cannone_profile_with(&cannone_options(), CANNONE_H8)
}

const CANNONE_H8: &str =
    include_str!("../../../../../../data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/shell-assembly-h8.tsv");

/// The same body meshed at 6 mm by the same generator (shell-assembly-h6.tsv: 64 parts, 11769 nodes, 17483 triangles,
/// 157 ties), the largest table on file.
const CANNONE_H6: &str =
    include_str!("../../../../../../data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/shell-assembly-h6.tsv");

/// The 6 mm table's assembly choices (`cannone-h6.assembly`): h8's, with the window to 5 kHz (the spruce pack's source
/// band) and the projection work its 11769 nodes need.
fn cannone_h6_options() -> super::super::shell::assembly::table::OptionsFile {
    super::super::shell::assembly::table::OptionsFile::parse(include_str!("../../cannone-h6.assembly")).unwrap()
}

fn cannone_profile_with(options: &super::super::shell::assembly::table::OptionsFile, table: &str)
    -> (Profile, super::super::shell::assembly::table::ConvertedTable)
{
    let converted = options.convert(table).unwrap();
    (Profile::parse(converted.profile.as_bytes()).unwrap(), converted)
}

/// A shell-assembly options file states each choice once: a missing, repeated, unknown or malformed record refuses, as do
/// an inverted window and another version.
#[test]
fn shell_assembly_options_refuse_missing_repeated_unknown_and_malformed_records() {
    use super::super::shell::assembly::table::OptionsFile;
    let text = include_str!("../../cannone-h8.assembly");
    assert!(OptionsFile::parse(text).is_ok());
    let refused = |t: &str| OptionsFile::parse(t).err().expect("must refuse").to_string();
    assert!(refused(&text.replace("window-hz 15 2000\n", "")).contains("need window-hz"));
    assert!(refused(&text.replace("bridge-direction 1 0 0\n", "")).contains("need bridge-direction"));
    assert!(refused(&format!("{text}reference-thickness-m 0.004\n")).contains("repeated"));
    assert!(refused(&format!("{text}pack maple violin-spruce-norway-viala2020-shell.fsmatpk\n")).contains("repeated"));
    assert!(refused(&format!("{text}bridge 0.01\n")).contains("needs 2 numbers"));
    assert!(refused(&format!("{text}bridge-component bridge top\n")).contains("malformed"));
    assert!(refused(&format!("{text}tie-budget 4\n")).contains("unknown"));
    assert!(refused(&text.replace("window-hz 15 2000", "window-hz 2000 15")).contains("0 < LO < HI"));
    assert!(refused(&text.replacen("options-v1", "options-v2", 1)).contains("not a frankensim-shell-assembly-options-v1"));
}

fn held_coating_pack() -> NormalizedPack {
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, MoistureBasis, QuantityKind, SemanticType, ValueForm}};
    // Manufactured passive films: the endpoints are not measured varnish data.
    let mut claims = ClaimSet::new();
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for (key, dims, low, high) in [
        ("density", Density::DIMS, 1000.0, 1200.0),
        ("varnish_e0", Pressure::DIMS, 2e8, 1.4e8),
        ("varnish_einf", Pressure::DIMS, 3e8, 2.1e8),
        ("varnish_alpha", Dims::NONE, 0.5, 0.5),
        ("varnish_beta", Dims::NONE, 0.5, 0.5),
        ("varnish_tau", Dims([0, 0, 1, 0, 0, 0]), 0.01, 0.02),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: PropertyValue::Curve { abscissa: "moisture_content".into(),
                abscissa_dims: Dims::NONE, knots: vec![(0.0, low), (0.2, high)], dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), 280.0, 320.0)
                .with_quantity("frequency", quantity(QuantityKind::Frequency(FrequencyConvention::Cyclic)), 0.01, 20_000.0)
                .with_quantity("moisture_content", quantity(QuantityKind::Moisture(MoistureBasis::DryMass)), 0.0, 0.2),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            provenance: Provenance { source: "manufactured held-state coating fixture".into(),
                license: "CC0-1.0".into(), artifact: None }, observations: vec![],
        }).unwrap();
    }
    NormalizedPack::new("synthetic-coating-environment", "synthetic-film-v1",
        fs_blake3::ContentHash([8; 32]), "CC0-1.0", claims, vec![], vec![]).unwrap()
}

#[test]
fn held_film_state_changes_original_laminate_and_coupled_pressure() {
    scoped(|cx, gate| {
        let source = held_coating_pack();
        let row = "coatings 1\ncoating body film.fsmatpk positive 0.00002 1000 0.3 0.001 0.35 0.9\n";
        let profile = Profile::parse(shell_text().replace("nodes 25\n",
            &format!("{row}nodes 25\n")).as_bytes()).unwrap();
        let prepare = |environment: Vec<(String, f64)>| -> Result<Ready> {
            let mut spec = shell_spec();
            if !environment.is_empty() { spec.material_states.insert("film.fsmatpk".into(), environment); }
            spec.compile_with_body(cx, gate, Some(profile.clone()), |name|
                Ok(if name == "film.fsmatpk" { source.clone() } else { pack(name, 1.0, None) }))
        };
        let environment = |value| vec![("moisture_content".into(), value)];
        let mut low = prepare(environment(0.0)).unwrap();
        let mut high = prepare(environment(0.2)).unwrap();
        let coating = |ready: &Ready| ready.performance.bank().source().bank()
            .material_shell().unwrap().specimen.coatings()[0].clone();
        let a = coating(&low);
        let b = coating(&high);
        assert_eq!(a.faces_m(), b.faces_m());
        assert_eq!(a.input().material.source().card_identity(), b.input().material.source().card_identity());
        assert!((b.mass_kg()/a.mass_kg() - 1.2).abs() < 1e-12);
        assert_ne!(a.reference_storage_pa(), b.reference_storage_pa());
        assert_eq!(b.input().material.law().relaxation_time_s(), 0.02);
        assert_eq!(b.input().material.source().state_coordinate("moisture_content"), Some(0.2));
        assert_eq!(b.input().material.source().state_coordinate("conditioning_relative_humidity"), Some(0.35));
        let modes = |r: &Ready| r.performance.bank().system().components()[1].modes().to_vec();
        assert_ne!(modes(&low), modes(&high));
        assert_ne!(low.performance.bank().source().bank().material_shell_image().unwrap().damping().matrix_per_s(),
            high.performance.bank().source().bank().material_shell_image().unwrap().damping().matrix_per_s());
        for name in ["wire.fsmatpk", "sleeve.fsmatpk", "body.fsmatpk", "film.fsmatpk"] {
            assert_eq!(low.catalog[name].pack.content_hash(), high.catalog[name].pack.content_hash());
        }
        let low_pressure = pressure(&mut low, 7);
        let high_pressure = pressure(&mut high, 31);
        assert_ne!(low_pressure, high_pressure);
        assert_eq!(high_pressure, pressure(&mut prepare(environment(0.2)).unwrap(), 1));
        for ready in [&low, &high] {
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(prepare(vec![]).is_err());
        assert!(prepare(environment(0.3)).is_err());
        assert!(prepare(vec![("conditioning_relative_humidity".into(), 0.2)]).is_err(),
            "held environment must not overwrite the coating's authored conditioning");
    });
}

/// A finish is an explicit layer on the selected side of an existing region.
/// Malformed geometry, tensor assumptions and bounds refuse before pack loading.
#[test]
fn shell_coatings_require_explicit_regions_sources_and_physical_assumptions() {
    let row = "coatings 1\ncoating body varnish.fsmatpk positive 0.00002 1000 0.3 0.001 0.35 0.25\n";
    let text = shell_text().replace("nodes 25\n", &format!("{row}nodes 25\n"));
    assert!(Profile::parse(text.as_bytes()).is_ok());
    for (old, new) in [
        ("coating body", "coating missing"),
        ("varnish.fsmatpk", "../varnish.fsmatpk"),
        ("varnish.fsmatpk", "@body"),
        ("positive 0.00002", "outside 0.00002"),
        ("positive 0.00002", "positive 0"),
        ("1000 0.3", "0 0.3"),
        ("1000 0.3", "1000 0.5"),
        ("0.001 0.35", "0 0.35"),
        ("0.35 0.25", "1.1 0.25"),
        ("0.35 0.25", "0.35 1"),
        ("0.35 0.25", "0.35 NaN"),
        ("coatings 1", "coatings 65"),
    ] {
        assert!(Profile::parse(text.replace(old, new).as_bytes()).is_err(), "accepted {new}");
    }
    assert!(Profile::parse(shell_text().replace("nodes 25\n", "coatings 0\nnodes 25\n").as_bytes()).is_ok());
}

#[test]
fn shell_coating_memory_requires_the_explicit_continuation_and_complete_budget() {
    let row = "coatings 1\ncoating body varnish.fsmatpk positive 0.00002 1000 0.3 0.001 0.35 0.25\n";
    let film = shell_text().replace("nodes 25\n", &format!("{row}nodes 25\n"));
    let text = film.replace("regions 1\n",
        "coating-memory prony-continuation 16 0.001 0.001 64\nregions 1\n");
    let profile = Profile::parse(text.as_bytes()).unwrap();
    assert_eq!(profile.material_memory_states(), 64);
    for (old, new) in [
        ("prony-continuation", "measured-broadband"),
        ("16 0.001", "0 0.001"),
        ("16 0.001", "33 0.001"),
        ("0.001 0.001 64", "0 0.001 64"),
        ("0.001 0.001 64", "0.001 1 64"),
        ("0.001 0.001 64", "0.001 NaN 64"),
        ("0.001 0.001 64", "0.001 0.001 0"),
        ("0.001 0.001 64", "0.001 0.001 3073"),
    ] {
        assert!(Profile::parse(text.replace(old, new).as_bytes()).is_err(), "accepted {new}");
    }
    assert!(Profile::parse(text.replace(row, "").as_bytes()).is_err());
    assert!(Profile::parse(text.replace("regions 1\n",
        "coating-memory prony-continuation 16 0.001 0.001 64\nregions 1\n").as_bytes()).is_err());
}

#[test]
fn cannone_memory_selection_keeps_the_same_film_and_reserves_one_assembly_budget() {
    use super::super::shell::assembly::table::OptionsFile;
    let memory = OptionsFile::parse(include_str!("../../cannone-h8-varnish-memory.assembly"))
        .unwrap().convert(CANNONE_H8).unwrap();
    let frozen = OptionsFile::parse(include_str!("../../cannone-h8-varnished.assembly"))
        .unwrap().convert(CANNONE_H8).unwrap();
    let profile = Profile::parse(memory.profile.as_bytes()).unwrap();
    assert_eq!(profile.material_memory_states(), 512);
    let stripped: String = memory.profile.lines().filter(|line| !line.starts_with("coating-memory "))
        .map(|line| format!("{line}\n")).collect();
    assert_eq!(stripped, frozen.profile);
    let text = include_str!("../../cannone-dominant-g3.material-strings");
    profile.preflight(&input::Spec::parse(text.as_bytes()).unwrap()).unwrap();
    let unprepared: String = text.lines().filter(|line| !line.starts_with("prepared-linear "))
        .map(|line| format!("{line}\n")).collect();
    assert!(profile.preflight(&input::Spec::parse(unprepared.as_bytes()).unwrap()).is_err());
}

/// The actual violin table carries 45 exterior components and 60 material
/// regions. An exterior selection covers them, without coating internal parts.
#[test]
fn cannone_varnish_selection_covers_only_declared_exterior_regions() {
    use super::super::shell::assembly::table::OptionsFile;
    let text = include_str!("../../cannone-h8-varnished.assembly");
    let options = OptionsFile::parse(text).unwrap();
    let coated = options.convert(CANNONE_H8).unwrap();
    let bare = cannone_options().convert(CANNONE_H8).unwrap();
    assert!(Profile::parse(coated.profile.as_bytes()).is_ok());
    assert_eq!(coated.length_m, bare.length_m);
    assert_eq!(coated.width_m, bare.width_m);
    assert_eq!(coated.bridge_snaps_m, bare.bridge_snaps_m);
    assert_eq!(coated.cavity, bare.cavity);
    assert_eq!(coated.profile.lines().filter(|line| line.starts_with("coatings ")).count(), 45);
    assert_eq!(coated.profile.lines().filter(|line| line.starts_with("coating ")).count(), 60);
    for component in ["bass-bar", "soundpost", "bridge"] {
        let section = coated.profile.split(&format!("component {component}\n")).nth(1).unwrap()
            .split("end-component").next().unwrap();
        assert!(!section.contains("coating "), "coated internal component {component}");
    }
    let stripped: String = coated.profile.lines().filter(|line| !line.starts_with("coating ")
        && !line.starts_with("coatings ")).map(|line| format!("{line}\n")).collect();
    assert_eq!(stripped, bare.profile, "adding a film changed another assembly choice");
    let unknown = text.replace("coating @exterior", "coating absent-part");
    assert!(OptionsFile::parse(&unknown).unwrap().convert(CANNONE_H8).is_err());
    let doubled = text.replace("coating @exterior", "coating rib");
    let options = OptionsFile::parse(&doubled).unwrap().convert(CANNONE_H8).unwrap();
    assert_eq!(options.profile.lines().filter(|line| line.starts_with("coating ")).count(), 16);
}

/// `cannone-dominant-g3.material-strings` is the il Cannone G3 note for the command line: the four Dominant string packs
/// on the rung's 0.383 m span and 0.328 m speaking length, the G bowed at 1 N and 0.2 m/s one tenth of its speaking length
/// from the bridge, 2 s at 144 kHz. Its body and cavity records are the ones this body's conversion gives, and it passes
/// preflight on that profile.
#[test]
fn the_cannone_dominant_g3_performance_matches_its_body_and_passes_preflight() {
    let text = include_str!("../../cannone-dominant-g3.material-strings");
    let (profile, converted) = cannone_profile();
    let [volume, neck_radius, neck_length] = converted.cavity.unwrap();
    let body = format!("body violin-spruce-norway-viala2020-shell.fsmatpk {} {} {} {CANNONE_BODY_MODES} 1\n",
        converted.length_m, converted.width_m, cannone_options().reference_thickness_m());
    assert!(text.contains(&body), "{body}");
    let cavity = format!("cavity physical-neck {volume} {neck_radius} {neck_length} flanged flanged 1\n");
    assert!(text.contains(&cavity), "{cavity}");
    for (k, string) in ["g", "d", "a", "e"].iter().enumerate() {
        let record = text.lines().filter(|l| l.starts_with("string ")).nth(k).unwrap();
        assert!(record.starts_with(&format!("string violin-string-dominant-{string}-pickering1985-section.fsmatpk ")));
    }
    profile.preflight(&input::Spec::parse(text.as_bytes()).unwrap()).unwrap();
}

/// The crown joints of `cannone-dominant-g3-crown-{compliant,stiff}.material-strings`, written by
/// `gen_crown_interface_packs.py`. Each Dominant string's lateral stick on the il Cannone bridge crown comes from Hertz
/// crossed cylinders and Mindlin's no-slip tangential stiffness k_t = 8 G* a (Johnson 1985, eqs. 4.53 and 7.42).
///
/// Both ends of the envelope of its STATED UNKNOWNs are covered:
/// - the maple orthotropy: G_RT with E_L, or G_LR with E_R;
/// - the crown edge radius: 0.1 mm, or half the CT bridge-top thickness;
/// - the unnamed wrap metal.
///
/// The test checks five things:
/// - each note is the committed note with only every string's scalar joint replaced, and the stiff end's audio clock
///   raised to 384 kHz over the same duration;
/// - re-derived here from the committed seeds, the CT tables and the note's own tensions, every resolved joint agrees to
///   2e-8 in stiffness and damping;
/// - every joint states the maple pack's constant-loss anchor, so its loss factor is share x eta at every frequency;
/// - every string's compliant end is the softer and the less lossy;
/// - the joints claim 0 to 2 pi 200 kHz (MagentaBadger 7268 (a)), and a declared state above that refuses.
#[test]
fn the_cannone_crown_joints_resolve_to_their_contact_mechanics_at_both_envelope_ends() {
    let scalar = |text: &str, name: &str| -> f64 {
        let f = text.lines().map(|l| l.split('\t').collect::<Vec<_>>())
            .find(|f| f[0] == "scalar" && f[3] == name).unwrap_or_else(|| panic!("no scalar {name}"));
        f[4].parse::<f64>().unwrap() * match f[5] { "MPa" => 1e6, "GPa" => 1e9, _ => 1.0 }
    };
    let maple = include_str!(
        "../../../../../../data/matdb/seed-v1/violin-maple-sycamore-sonderegger2013-shell/properties.tsv");
    let aluminium = include_str!("../../../../../../data/matdb/seed-v1/aluminum-2024-t3-sheet-mil-hdbk-5j/properties.tsv");
    let silver = include_str!("../../../../../../data/matdb/seed-v1/silver-pure-smith-fickett-nist1995/properties.tsv");
    let silver_shear = include_str!("../../../../../../data/matdb/seed-v1/silver-pure-wrap-nist1995/properties.tsv");
    let wire = include_str!("../../../../../../data/matdb/seed-v1/music-wire-spring-moduli-fuchs-1968/properties.tsv");
    let contacts = include_str!(
        "../../../../../../data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/string-contacts.tsv");
    let dimensions = include_str!(
        "../../../../../../data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/dimensions.tsv");
    let base = include_str!("../../cannone-dominant-g3.material-strings");
    let notes = [("compliant", include_str!("../../cannone-dominant-g3-crown-compliant.material-strings")),
        ("stiff", include_str!("../../cannone-dominant-g3-crown-stiff.material-strings"))];
    let packs: [(&str, &[u8]); 4] = [
        ("violin-cannone-crown-dominant-e-interface.fsintpk",
            include_bytes!("../../packs/violin-cannone-crown-dominant-e-interface.fsintpk")),
        ("violin-cannone-crown-dominant-a-interface.fsintpk",
            include_bytes!("../../packs/violin-cannone-crown-dominant-a-interface.fsintpk")),
        ("violin-cannone-crown-dominant-d-interface.fsintpk",
            include_bytes!("../../packs/violin-cannone-crown-dominant-d-interface.fsintpk")),
        ("violin-cannone-crown-dominant-g-interface.fsintpk",
            include_bytes!("../../packs/violin-cannone-crown-dominant-g-interface.fsintpk")),
    ];
    let loader = |name: &str| -> Result<NormalizedInterfacePack> {
        let (_, bytes) = packs.iter().find(|(n, _)| *n == name).ok_or("unknown crown pack")?;
        Ok(NormalizedInterfacePack::from_bytes(bytes)?)
    };

    // The contact, re-derived from the same committed sources (gen_crown_interface_packs.py's derivation).
    let nu_maple = scalar(maple, "poisson_ratio_12");
    let crown = [(scalar(maple, "young_modulus_1"), scalar(maple, "shear_modulus_23")),
        (scalar(maple, "young_modulus_2"), scalar(maple, "shear_modulus_12"))];
    let tau = [("plane_stress_viscosity_66", "shear_modulus_12"), ("transverse_shear_viscosity_13", "shear_modulus_13"),
        ("transverse_shear_viscosity_23", "shear_modulus_23")]
        .map(|(v, g)| scalar(maple, v) / scalar(maple, g)).into_iter().fold(0.0_f64, f64::max);
    let crown_thickness = dimensions.lines().map(|l| l.split('\t').collect::<Vec<_>>())
        .find(|f| f.len() > 4 && f[0] == "bridge" && f[2] == "thickness-z-min-above-y-minus-52" && f[4] == "m")
        .map(|f| f[3].parse::<f64>().unwrap()).unwrap();
    let edge_radius = [1.0e-4, 0.5 * crown_thickness];
    let al = (scalar(aluminium, "young_modulus"), scalar(aluminium, "poisson_ratio"), scalar(aluminium, "shear_modulus"));
    let ag = (scalar(silver, "young_modulus"), scalar(silver, "poisson_ratio"), scalar(silver_shear, "shear_modulus"));
    let (e_w, g_w) = (scalar(wire, "young_modulus"), scalar(wire, "shear_modulus"));
    let steel = (e_w, e_w / (2.0 * g_w) - 1.0, g_w);
    let joint = |string: &str, end: usize, radius: f64, tension: f64| -> (f64, f64) {
        let f: Vec<&str> = contacts.lines().map(|l| l.split('\t').collect::<Vec<_>>())
            .find(|f| f.len() == 7 && f[0] == string && f[1] == "bridge").unwrap();
        let theta = (-f[5].parse::<f64>().unwrap()).atan() + f[6].parse::<f64>().unwrap().atan();
        let load = 2.0 * tension * (theta / 2.0).sin();
        let (e_s, nu_s, g_s) = [al, if string == "e" { steel } else { ag }][end];
        let (e_b, g_b) = crown[end];
        let e_star = 1.0 / ((1.0 - nu_s * nu_s) / e_s + (1.0 - nu_maple * nu_maple) / e_b);
        let (maple_shear, metal_shear) = ((2.0 - nu_maple) / g_b, (2.0 - nu_s) / g_s);
        let a = (3.0 * load * (radius * edge_radius[end]).sqrt() / (4.0 * e_star)).cbrt();
        let k = 8.0 * a / (maple_shear + metal_shear);
        (k, maple_shear / (maple_shear + metal_shear) * tau * k)
    };

    let mut resolved = Vec::new();
    for (end, (label, text)) in notes.iter().enumerate() {
        let (changed, original) = (text.lines().collect::<Vec<_>>(), base.lines().collect::<Vec<_>>());
        assert_eq!(changed.len(), original.len());
        for (c, o) in changed.iter().zip(&original) {
            if end == 1 && o.starts_with("audio ") {
                // The stiff end renders at 384 kHz over the same duration (MagentaBadger 7268 (b)).
                let (c, o): (Vec<&str>, Vec<&str>) = (c.split(' ').collect(), o.split(' ').collect());
                let seconds = |f: &[&str]| f[2].parse::<f64>().unwrap() / f[1].parse::<f64>().unwrap();
                assert_eq!((c[1], &c[3..], seconds(&c)), ("384000", &o[3..], seconds(&o)));
                continue;
            }
            if !o.starts_with("string ") {
                assert_eq!(c, o);
                continue;
            }
            let (c, o): (Vec<&str>, Vec<&str>) = (c.split(' ').collect(), o.split(' ').collect());
            assert_eq!((&c[..9], &o[9..]), (&o[..9], &["200000", "0.1"][..]), "{label}");
            assert_eq!((c.len(), c[9]), (14, "interface"));
            assert_eq!(c[12], format!("crown_lateral_stiffness_per_area_{label}"));
            assert_eq!(c[13], format!("crown_lateral_viscosity_per_area_{label}"));
        }
        let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
        let gate = CancelGate::new_clock_free();
        let claimed = core::f64::consts::TAU * 200e3;
        assert!(spec.bridge_interfaces.load(spec.temperature, 1.000_001 * claimed, &gate, loader).is_err(),
            "a solid state above the joints' 2 pi 200 kHz must refuse");
        spec.bridge_interfaces.load(spec.temperature, spec.omega, &gate, loader).unwrap();
        let bindings = spec.bridge_interfaces.bindings();
        assert_eq!(bindings.len(), 4);
        for binding in bindings {
            let row: Vec<&str> = changed.iter().filter(|l| l.starts_with("string ")).nth(binding.string).unwrap()
                .split(' ').collect();
            let string = &row[1]["violin-string-dominant-".len()..][..1];
            assert_eq!(row[10], format!("violin-cannone-crown-dominant-{string}-interface.fsintpk"));
            let (k, d) = joint(string, end, row[3].parse().unwrap(), row[5].parse().unwrap());
            // The source per area and the row's area are each printed to 9 significant digits.
            let law = &binding.interface;
            assert_eq!(law.constant_loss_reference_omega_rad_s(), Some(scalar(maple, "constant_loss_reference_omega")));
            let (lo, hi) = law.omega_band_rad_s();
            assert!(lo == 0.0 && (hi / claimed - 1.0).abs() < 1e-8, "{label} {string}: band {lo}..{hi}");
            assert!((law.stiffness_n_m() / k - 1.0).abs() < 2e-8, "{label} {string}: {} vs {k}", law.stiffness_n_m());
            assert!((law.damping_n_s_m() / d - 1.0).abs() < 2e-8, "{label} {string}: {} vs {d}", law.damping_n_s_m());
            resolved.push((end, binding.string, law.stiffness_n_m(), law.damping_n_s_m()));
        }
    }
    for &(_, string, k, d) in resolved.iter().filter(|r| r.0 == 0) {
        let &(_, _, k_stiff, d_stiff) = resolved.iter().find(|r| r.0 == 1 && r.1 == string).unwrap();
        assert!(k < k_stiff && d < d_stiff, "string {string}: compliant {k} {d}, stiff {k_stiff} {d_stiff}");
    }
}

/// Both il Cannone tables are admitted, and their four bridges land on bridge nodes within 0.6 mm of the string
/// contacts and push along x. The 6 mm body's 64 components and 157 ties are exactly the assembly caps: one more
/// component or one more tie refuses. Its rib set (2904 nodes, 5082 triangles) and its 11769 nodes and 17483 triangles
/// in all fit the node and triangle caps. With the violin-sized strings (synthetic steel, `violin_performance`), each
/// table's own cavity row as a physical neck and 128 body modes, each passes preflight: every check that runs before
/// pack I/O and the reduction.
#[test]
fn the_cannone_violin_body_at_6_mm_fills_the_assembly_caps_and_both_tables_pass_preflight() {
    let records = |text: &str, name: &str| text.lines().filter(|l| l.split(' ').next() == Some(name)).count();
    let tables = [(cannone_options(), CANNONE_H8, (48, 119, 4)), (cannone_h6_options(), CANNONE_H6, (64, 157, 4))];
    for (options, table, counts) in tables {
        let (profile, converted) = cannone_profile_with(&options, table);
        let text = &converted.profile;
        assert_eq!((records(text, "component"), records(text, "edge-tie"), records(text, "bridge")), counts);
        assert!(text.lines().filter(|l| l.starts_with("bridge bridge ")).all(|l| l.ends_with(" 1 0 0")));
        assert!(converted.bridge_snaps_m.iter().all(|&d| d < 0.000_6), "bridge snaps {:?} m", converted.bridge_snaps_m);
        profile.preflight(&cannone_inputs_with(&options, table, 96, 128).1).unwrap();
    }
    let text = &cannone_profile_with(&cannone_h6_options(), CANNONE_H6).1.profile;
    let refusal = |text: &str| Profile::parse(text.as_bytes()).err().expect("one past the cap must refuse").to_string();
    // One more component: the bridge again, under a new name.
    let (start, end) = (text.find("component bridge\n").unwrap(), text.find("face-ties").unwrap());
    let part_end = start + text[start..end].find("end-component\n").unwrap() + "end-component\n".len();
    let copy = text[start..part_end].replacen("component bridge\n", "component bridge-copy\n", 1);
    let one_more_component = format!("{}{copy}{}", &text[..part_end], &text[part_end..])
        .replacen("components 64\n", "components 65\n", 1);
    assert!(refusal(&one_more_component).contains("count exceeds its explicit range"));
    // One more tie: the first edge tie again.
    let lines: Vec<&str> = text.lines().collect();
    let first = lines.iter().position(|l| l.starts_with("edge-tie ")).unwrap();
    let f: Vec<&str> = lines[first].split(' ').collect();
    let block = 1 + f[5].parse::<usize>().unwrap() + f[6].parse::<usize>().unwrap();
    let one_more_tie = [&lines[..first + block], &lines[first..]].concat().join("\n")
        .replacen("edge-ties 157\n", "edge-ties 158\n", 1);
    assert!(refusal(&one_more_tie).contains("count exceeds its explicit range"));
}

/// The violin-sized G3 rung's strings, bow and schedule (`violin_performance`: synthetic steel, NOT violin strings) on the
/// il Cannone body (`cannone_profile`) with `modes` retained body modes and the table's own cavity row (2.061 L through
/// one equivalent neck of the two f-holes, both ends flanged). The body, the cavity and the strings' directions
/// (`CANNONE_DIRECTIONS`) differ from that rung.
fn cannone_inputs(samples: u64, modes: usize) -> (Profile, input::Spec) {
    cannone_inputs_with(&cannone_options(), CANNONE_H8, samples, modes)
}

/// `cannone_inputs` on one il Cannone table converted with its options.
fn cannone_inputs_with(options: &super::super::shell::assembly::table::OptionsFile, table: &str, samples: u64,
    modes: usize) -> (Profile, input::Spec)
{
    let (profile, converted) = cannone_profile_with(options, table);
    let [volume, neck_radius, neck_length] = converted.cavity.unwrap();
    let performance = violin_performance(samples, &format!("body violin-spruce-norway-viala2020-shell.fsmatpk {} {} {} \
        {modes} 1", converted.length_m, converted.width_m, options.reference_thickness_m()),
        &format!("cavity physical-neck {volume} {neck_radius} {neck_length} flanged flanged 1"), CANNONE_DIRECTIONS);
    (profile, input::Spec::parse(performance.as_bytes()).unwrap())
}

/// Every assembled il Cannone mode between 15 Hz and 2 kHz: the window holds 30 (a request for 128 refuses, naming 30).
const CANNONE_BODY_MODES: usize = 30;

/// The two violin shell packs, built by `xtask matdb-pack` from data/matdb/seed-v1, and the synthetic steel elsewhere.
fn cannone_pack(name: &str) -> NormalizedPack {
    let bytes: &[u8] = match name {
        "violin-spruce-norway-viala2020-shell.fsmatpk" =>
            include_bytes!("../../packs/violin-spruce-norway-viala2020-shell.fsmatpk"),
        "violin-maple-sycamore-sonderegger2013-shell.fsmatpk" =>
            include_bytes!("../../packs/violin-maple-sycamore-sonderegger2013-shell.fsmatpk"),
        _ => return pack(name, 1.0, None),
    };
    NormalizedPack::from_bytes(bytes).unwrap()
}

/// Each retained body mode as rendered: frequency [Hz], and its Q from the diagonal of the modal viscosity matrix. Both
/// packs state their Kelvin-Voigt reference (`constant_loss_reference_omega`, 2 pi 200 rad/s), so every mode's damping
/// carries their loss factors at its own frequency: Q is about one over the loss factor of the mode's material mix.
fn cannone_body_modes(r: &Ready) -> String {
    let image = r.performance.bank().source().bank().material_shell_assembly_image().unwrap();
    let damping = image.damping().matrix_per_s();
    let n = image.eigenvalue_intervals().len();
    image.eigenvalue_intervals().iter().enumerate().map(|(i, &(lo, hi))| {
        let omega = (0.5 * (lo + hi)).sqrt();
        let (f, q) = (omega / (2.0 * std::f64::consts::PI), omega / damping[i * n + i]);
        format!("{f:.1} Hz Q {q:.1}")
    }).collect::<Vec<_>>().join(", ")
}

/// The il Cannone G3 rung (`cannone_inputs`) for a quarter second: the wall time per audio second with what it was
/// measured on, the time of construction (assembled shell modes, projection, prepared network) and render, every
/// retained body mode, and the pressure bits.
#[test]
#[ignore = "il Cannone violin G3 cost rung: run with --ignored --nocapture"]
fn the_cannone_violin_g3_rung_reports_its_body_modes_and_wall_time_per_audio_second() {
    let (profile, spec) = cannone_inputs(12_000, CANNONE_BODY_MODES);
    scoped(|cx, gate| {
        let constructing = std::time::Instant::now();
        let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(cannone_pack(name))).unwrap();
        let built = constructing.elapsed();
        eprintln!("il Cannone G3 rung: {}", cost_context(&r.spec));
        eprintln!("il Cannone G3 rung: body modes {}", cannone_body_modes(&r));
        let source = r.performance.bank().source().bank();
        assert_eq!(source.material_shell_assembly_image().unwrap().eigenvalue_intervals().len(), CANNONE_BODY_MODES);
        assert_eq!(source.strings().len(), 4);
        assert!(source.material_cavity().is_some());
        let rendering = std::time::Instant::now();
        let samples = pressure(&mut r, 128);
        let rendered = rendering.elapsed();
        assert_eq!(samples.len(), 12_000);
        assert!(samples.iter().all(|p| p.is_finite()));
        assert!(samples.iter().any(|p| p.abs() > 1e-14), "the bowed body radiated nothing");
        let per_second = rendered.as_secs_f64() / (samples.len() as f64 / 48_000.0);
        eprintln!("il Cannone G3 rung: constructed in {built:.1?}, {} samples rendered in {rendered:.1?} \
            ({per_second:.1} s render per audio second; a 2 s note with its construction takes {:.1} s per audio second)",
            samples.len(), (built.as_secs_f64() + 2.0 * per_second) / 2.0);
        eprintln!("il Cannone G3 rung: pressure bits fnv1a {:#018x}, peak {:.3e} Pa", pressure_bits(&samples),
            samples.iter().fold(0.0_f64, |m, p| m.max(p.abs())));
    });
}

/// The il Cannone G3 rung bowed for 2 s, for listening and for comparison with the violin-sized rung's 2 s render on the
/// clamped plate (the same strings, bow and schedule). FRANKENSIM_CANNONE_RUNG_OUT names an existing directory, which
/// receives `cannone-g3-rung-2s.f64le` (every pressure sample [Pa], exact) and `cannone-g3-rung-2s.wav` (PCM16 whose full
/// scale is the smallest power of two in pascals above the peak; no normalization or mastering, no clip).
#[test]
#[ignore = "two-second il Cannone render for listening: set FRANKENSIM_CANNONE_RUNG_OUT, run with --ignored --nocapture"]
fn the_cannone_violin_g3_rung_renders_two_seconds_for_listening() {
    let out = std::path::PathBuf::from(std::env::var_os("FRANKENSIM_CANNONE_RUNG_OUT")
        .expect("FRANKENSIM_CANNONE_RUNG_OUT names the output directory"));
    let (profile, spec) = cannone_inputs(96_000, CANNONE_BODY_MODES);
    scoped(|cx, gate| {
        let constructing = std::time::Instant::now();
        let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(cannone_pack(name))).unwrap();
        let built = constructing.elapsed();
        eprintln!("il Cannone G3 2 s: body modes {}", cannone_body_modes(&r));
        let rendering = std::time::Instant::now();
        let samples = pressure(&mut r, 128);
        let rendered = rendering.elapsed();
        assert_eq!(samples.len(), 96_000);
        assert!(samples.iter().all(|p| p.is_finite()));
        let peak = samples.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
        assert!(peak > 0.0, "the bowed body radiated nothing");
        let full_scale = 2f64.powi(peak.log2().floor() as i32 + 1);
        let (wav, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&samples, 48_000, full_scale).unwrap();
        assert_eq!(clips, 0);
        std::fs::write(out.join("cannone-g3-rung-2s.wav"), wav).unwrap();
        std::fs::write(out.join("cannone-g3-rung-2s.f64le"),
            samples.iter().flat_map(|p| p.to_le_bytes()).collect::<Vec<u8>>()).unwrap();
        eprintln!("il Cannone G3 2 s: constructed in {built:.1?}, {} samples rendered in {rendered:.1?} \
            ({:.1} s per audio second with construction)", samples.len(),
            (built + rendered).as_secs_f64() / (samples.len() as f64 / 48_000.0));
        eprintln!("il Cannone G3 2 s: {}", cost_context(&r.spec));
        eprintln!("il Cannone G3 2 s: peak {peak:.4e} Pa, WAV full scale {full_scale:e} Pa, pressure bits fnv1a {:#018x} \
            (first 12000 samples {:#018x})", pressure_bits(&samples), pressure_bits(&samples[..12_000]));
    });
}

/// The violin-sized G3 rung (above) for a quarter second: the wall time per audio second with what it was measured on, the
/// time of each stage (construction: shell modes, projection and the prepared network; render), and the pressure bits.
#[test]
#[ignore = "violin-sized cost rung of bead .358: run with --ignored --nocapture"]
fn the_violin_sized_g3_rung_reports_its_wall_time_per_audio_second() {
    let (profile, spec) = violin_inputs(12_000);
    scoped(|cx, gate| {
        let constructing = std::time::Instant::now();
        let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(if name == "body.fsmatpk" {
            pack_with_shear(name, 1.0, 1.0, None) } else { pack(name, 1.0, None) })).unwrap();
        let built = constructing.elapsed();
        eprintln!("violin G3 rung: {}", cost_context(&r.spec));
        let source = r.performance.bank().source().bank();
        assert_eq!(source.material_shell_image().unwrap().physical_modes().len(), 128);
        assert_eq!(source.strings().len(), 4);
        assert!(source.material_cavity().is_some());
        let rendering = std::time::Instant::now();
        let samples = pressure(&mut r, 128);
        let rendered = rendering.elapsed();
        assert_eq!(samples.len(), 12_000);
        assert!(samples.iter().all(|p| p.is_finite()));
        assert!(samples.iter().any(|p| p.abs() > 1e-14), "the bowed body radiated nothing");
        let audio_s = samples.len() as f64 / 48_000.0;
        let per_second = rendered.as_secs_f64() / audio_s;
        eprintln!("violin G3 rung: constructed in {built:.1?}, {} samples rendered in {rendered:.1?} \
            ({per_second:.1} s render per audio second; a 2 s note with its construction takes {:.1} s per audio second)",
            samples.len(), (built.as_secs_f64() + 2.0 * per_second) / 2.0);
        eprintln!("violin G3 rung: {}", cost_context(&r.spec));
        eprintln!("violin G3 rung: pressure bits fnv1a {:#018x}, peak {:.3e} Pa", pressure_bits(&samples),
            samples.iter().fold(0.0_f64, |m, p| m.max(p.abs())));
    });
}

/// The violin-sized G3 rung bowed for 2 s, for listening, P2 grading and P8's blind A/B (bead .358). It is synthetic
/// throughout; see `violin_inputs`. FRANKENSIM_VIOLIN_RUNG_OUT names an existing directory, which receives:
/// - `violin-g3-rung-2s.f64le`: every pressure sample [Pa], exact;
/// - `violin-g3-rung-2s.wav`: PCM16 whose full scale is the smallest power of two in pascals above the peak. The full scale
///   is printed. There is no normalization or mastering, and the WAV has no clip.
#[test]
#[ignore = "two-second violin-sized render for listening: set FRANKENSIM_VIOLIN_RUNG_OUT, run with --ignored --nocapture"]
fn the_violin_sized_g3_rung_renders_two_seconds_for_listening() {
    let out = std::path::PathBuf::from(std::env::var_os("FRANKENSIM_VIOLIN_RUNG_OUT")
        .expect("FRANKENSIM_VIOLIN_RUNG_OUT names the output directory"));
    let (profile, spec) = violin_inputs(96_000);
    scoped(|cx, gate| {
        let constructing = std::time::Instant::now();
        let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(if name == "body.fsmatpk" {
            pack_with_shear(name, 1.0, 1.0, None) } else { pack(name, 1.0, None) })).unwrap();
        let built = constructing.elapsed();
        let rendering = std::time::Instant::now();
        let samples = pressure(&mut r, 128);
        let rendered = rendering.elapsed();
        assert_eq!(samples.len(), 96_000);
        assert!(samples.iter().all(|p| p.is_finite()));
        let peak = samples.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
        assert!(peak > 0.0, "the bowed body radiated nothing");
        let full_scale = 2f64.powi(peak.log2().floor() as i32 + 1);
        let (wav, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&samples, 48_000, full_scale).unwrap();
        assert_eq!(clips, 0);
        std::fs::write(out.join("violin-g3-rung-2s.wav"), wav).unwrap();
        std::fs::write(out.join("violin-g3-rung-2s.f64le"),
            samples.iter().flat_map(|p| p.to_le_bytes()).collect::<Vec<u8>>()).unwrap();
        eprintln!("violin G3 2 s: constructed in {built:.1?}, {} samples rendered in {rendered:.1?} \
            ({:.1} s per audio second with construction)", samples.len(),
            (built + rendered).as_secs_f64() / (samples.len() as f64 / 48_000.0));
        eprintln!("violin G3 2 s: {}", cost_context(&r.spec));
        // A longer horizon must not change earlier samples: the first quarter second is compared against the bits the cost
        // rung prints for its own quarter-second render.
        eprintln!("violin G3 2 s: peak {peak:.4e} Pa, WAV full scale {full_scale:e} Pa, pressure bits fnv1a {:#018x} \
            (first 12000 samples {:#018x})", pressure_bits(&samples), pressure_bits(&samples[..12_000]));
    });
}

/// The six strings of the guitar E4 rung, E2 A2 D3 G3 B3 E4: pitch [Hz] and radius [m] of plain strings in the synthetic
/// steel of these tests, NOT guitar strings (the nylon and wound-bass packs in data/matdb are literature facts, not the
/// string-section properties this route reads). Each speaking length's tension is derived: 70 N each, within 0.4 N.
const GUITAR_STRINGS: [(f64, f64); 6] =
    [(82.41, 0.000_5), (110.0, 0.000_375), (146.83, 0.000_28), (196.0, 0.000_21), (246.94, 0.000_166), (329.63, 0.000_125)];

/// Every guitar string keeps every bending and twist mode of its full span below this frequency [Hz]. It is half the
/// violin rung's, set when the prepared ceiling was 512 coordinates: six strings at 10 kHz with the 128 body modes
/// exceeded it.
const GUITAR_STRING_BAND_HZ: f64 = 5_000.0;

/// The Esteso 1964 guitar body (data/instruments/geometry/guitar-esteso-1964-zimmer2024/shell-assembly-h12.tsv, 19 parts,
/// 5805 nodes) written as a v4 assembled profile, with its three FPL-GTR-282 shell packs. The spruce top is the performance
/// body pack. The window holds the body modes from 15 Hz (above the free body's rigid motion) to 2 kHz, the packs' declared
/// angular-frequency validity.
fn esteso_profile(strings: &[(f64, f64)]) -> (Profile, super::super::shell::assembly::table::ConvertedTable) {
    let two_pi = 2.0 * std::f64::consts::PI;
    let converted = super::super::shell::assembly::table::convert(
        include_str!("../../../../../../data/instruments/geometry/guitar-esteso-1964-zimmer2024/shell-assembly-h12.tsv"),
        &super::super::shell::assembly::table::TableOptions {
            window: ((two_pi * 15.0).powi(2), (two_pi * 2000.0).powi(2)),
            projection_work: MAX_SHELL_PROJECTION_WORK,
            packs: &[("spruce", "@body"), ("rosewood", "guitar-rosewood-brazilian-fpl-gtr282-bremaud2011-shell.fsmatpk"),
                ("cedar", "guitar-cedar-western-red-fpl-gtr282-shell.fsmatpk")],
            reference_thickness_m: ESTESO_REFERENCE_THICKNESS_M,
            bridge_component: "top",
            bridges: strings,
            bridge_direction: [0.0, 0.0, 1.0],
        }).unwrap();
    (Profile::parse(converted.profile.as_bytes()).unwrap(), converted)
}

/// The performance body record's thickness [m]; every part's own thickness is written as a ratio to it.
const ESTESO_REFERENCE_THICKNESS_M: f64 = 0.0025;

/// The three guitar shell packs, built by `xtask matdb-pack` from data/matdb/seed-v1, and the synthetic steel elsewhere.
fn esteso_pack(name: &str) -> NormalizedPack {
    let bytes: &[u8] = match name {
        "guitar-spruce-sitka-fpl-gtr282-shell.fsmatpk" =>
            include_bytes!("../../packs/guitar-spruce-sitka-fpl-gtr282-shell.fsmatpk"),
        "guitar-rosewood-brazilian-fpl-gtr282-bremaud2011-shell.fsmatpk" =>
            include_bytes!("../../packs/guitar-rosewood-brazilian-fpl-gtr282-bremaud2011-shell.fsmatpk"),
        "guitar-cedar-western-red-fpl-gtr282-shell.fsmatpk" =>
            include_bytes!("../../packs/guitar-cedar-western-red-fpl-gtr282-shell.fsmatpk"),
        _ => return pack(name, 1.0, None),
    };
    NormalizedPack::from_bytes(bytes).unwrap()
}

/// The guitar E4 cost rung of bead .358: the Esteso body (`esteso_profile`) with 128 retained modes and its own cavity row
/// (13.04 L through the 43.5 mm sound hole, both ends flanged), and six strings (`GUITAR_STRINGS`) of 0.65 m speaking length
/// with a 20 mm afterlength. Each meets its bridge foot at the top node nearest its saddle point. The table's saddle line
/// x = 334.54 mm carries the strings 11.7 mm apart across the 58.5 mm of its outer drive points, the bass E on the -y side.
/// Each string has two planes and its twist and keeps every mode below `GUITAR_STRING_BAND_HZ`. String 5 (E4) is plucked:
/// held by 1 N at a fifth of its speaking length from the bridge in static equilibrium, then released at sample 0. No
/// cylinder gas loading (the 512-coordinate ceiling).
fn guitar_inputs(samples: u64) -> (Profile, input::Spec, super::super::shell::assembly::table::ConvertedTable) {
    use std::f64::consts::PI;
    let saddles: Vec<(f64, f64)> = (0..6).map(|k| (0.334_54, -0.029_25 + 0.0117 * k as f64)).collect();
    let (profile, converted) = esteso_profile(&saddles);
    let (speaking, span) = (0.65_f64, 0.67_f64);
    let bridge = speaking / span;
    let (density, young, shear) = (7800.0, 200e9, 80e9);
    let mut strings = String::new();
    for &(pitch, radius) in &GUITAR_STRINGS {
        let tension = density * PI * radius * radius * (2.0 * speaking * pitch).powi(2);
        let fundamental = 2.0 * speaking * pitch / (2.0 * span);
        let inharmonicity = PI * PI * young * PI * radius.powi(4) / 4.0 / (tension * span * span);
        let bending = (1..=input::MAX_MODES)
            .take_while(|&n| n as f64 * fundamental * (1.0 + inharmonicity * (n * n) as f64).sqrt() < GUITAR_STRING_BAND_HZ)
            .count();
        let twist = (GUITAR_STRING_BAND_HZ * 2.0 * span / (shear / density).sqrt()) as usize;
        strings.push_str(&format!("string wire.fsmatpk {span} {radius} tension {tension} {bending} {twist} {bridge} 3000 0.1\n\
            section bonded {} 1\nlayer sleeve.fsmatpk {radius}\npolarizations two 0.6 0.8 1 0\n", radius * 2.0 / 3.0));
    }
    let [volume, neck_radius, neck_length] = converted.cavity.unwrap();
    let pluck = bridge * 0.8;
    let performance = format!("frankensim-material-strings-v3\naudio 48000 {samples} 0.01 128\nstate 293.15 1000\n\
        gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2\n\
        body guitar-spruce-sitka-fpl-gtr282-shell.fsmatpk {} {} {ESTESO_REFERENCE_THICKNESS_M} 128 1\nstrings 6\n{strings}\
        prepared-linear 512 8589934592 1000000000000000\n\
        cavity physical-neck {volume} {neck_radius} {neck_length} flanged flanged 1\n\
        friction stribeck 0.8 0.4 0.1 0\nforces preload 1 1\nport 5 {pluck} centreline 1\nforce 0 0 0\n\
        schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\ntrack\tbow\tbow\t0\n\
        initial\tbow\t0e0\t0e0\t2e-1\nevents\t0\n", converted.length_m, converted.width_m);
    let spec = input::Spec::parse(performance.as_bytes()).unwrap();
    (profile, spec, converted)
}

/// The guitar E4 rung (above) for a quarter second: the wall time per audio second with what it was measured on, the time of
/// each stage, the bridge snap distances and the pressure bits.
#[test]
#[ignore = "guitar E4 cost rung of bead .358: run with --ignored --nocapture"]
fn the_esteso_guitar_e4_rung_reports_its_wall_time_per_audio_second() {
    let (profile, spec, converted) = guitar_inputs(12_000);
    eprintln!("guitar E4 rung: Esteso h12 as a v4 profile of {} bytes, body {:.4} x {:.4} m; bridge snaps {:?} mm",
        converted.profile.len(), converted.length_m, converted.width_m,
        converted.bridge_snaps_m.iter().map(|d| (d * 1e4).round() / 10.0).collect::<Vec<_>>());
    scoped(|cx, gate| {
        let constructing = std::time::Instant::now();
        let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(esteso_pack(name))).unwrap();
        let built = constructing.elapsed();
        eprintln!("guitar E4 rung: {}", cost_context(&r.spec));
        let source = r.performance.bank().source().bank();
        assert_eq!(source.strings().len(), 6);
        assert!(source.material_cavity().is_some());
        let rendering = std::time::Instant::now();
        let samples = pressure(&mut r, 128);
        let rendered = rendering.elapsed();
        assert_eq!(samples.len(), 12_000);
        assert!(samples.iter().all(|p| p.is_finite()));
        assert!(samples.iter().any(|p| p.abs() > 1e-14), "the plucked body radiated nothing");
        let per_second = rendered.as_secs_f64() / (samples.len() as f64 / 48_000.0);
        eprintln!("guitar E4 rung: constructed in {built:.1?}, {} samples rendered in {rendered:.1?} \
            ({per_second:.1} s render per audio second; a 2 s note with its construction takes {:.1} s per audio second)",
            samples.len(), (built.as_secs_f64() + 2.0 * per_second) / 2.0);
        eprintln!("guitar E4 rung: {}", cost_context(&r.spec));
        eprintln!("guitar E4 rung: pressure bits fnv1a {:#018x}, peak {:.3e} Pa", pressure_bits(&samples),
            samples.iter().fold(0.0_f64, |m, p| m.max(p.abs())));
    });
}

/// Cost probe, not a correctness check (`-- --ignored --nocapture`): prints how the prepared route's construction time and time
/// per sample grow with the retained body modes of the same plate, so the growth is measured rather than assumed.
#[test]
#[ignore = "timing probe: run with --ignored --nocapture"]
fn prepared_body_cost_against_retained_modes() {
    for modes in [10, 20, 40, 70] {
        let (profile, spec) = ceiling_inputs(modes, 28);
        scoped(|cx, gate| {
            let constructing = std::time::Instant::now();
            let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(if name == "body.fsmatpk" {
                pack_with_shear(name, 1.0, 1.0, None) } else { pack(name, 1.0, None) })).unwrap();
            let built = constructing.elapsed();
            let rendering = std::time::Instant::now();
            let samples = pressure(&mut r, 7);
            eprintln!("cost probe: {modes} body modes: constructed in {built:.1?}, {:.1?} per sample",
                rendering.elapsed() / samples.len() as u32);
        });
    }
}
