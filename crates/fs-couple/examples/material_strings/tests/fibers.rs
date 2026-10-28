//! A material pack and geometry record reach the existing physical WAV path.
use super::*;
use fs_material::section::axial::YOUNG_VISCOSITY_PROPERTY;

fn fiber_pack(rho: f64, omit_axial_loss: bool) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, rho),
        ("young_modulus", Pressure::DIMS, 3e9),
        (ETA, DynViscosity::DIMS, 50.0),
        (YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS, 50.0),
    ] {
        if omit_axial_loss && key == YOUNG_VISCOSITY_PROPERTY {
            continue;
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, 200_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured fiber integration fixture; not measured horsehair".into(),
                    license: "CC0-1.0".into(),
                    artifact: Some(fs_blake3::ContentHash([11; 32])),
                },
            })
            .unwrap();
    }
    let pack = NormalizedPack::new(
        "synthetic-fiber",
        "fixture",
        fs_blake3::ContentHash([11; 32]),
        "manufactured freely redistributable data",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

fn fiber_spec() -> input::Spec {
    input::Spec::parse(
        include_str!("../fiber-ribbon.material-strings")
            .replace("48000 4096", "48000 256")
            .as_bytes(),
    )
    .unwrap()
}

fn ready(cx: &Cx<'_>, gate: &CancelGate, density: f64) -> Ready {
    fiber_spec()
        .compile(cx, gate, |name| {
            Ok(if name == "fiber.fsmatpk" {
                fiber_pack(density, false)
            } else {
                pack(name, 7800.0, 200_000.0, None)
            })
        })
        .unwrap()
}

fn accelerating_fiber_text() -> String {
    include_str!("../fiber-ribbon.material-strings")
        .replace("48000 4096", "48000 256")
        .replace("initial\tbow\t1e-2\t", "initial\tbow\t0\t")
        .replace("forces retain", concat!(
            "bow_fiber_program 0 0.1 20 0.01 1048576 3\n",
            "fiber_speed 64 0.01\nfiber_speed 128 -0.01\nfiber_speed 192 0\n",
            "forces retain"))
}

#[test]
fn source_declared_fiber_support_acceleration_reversal_and_material_reach_wav() {
    scoped(|cx,gate| {
        let make=|density| input::Spec::parse(accelerating_fiber_text().as_bytes()).unwrap()
            .compile(cx,gate,|name| Ok(if name=="fiber.fsmatpk" {fiber_pack(density,false)}
                else {pack(name,7800.0,200_000.0,None)})).unwrap();
        let mut actual=make(1300.0);
        let mut dense=make(5200.0);
        assert_eq!(actual.performance.fiber_support_program().unwrap().knots.len(),3);
        let expected=pressure(&mut actual,37);
        let changed=pressure(&mut dense,128);
        assert!(expected.iter().any(|v|v.abs()>1e-12));
        assert!(expected.iter().zip(changed).any(|(a,b)|(a-b).abs()>1e-12));
        let frame=actual.performance.fiber_support_frame().unwrap();
        assert_eq!(frame.sample,256);
        assert_eq!(frame.velocity_m_s,0.0);
        assert!(frame.laboratory_mechanical_energy_j>0.0);
        let (encoded,clips)=fs_couple::pcm_wav::encode_pcm16_wav(&expected,48000,0.01).unwrap();
        let mut streamed=make(1300.0);
        let (output,summary)=streamed.write(Cursor::new(Vec::<u8>::new()),gate).unwrap();
        assert_eq!(output.into_inner(),encoded);
        assert_eq!(summary.clipped_samples,clips as u64);
        assert_eq!(actual.performance.fiber_support_frame(),streamed.performance.fiber_support_frame());
    });
}

#[test]
fn unsupported_support_motion_refuses_before_source_io() {
    let text=accelerating_fiber_text();
    for (from,to) in [
        ("bow_fiber_program 0 0.1 20 0.01 1048576 3", "bow_fiber_program 0 0.1 1 0.01 1048576 3"),
        ("fiber_speed 128 -0.01", "fiber_speed 64 -0.01"),
        ("fiber_speed 192 0", "fiber_speed 257 0"),
        ("initial\tbow\t0\t", "initial\tbow\t1e-2\t"),
        ("bow_fiber_program 0", "bow_fiber_motion 1 1000 0.001 1048576\nbow_fiber_program 0"),
        ("bow_fibers fiber.fsmatpk 0.00004 0.6 180 40 0.4 2 2\n", ""),
    ] {
        assert!(input::Spec::parse(text.replace(from,to).as_bytes()).is_err(),"{from}");
    }
}

#[test]
fn material_fiber_input_changes_the_actual_shared_pressure_and_wav() {
    scoped(|cx, gate| {
        let mut light = ready(cx, gate, 1300.0);
        let mut dense = ready(cx, gate, 5200.0);
        assert_eq!(light.catalog.len(), 3);
        assert_eq!(light.spec.bytes, dense.spec.bytes);
        assert_eq!(
            light.performance.bank().system().components()[0].modes(),
            dense.performance.bank().system().components()[0].modes()
        );
        assert_eq!(
            light
                .performance
                .bank()
                .source()
                .fiber_ribbon()
                .unwrap()
                .model()
                .resolution()
                .cells,
            3
        );
        let expected = pressure(&mut light, 37);
        let changed = pressure(&mut dense, 128);
        assert!(
            expected
                .iter()
                .zip(changed)
                .any(|(a, b)| (a - b).abs() > 1e-12)
        );
        let (encoded, clips) =
            fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48000, 0.01).unwrap();
        let mut actual = ready(cx, gate, 1300.0);
        let (output, summary) = actual.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(output.into_inner(), encoded);
        assert_eq!(summary.samples, 256);
        assert_eq!(summary.clipped_samples, clips as u64);
    });
}

#[test]
fn missing_fiber_loss_and_unresolved_or_accelerating_input_refuse() {
    scoped(|cx, gate| {
        assert!(
            fiber_spec()
                .compile(cx, gate, |name| Ok(if name == "fiber.fsmatpk" {
                    fiber_pack(1300.0, true)
                } else {
                    pack(name, 7800.0, 200_000.0, None)
                }))
                .is_err()
        );
        let mut unresolved = fiber_spec();
        unresolved.bow_fibers.as_mut().unwrap().axial_modes = 128;
        assert!(
            unresolved
                .compile(cx, gate, |name| Ok(if name == "fiber.fsmatpk" {
                    fiber_pack(1300.0, false)
                } else {
                    pack(name, 7800.0, 200_000.0, None)
                }))
                .is_err()
        );
    });
    let text = include_str!("../fiber-ribbon.material-strings");
    for (from, to) in [
        ("bow_ribbon 0.01 3\n", ""),
        ("0.6 180 40", "0.6 0 40"),
        ("40 0.4 2 2", "40 1.0 2 2"),
        ("40 0.4 2 2", "40 0.4 0 2"),
    ] {
        assert!(input::Spec::parse(text.replace(from, to).as_bytes()).is_err());
    }
    let schedule = GestureSchedule::try_new(
        700,
        vec![GestureTrack {
            id: "bow".into(),
            target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow {
                velocity_m_per_s: 0.01,
                normal_force_n: 0.01,
                station: 0.2,
            },
            events: vec![GestureEvent {
                time_s: 0.001,
                transition_s: 0.001,
                value: GestureValue::Bow {
                    velocity_m_per_s: 0.02,
                    normal_force_n: 0.01,
                    station: 0.2,
                },
            }],
        }],
    )
    .unwrap();
    let mut bytes =
        format!("{}\nschedule\n", text.split_once("\nschedule\n").unwrap().0).into_bytes();
    bytes.extend(schedule.to_canonical_bytes());
    assert!(input::Spec::parse(&bytes).is_err());
}


fn contact_spec() -> input::Spec {
    input::Spec::parse(
        include_str!("../fiber-contact.material-strings")
            .replace("48000 4096", "48000 256").as_bytes(),
    ).unwrap()
}

fn contact_ready(
    cx: &Cx<'_>, gate: &CancelGate, tension: f64, profile: &str,
) -> Ready {
    let mut spec = contact_spec();
    spec.bow_fibers.as_mut().unwrap().tension_n = tension;
    spec.compile_with_profiles(
        cx, gate, None, Some(crate::stops::Profile::parse(profile.as_bytes()).unwrap()),
        |name| Ok(if name == "fiber.fsmatpk" {
            fiber_pack(1300.0, false)
        } else {
            pack(name, 7800.0, 200_000.0, None)
        }),
    ).unwrap()
}

#[test]
fn fiber_normal_contact_reaches_the_same_scheduled_pressure_and_pcm_history() {
    scoped(|cx, gate| {
        let profile = include_str!("../fiber-contact.stops");
        let mut run = contact_ready(cx, gate, 40.0, profile);
        let source = run.performance.bank().source();
        assert_eq!(source.bow_contact_loads().len(), 3);
        let first = source.fiber_ribbon().unwrap().component_start();
        for (cell, (bow, stop)) in source.bows().iter().zip(source.stops()).enumerate() {
            assert_eq!(bow.location, stop.location);
            assert_eq!(stop.body.as_ref().unwrap().component + run.spec.strings.len(), first + cell);
        }
        assert!(run.performance.events().iter().all(|event| event.drive.normal_force_n == 0.0));
        let mut output = Vec::new();
        let mut loaded = false;
        for _ in 0..256 {
            let mut sample = [0.0];
            run.performance.block(&mut sample).unwrap();
            output.push(sample[0]);
            let system = run.performance.bank().system();
            let frame = system.last_frame().unwrap();
            for (port, normal) in frame.ports.iter().zip(system.last_obstacle_frames().unwrap()) {
                // The solver's actual reaction must arrive without being divided
                // by cell count a second time or replaced by a gesture load.
                assert_eq!(port.normal_force_n, normal.normal_force_n);
                loaded |= port.normal_force_n > 0.0;
                if port.normal_force_n == 0.0 { assert_eq!(port.traction_n, 0.0); }
            }
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(loaded && output.iter().any(|p| p.abs() > 1e-12));
        assert!(run.performance.bank().system().components()[first..].iter().any(|cell|
            cell.states()[2..].iter().any(|state| state.displacement_m_sqrt_kg != 0.0)));
        let mut block = contact_ready(cx, gate, 40.0, profile);
        assert_eq!(pressure(&mut block, 37), output);
        let mut wav = contact_ready(cx, gate, 40.0, profile);
        let (bytes, summary) = wav.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let (expected, clips) = fs_couple::pcm_wav::encode_pcm16_wav(
            &output, wav.spec.rate, wav.spec.full_scale,
        ).unwrap();
        assert_eq!(bytes.into_inner(), expected);
        assert_eq!(summary.clipped_samples, clips as u64);
    });
}

#[test]
fn moving_fiber_contact_patch_uses_source_materials_and_the_original_wav_path() {
    scoped(|cx, gate| {
        let text = include_str!("../fiber-contact.material-strings")
            .replace("48000 4096", "48000 256")
            .replace("forces retain 0 0", concat!(
                "bow_station_motion 2 0.00003 1048576 2\n",
                "bow_station_knot 64 0.202\n",
                "bow_station_knot 128 0.201\n",
                "forces retain 0 0"));
        let profile = include_str!("../fiber-contact.stops");
        let build = |density| input::Spec::parse(text.as_bytes()).unwrap().compile_with_profiles(
            cx, gate, None, Some(crate::stops::Profile::parse(profile.as_bytes()).unwrap()),
            |name| Ok(if name == "fiber.fsmatpk" { fiber_pack(density, false) }
                else { pack(name, 7800.0, 200_000.0, None) }),
        ).unwrap();
        let mut moving = build(1300.0);
        let first_fiber = moving.performance.bank().source().fiber_ribbon().unwrap().component_start();
        let initial:Vec<_> = moving.performance.bank().source().bows().iter()
            .map(|bow| bow.location.station_fraction).collect();
        let mut output = Vec::new();
        let (mut loaded, mut attachment_work) = (false, 0.0);
        for _ in 0..256 {
            let mut sample = [0.0];
            moving.performance.block(&mut sample).unwrap();
            output.push(sample[0]);
            let bank = moving.performance.bank();
            let bow_stations = bank.bow_station_fractions().unwrap();
            let stop_stations = bank.stop_station_fractions().unwrap();
            let system = bank.system();
            let frame = system.last_frame().unwrap();
            for binding in bank.source().bow_contact_loads() {
                assert_eq!(bow_stations[binding.bow], stop_stations[binding.stop]);
                let normal = system.last_obstacle_frames().unwrap()[binding.stop];
                assert_eq!(frame.ports[binding.bow].normal_force_n, normal.normal_force_n);
                attachment_work += normal.attachment_work_j.abs();
                loaded |= normal.normal_force_n > 0.0;
            }
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        let bank = moving.performance.bank();
        assert!(loaded && attachment_work > 0.0 && output.iter().any(|p| p.abs() > 1e-12));
        for (&station, &original) in bank.bow_station_fractions().unwrap().iter().zip(&initial) {
            assert!((station-original-0.001).abs() < 1e-15);
        }
        assert!(bank.system().components()[first_fiber..].iter().any(|fiber|
            fiber.states().iter().any(|s| s.displacement_m_sqrt_kg != 0.0)));
        // Changing only the fiber material reaches the same continuing coupled
        // pressure path. The stationary and migrating contact are distinguishable.
        let dense = pressure(&mut build(5200.0), 128);
        let fixed = pressure(&mut contact_ready(cx, gate, 40.0, profile), 128);
        assert!(output.iter().zip(dense).any(|(a,b)| (a-b).abs() > 1e-12));
        assert!(output.iter().zip(fixed).any(|(a,b)| (a-b).abs() > 1e-12));
        assert_eq!(output, pressure(&mut build(1300.0), 37));
        let mut wav = build(1300.0);
        let (bytes, summary) = wav.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let (expected, clips) = fs_couple::pcm_wav::encode_pcm16_wav(
            &output, wav.spec.rate, wav.spec.full_scale,
        ).unwrap();
        assert_eq!(bytes.into_inner(), expected);
        assert_eq!(summary.clipped_samples, clips as u64);
    });
}

#[test]
fn accelerating_fibers_keep_solved_normal_loads_through_input_and_pcm_output() {
    scoped(|cx, gate| {
        let text = include_str!("../fiber-contact.material-strings")
            .replace("48000 4096", "48000 256")
            .replace("initial\tbow\t1e-2\t", "initial\tbow\t0\t")
            .replace("forces retain", concat!(
                "bow_fiber_program 0 0.1 20 0.01 1048576 3\n",
                "fiber_speed 64 0.01\nfiber_speed 128 -0.01\nfiber_speed 192 0\n",
                "forces retain"));
        let profile = include_str!("../fiber-contact.stops");
        let build = |density| input::Spec::parse(text.as_bytes()).unwrap().compile_with_profiles(
            cx, gate, None, Some(crate::stops::Profile::parse(profile.as_bytes()).unwrap()),
            |name| Ok(if name == "fiber.fsmatpk" { fiber_pack(density, false) }
                else { pack(name, 7800.0, 200_000.0, None) }),
        ).unwrap();
        let mut run = build(1300.0);
        let mut output = Vec::new();
        let mut normal_work = 0.0;
        let mut loaded = false;
        for _ in 0..256 {
            let mut sample = [0.0];
            run.performance.block(&mut sample).unwrap();
            output.push(sample[0]);
            normal_work += run.performance.fiber_support_frame().unwrap().normal_attachment_work_j.abs();
            let system = run.performance.bank().system();
            for (bow, normal) in system.last_frame().unwrap().ports.iter()
                .zip(system.last_obstacle_frames().unwrap()) {
                assert_eq!(bow.normal_force_n, normal.normal_force_n);
                loaded |= bow.normal_force_n > 0.0;
            }
        }
        assert!(loaded && normal_work > 0.0 && output.iter().any(|x| x.abs() > 1e-12));
        assert_eq!(output, pressure(&mut build(1300.0), 37));
        assert!(output.iter().zip(pressure(&mut build(5200.0), 128)).any(|(a, b)| (a-b).abs() > 1e-12));
        let mut wav = build(1300.0);
        let (bytes, summary) = wav.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let (expected, clips) = fs_couple::pcm_wav::encode_pcm16_wav(
            &output, wav.spec.rate, wav.spec.full_scale,
        ).unwrap();
        assert_eq!(bytes.into_inner(), expected);
        assert_eq!(summary.clipped_samples, clips as u64);
    });
}

#[test]
fn fiber_tension_and_surface_clearance_change_load_and_sound_without_rewriting_physics() {
    scoped(|cx, gate| {
        let profile = include_str!("../fiber-contact.stops");
        let mut loose = contact_ready(cx, gate, 10.0, profile);
        let mut taut = contact_ready(cx, gate, 100.0, profile);
        let first = loose.performance.bank().source().fiber_ribbon().unwrap().component_start();
        // Tension changes the transverse family, not axial elasticity, the
        // original string/body modes, or the acoustic observation.
        for (a, b) in loose.performance.bank().system().components()[..first].iter()
            .zip(&taut.performance.bank().system().components()[..first]) {
            assert_eq!(a.modes(), b.modes());
        }
        assert_eq!(
            &loose.performance.bank().system().components()[first].modes()[..2],
            &taut.performance.bank().system().components()[first].modes()[..2],
        );
        let (mut normal_change, mut pressure_change) = (0.0, 0.0);
        for _ in 0..256 {
            let (mut a, mut b) = ([0.0], [0.0]);
            loose.performance.block(&mut a).unwrap();
            taut.performance.block(&mut b).unwrap();
            pressure_change += (a[0] - b[0]).abs();
            let a = loose.performance.bank().system().last_frame().unwrap();
            let b = taut.performance.bank().system().last_frame().unwrap();
            normal_change += a.ports.iter().zip(&b.ports)
                .map(|(a, b)| (a.normal_force_n - b.normal_force_n).abs()).sum::<f64>();
        }
        assert!(normal_change > 1e-10 && pressure_change > 1e-12);
        // A positive physical clearance with both bodies at rest has no normal
        // force and hence no bow excitation, despite the unchanged bow speed.
        let mut separated = contact_ready(cx, gate, 40.0, &profile.replace("-0.000001", "0.00003"));
        assert!(pressure(&mut separated, 128).iter().all(|value| *value == 0.0));
        assert!(separated.performance.bank().system().last_frame().unwrap().ports.iter()
            .all(|port| port.normal_force_n == 0.0 && port.traction_n == 0.0));
    });
}

#[test]
fn incompatible_fiber_contact_geometry_and_motion_refuse_before_material_io() {
    let profile = include_str!("../fiber-contact.stops");
    scoped(|cx, gate| {
        for (text, spec) in [
            (profile.replace("fiber:0", "fiber:2"), contact_spec()),
            (profile.replace("ribbon:2", "ribbon:3"), contact_spec()),
            (profile.replace("bow-contact 0 centre\n", ""), contact_spec()),
            (format!("{profile}motion 0.001 0.1 1\nknot 128 -0.000001 0 0\n"), contact_spec()),
            (profile.to_owned(), {
                let mut spec = contact_spec(); spec.bow_fibers = None; spec
            }),
            (profile.to_owned(), {
                let text = include_str!("../fiber-contact.material-strings")
                    .replace("1e-2\t0e0\t2e-1", "1e-2\t1e-3\t2e-1");
                input::Spec::parse(text.as_bytes()).unwrap()
            }),
        ] {
            let profile = crate::stops::Profile::parse(text.as_bytes()).unwrap();
            let mut calls = 0;
            let result = spec.compile_with_profiles(cx, gate, None, Some(profile), |name| {
                calls += 1;
                Ok(pack(name, 7800.0, 200_000.0, None))
            });
            assert!(result.is_err());
            assert_eq!(calls, 0);
        }
    });
}

fn nonlinear_fiber_text(axial: bool) -> String {
    let source = if axial {
        include_str!("../longitudinal-cylinder-bow.material-strings")
    } else {
        include_str!("../nonlinear-cylinder-bow.material-strings")
    };
    source.replace("48000 4096", "48000 129")
        .replace("prepared-linear 64 67108864 200000000000\n", "")
        .replace("cylinder-gas-loading 0.001 30 4", "cylinder-gas-loading 0.001 30 2")
        .replace("polarizations two 0.6 0.8 1 0\n",
            "polarizations two 0.6 0.8 1 0\nbending-loss prony 1 fiber_test_delta fiber_test_time\n")
        .replace("friction coulomb 0.8 0.4\n",
            "friction coulomb 0.8 0.4\nbow_ribbon 0.01 3\nbow_fibers fiber.fsmatpk 0.00004 0.6 180 40 0.4 2 2\n")
}

fn nonlinear_pack(name: &str, fiber_density: f64) -> NormalizedPack {
    if name == "fiber.fsmatpk" { return fiber_pack(fiber_density, false); }
    let layered = matches!(name, "wire.fsmatpk" | "sleeve.fsmatpk");
    let original = pack(name, 7800.0, 200_000.0, layered.then_some(ETA));
    if !layered { return original; }
    let mut claims = original.claims().clone();
    for (key, dims, value) in [
        ("equilibrium_young_modulus", Pressure::DIMS, 200e9),
        ("fiber_test_delta", Pressure::DIMS, 20e9),
        ("fiber_test_time", Dims([0, 0, 1, 0, 0, 0]), 0.001),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0)
                .with("omega", 0.0, 200_000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "manufactured nonlinear fiber composition".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    NormalizedPack::new(original.pack_id(), "nonlinear-fiber-fixture", original.source_artifact(),
        original.redistribution_terms(), claims, vec![], vec![]).unwrap()
}

#[test]
fn fibers_share_nonlinear_strain_layer_histories_gas_and_original_pcm() {
    scoped(|cx, gate| {
        for axial in [false, true] {
            let text = nonlinear_fiber_text(axial);
            let build = |rho| input::Spec::parse(text.as_bytes()).unwrap()
                .compile(cx, gate, |name| Ok(nonlinear_pack(name, rho))).unwrap();
            let mut direct = build(1300.0);
            let bank = direct.performance.bank();
            let first_fiber = bank.source().fiber_ribbon().unwrap().component_start();
            let source = bank.source().bank();
            assert!(source.cylinder_gas_loading().is_some());
            assert!(source.material_cavity().is_some());
            assert_eq!(source.geometric_stretching().is_some(), !axial);
            assert_eq!(source.longitudinal_waves().is_some(), axial);
            assert_eq!(bank.system().material_memory().len(), 4);
            let axial_component = source.longitudinal_waves().map(|waves|
                waves.strings()[0].axial_component());
            if let Some(index) = axial_component { assert!(index < first_fiber); }
            let mut dense = build(5200.0);
            for (a, b) in bank.system().components()[..first_fiber].iter()
                .zip(&dense.performance.bank().system().components()[..first_fiber]) {
                assert_eq!(a.modes(), b.modes(), "fiber substitution changed the original string/body basis");
            }
            let expected = pressure(&mut direct, 1);
            assert!(expected.iter().any(|p| p.abs() > 1e-15));
            assert_eq!(expected, pressure(&mut build(1300.0), 37));
            assert!(expected.iter().zip(pressure(&mut dense, 128))
                .any(|(a, b)| (a-b).abs() > 1e-15));
            let system = direct.performance.bank().system();
            let frame = system.last_frame().unwrap();
            assert!(frame.nonlinear.is_some(), "appending fibers discarded the strain law");
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert!(system.material_memory().iter().any(|z| *z != 0.0));
            assert!(system.components()[first_fiber..].iter().any(|fiber|
                fiber.states()[..2].iter().any(|q| q.displacement_m_sqrt_kg != 0.0)));
            if let Some(index) = axial_component {
                assert!(system.components()[index].states().iter()
                    .any(|q| q.displacement_m_sqrt_kg != 0.0));
            } else {
                assert!(frame.nonlinear.unwrap().nonlinear_energy_j > 0.0);
            }
            let (wav, summary) = build(1300.0).write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
            let encoded = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48_000, 0.01).unwrap();
            assert_eq!(wav.into_inner(), encoded.0);
            assert_eq!(summary.samples, 129);
        }
    });
}

#[test]
fn fiber_normal_reactions_and_string_stretching_share_one_work_balance() {
    scoped(|cx, gate| {
        let text = include_str!("../fiber-contact.material-strings")
            .replace("48000 4096", "48000 129")
            .replace("cavity neck", "geometric-stretching 0.2 64 100000 1e-11 1e-10\ncavity neck");
        let build = || input::Spec::parse(text.as_bytes()).unwrap().compile_with_profiles(
            cx, gate, None,
            Some(crate::stops::Profile::parse(include_str!("../fiber-contact.stops").as_bytes()).unwrap()),
            |name| Ok(if name == "fiber.fsmatpk" { fiber_pack(1300.0, false) }
                else { pack(name, 7800.0, 200_000.0, None) }),
        ).unwrap();
        let mut direct = build();
        let initial = direct.performance.bank().system().total_energy_j().unwrap();
        assert!(initial > 0.0, "the original elastic interference must remain stored");
        let (mut work, mut loss, mut stretch) = (0.0, 0.0, 0.0_f64);
        let mut output = Vec::new();
        for _ in 0..129 {
            let mut sample = [0.0]; direct.performance.block(&mut sample).unwrap();
            output.push(sample[0]);
            let system = direct.performance.bank().system();
            let frame = system.last_frame().unwrap();
            for (port, contact) in frame.ports.iter().zip(system.last_obstacle_frames().unwrap()) {
                assert_eq!(port.normal_force_n, contact.normal_force_n);
                if contact.normal_force_n == 0.0 { assert_eq!(port.traction_n, 0.0); }
            }
            work += frame.external_work_j + frame.drive_work_j + frame.obstacle_work_j;
            loss += frame.network_dissipation_j + frame.friction_dissipation_j + frame.obstacle_dissipation_j;
            stretch = stretch.max(frame.nonlinear.unwrap().nonlinear_energy_j);
            assert!((system.total_energy_j().unwrap() + loss - work - initial).abs() < 1e-10);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(stretch > 0.0 && output.iter().any(|p| p.abs() > 1e-12));
        assert_eq!(output, pressure(&mut build(), 37));
    });
}

#[test]
fn explicit_fiber_strain_source_reaches_contact_pressure_and_the_same_wav() {
    scoped(|cx, gate| {
        let text = include_str!("../fiber-strain-contact.material-strings")
            .replace("48000 4096", "48000 129");
        let build = || input::Spec::parse(text.as_bytes()).unwrap().compile_with_profiles(
            cx, gate, None,
            Some(crate::stops::Profile::parse(include_str!("../fiber-contact.stops").as_bytes()).unwrap()),
            |name| Ok(if name == "fiber.fsmatpk" { fiber_pack(1300.0, false) }
                else { pack(name, 7800.0, 200_000.0, None) }),
        ).unwrap();
        let mut run = build();
        let source = run.performance.bank().source();
        assert!(source.bank().geometric_stretching().is_some());
        assert!(source.bank().material_cavity().is_some());
        let fibers = source.fiber_ribbon().unwrap();
        let first = fibers.component_start();
        let law = fibers.model().strain().unwrap();
        let declared = run.spec.bow_fibers.as_ref().unwrap().strain.unwrap();
        assert_eq!(law.limits(), declared.0);
        assert_eq!(law.config(), declared.1);
        let mut output = Vec::new();
        for _ in 0..129 {
            let mut sample = [0.0];
            run.performance.block(&mut sample).unwrap();
            output.push(sample[0]);
            let frame = run.performance.bank().system().last_frame().unwrap();
            assert!(frame.nonlinear.is_some());
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        let bank = run.performance.bank();
        let q: Vec<_> = bank.system().components()[first].states().iter()
            .map(|state| state.displacement_m_sqrt_kg).collect();
        let observed = bank.source().fiber_ribbon().unwrap().model().strain().unwrap()
            .observe(&q).unwrap();
        assert!(observed.mean_dynamic_strain > 0.0 && observed.strain_energy_j > 0.0);
        assert!(output.iter().any(|p| p.abs() > 1e-12));
        let (wav, summary) = build().write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let (expected, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&output, 48000, 0.01).unwrap();
        assert_eq!(wav.into_inner(), expected);
        assert_eq!(summary.samples, 129);
        assert_eq!(summary.clipped_samples, clips as u64);
    });
}

#[test]
fn fiber_strain_domain_resolution_and_joint_config_refuse_before_material_io() {
    let text = include_str!("../fiber-strain-contact.material-strings");
    for (from, to) in [
        ("bow_fiber_strain 0.2 0.03", "bow_fiber_strain 0.4 0.03"),
        ("bow_fiber_strain 0.2 0.03", "bow_fiber_strain 0.2 0.2"),
        ("bow_fiber_strain 0.2 0.03 64", "bow_fiber_strain 0.2 0.03 0"),
        ("bow_fiber_strain 0.2 0.03 64", "bow_fiber_strain 0.2 0.03 32"),
        ("40 0.4 2 2", "40 0.4 33 2"),
        ("bow_fibers fiber.fsmatpk 0.00004 0.6 180 40 0.4 2 2\n", ""),
    ] {
        assert!(input::Spec::parse(text.replace(from, to).as_bytes()).is_err());
    }
    scoped(|cx, gate| {
        let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
        spec.bow_fibers.as_mut().unwrap().strain.as_mut().unwrap().1.max_iterations = 63;
        let mut reads = 0;
        let result = spec.compile(cx, gate, |_| {
            reads += 1;
            Ok(fiber_pack(1300.0, false))
        });
        assert!(result.is_err());
        assert_eq!(reads, 0);
    });
}

#[test]
fn conservative_fiber_and_string_strain_share_force_loaded_initialization() {
    scoped(|cx,gate| {
        let text=include_str!("../fiber-strain-contact.material-strings")
            .replace("48000 4096","48000 129")
            .replace("forces retain 0 0","forces preload 1 1\nport 0 0.2 centreline 0.01\nforce 0 0 0");
        let build=||input::Spec::parse(text.as_bytes()).unwrap().compile(cx,gate,|name|Ok(
            if name=="fiber.fsmatpk" {fiber_pack(1300.0,false)} else {pack(name,7800.0,200_000.0,None)}
        )).unwrap();
        let mut run=build();let bank=run.performance.bank();
        assert!(bank.source().bank().geometric_stretching().is_some());
        assert!(bank.source().fiber_ribbon().unwrap().model().strain().is_some());
        let report=bank.system().nonlinear_preload_diagnostics().unwrap();
        assert!(report.force_residual<=report.force_tolerance);
        assert_eq!(run.performance.preload_energy_j(),Some(bank.system().total_energy_j().unwrap()));
        let expected=pressure(&mut run,37);
        assert!(expected.iter().any(|p|p.abs()>1e-12));
        assert_eq!(expected,pressure(&mut build(),128));
    });
}

fn moving_fiber_spec() -> input::Spec {
    let text=include_str!("../fiber-ribbon.material-strings")
        .replace("48000 4096","48000 256")
        .replace("\nforces ","\nbow_fiber_motion 1 1000 0.001 1048576\nforces ");
    let schedule=GestureSchedule::try_new(48000,vec![GestureTrack {
        id:"bow".into(),target:GestureTarget::BowStroke {string:0},
        initial:GestureValue::Bow {velocity_m_per_s:0.01,normal_force_n:0.01,station:0.2},
        events:vec![GestureEvent {time_s:0.001,transition_s:0.001,
            value:GestureValue::Bow {velocity_m_per_s:-0.01,normal_force_n:0.01,station:0.2}}],
    }]).unwrap();
    let mut bytes=format!("{}\nschedule\n",text.split_once("\nschedule\n").unwrap().0).into_bytes();
    bytes.extend(schedule.to_canonical_bytes());
    input::Spec::parse(&bytes).unwrap()
}
fn moving_fiber_ready(cx:&Cx<'_>,gate:&CancelGate) -> Ready {
    moving_fiber_spec().compile(cx,gate,|name| Ok(if name=="fiber.fsmatpk" {
        fiber_pack(1300.0,false)
    } else { pack(name,7800.0,200_000.0,None) })).unwrap()
}
#[test]
fn native_fiber_acceleration_and_reversal_reach_the_same_pressure_and_pcm_writer() {
    scoped(|cx,gate| {
        let mut run=moving_fiber_ready(cx,gate);
        let original=run.spec.schedule.to_canonical_bytes();
        let source=run.performance.bank().source().fiber_ribbon().unwrap().station_fraction();
        let expected=pressure(&mut run,37);
        assert!(expected.iter().any(|p|p.abs()>1e-12));
        let state=run.performance.bank().fiber_support_kinematics().unwrap();
        assert_eq!(state.velocity_m_s,-0.01);
        assert_ne!(state.station_fraction,source);
        assert_eq!(run.performance.schedule().to_canonical_bytes(),original);
        let mut encoded_run=moving_fiber_ready(cx,gate);
        let (encoded,summary)=encoded_run.write(Cursor::new(Vec::new()),gate).unwrap();
        let (direct,clips)=fs_couple::pcm_wav::encode_pcm16_wav(&expected,48000,0.01).unwrap();
        assert_eq!(encoded.into_inner(),direct);
        assert_eq!(summary.samples,256);
        assert_eq!(summary.clipped_samples,clips as u64);
    });
}
#[test]
fn source_motion_does_not_admit_untransported_normal_or_thermal_hair_patches() {
    scoped(|cx,gate| {
        for translating_on_string in [false,true] {
            let mut contact=contact_spec();
            if translating_on_string {
                let text=include_str!("../fiber-contact.material-strings")
                    .replace("48000 4096","48000 256")
                    .replace("\nforces ",concat!("\nbow_station_motion 2 0.00003 1048576 2\n",
                        "bow_station_knot 64 0.202\nbow_station_knot 128 0.201\nforces "));
                contact=input::Spec::parse(text.as_bytes()).unwrap();
            }
            contact.bow_fibers.as_mut().unwrap().motion=moving_fiber_spec().bow_fibers.unwrap().motion;
            let parsed=stops::Profile::parse(include_str!("../fiber-contact.stops").as_bytes());
            // String-side patch travel is admissible, but does not transport a
            // normal footprint along accelerating hair. Both orders retain the
            // actual source profile instead of a fabricated modal map.
            assert!(parsed.is_ok());
            let result=contact.compile_with_profiles(cx,gate,None,Some(parsed.unwrap()),|name|Ok(if name=="fiber.fsmatpk" {
                fiber_pack(1300.0,false)
            } else {pack(name,7800.0,200_000.0,None)}));
            assert!(result.is_err());
        }
    });
    let text=include_str!("../fiber-ribbon.material-strings");
    for limits in ["1 0 0.001 1048576","1 1000 0.001 0","NaN 1000 0.001 1048576"] {
        let text=text.replace("\nforces ",&format!("\nbow_fiber_motion {limits}\nforces "));
        assert!(input::Spec::parse(text.as_bytes()).is_err());
    }
}

#[test]
fn fiber_control_ticks_activate_targets_for_the_following_mechanical_endpoint() {
    scoped(|cx,gate| {
        let mut run=moving_fiber_ready(cx,gate);
        for sample in 0..97 {
            let GestureValue::Bow {velocity_m_per_s,..}=run.spec.schedule.sample_value("bow",sample).unwrap()
                else {panic!("original bow track")};
            run.performance.block(&mut [0.0]).unwrap();
            let kinematics=run.performance.bank().fiber_support_kinematics().unwrap();
            assert_eq!(kinematics.velocity_m_s,velocity_m_per_s);
            assert_eq!(run.performance.bank().system().samples_rendered(),sample+1);
            if sample==48 {assert_eq!(kinematics.velocity_m_s,0.01);}
            if sample==95 {assert!(kinematics.velocity_m_s> -0.01);}
            if sample==96 {
                assert_eq!(kinematics.velocity_m_s,-0.01);
                // Forty-nine initial intervals, then a 48-interval linear
                // reversal whose integrated support displacement is zero.
                assert!((kinematics.displacement_m-49.0/48000.0*0.01).abs()<1e-17);
            }
        }
    });
}
