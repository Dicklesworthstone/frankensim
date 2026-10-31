//! A material pack and geometry record reach the existing physical WAV path.
use super::*;
use fs_material::section::axial::YOUNG_VISCOSITY_PROPERTY;

fn fiber_pack(rho: f64, omit_axial_loss: bool) -> NormalizedPack {
    fiber_elastic_pack(rho, 3e9, omit_axial_loss)
}

fn fiber_elastic_pack(rho: f64, young: f64, omit_axial_loss: bool) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, rho),
        ("young_modulus", Pressure::DIMS, young),
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

fn memory_text() -> String {
    include_str!("../fiber-ribbon.material-strings")
        .replace("48000 4096", "48000 256")
        .replace("forces retain", concat!(
            "bow_fiber_memory 24 young_modulus 2 fiber_delta_fast fiber_tau_fast ",
            "fiber_delta_slow fiber_tau_slow\nforces retain"))
}

fn memory_pack(delta: f64, omit_time: bool) -> NormalizedPack {
    let mut claims = fiber_pack(1300.0, false).claims().clone();
    for (key, dims, value) in [
        ("fiber_delta_fast", Pressure::DIMS, delta),
        ("fiber_tau_fast", fs_qty::Time::DIMS, 0.0004),
        ("fiber_delta_slow", Pressure::DIMS, 0.5 * delta),
        ("fiber_tau_slow", fs_qty::Time::DIMS, 0.006),
    ] {
        if omit_time && key == "fiber_tau_slow" {
            continue;
        }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0).with("omega", 0.0, 200_000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured isotropic fiber relaxation; not measured horsehair".into(),
                license: "CC0-1.0".into(),
                artifact: Some(fs_blake3::ContentHash([12; 32])),
            },
        }).unwrap();
    }
    let pack = NormalizedPack::new("synthetic-relaxing-fiber", "fixture",
        fs_blake3::ContentHash([12; 32]), "manufactured freely redistributable data",
        claims, vec![], vec![]).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

#[test]
fn declared_fiber_spectrum_reaches_contact_memory_pressure_and_streamed_wav() {
    scoped(|cx, gate| {
        let make = |delta| input::Spec::parse(memory_text().as_bytes()).unwrap()
            .compile(cx, gate, |name| Ok(if name == "fiber.fsmatpk" {
                memory_pack(delta, false)
            } else { pack(name, 7800.0, 200_000.0, None) })).unwrap();
        let mut soft = make(0.6e9);
        let mut stiff = make(3e9);
        assert_eq!(soft.spec.bytes, stiff.spec.bytes);
        let model = soft.performance.bank().source().fiber_ribbon().unwrap().model();
        let memory = model.material_relaxation().unwrap();
        assert_eq!(memory.memory_states(), 24);
        assert_eq!(memory.source().properties()[0].branches.len(), 2);
        assert_eq!(memory.source().branches()[0].relaxation_time_s, 0.0004);
        assert_eq!(memory.source().branches()[1].relaxation_time_s, 0.006);
        let expected = pressure(&mut soft, 37);
        let changed = pressure(&mut stiff, 128);
        assert!(expected.iter().any(|v| v.abs() > 1e-12));
        assert!(expected.iter().zip(changed).any(|(a, b)| (a - b).abs() > 1e-12));
        let last = soft.performance.bank().system().last_frame().unwrap();
        assert!(last.energy_residual_j.abs() <= last.energy_tolerance_j);
        let (encoded, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48000, 0.01).unwrap();
        let mut streamed = make(0.6e9);
        let (output, summary) = streamed.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(output.into_inner(), encoded);
        assert_eq!(summary.clipped_samples, clips as u64);
    });
}

#[test]
fn missing_fiber_spectrum_or_insufficient_history_budget_cannot_render() {
    scoped(|cx, gate| {
        for (text, omit_time) in [
            (memory_text(), true),
            (memory_text().replace("bow_fiber_memory 24", "bow_fiber_memory 23"), false),
        ] {
            let result = input::Spec::parse(text.as_bytes()).unwrap()
                .compile(cx, gate, |name| Ok(if name == "fiber.fsmatpk" {
                    memory_pack(0.6e9, omit_time)
                } else { pack(name, 7800.0, 200_000.0, None) }));
            assert!(result.is_err());
        }
    });
    for (from, to) in [
        ("bow_fiber_memory 24", "bow_fiber_memory 0"),
        ("bow_fiber_memory 24", "bow_fiber_memory 257"),
        ("fiber_delta_slow fiber_tau_slow", "fiber_delta_fast fiber_tau_fast"),
    ] {
        assert!(input::Spec::parse(memory_text().replace(from, to).as_bytes()).is_err(), "{to}");
    }
}

fn strain_memory_text() -> String {
    include_str!("../fiber-strain-contact.material-strings")
        .replace("48000 4096", "48000 129")
        .replace("bow_fiber_strain", concat!(
            "bow_fiber_memory 30 young_modulus 2 fiber_delta_fast fiber_tau_fast ",
            "fiber_delta_slow fiber_tau_slow\nbow_fiber_strain"))
}

#[test]
fn selected_fiber_memory_acts_on_full_strain_through_contact_pressure_and_wav() {
    scoped(|cx, gate| {
        let text = strain_memory_text();
        let make = |delta| input::Spec::parse(text.as_bytes()).unwrap().compile_with_profiles(
            cx, gate, None,
            Some(crate::stops::Profile::parse(include_str!("../fiber-contact.stops").as_bytes()).unwrap()),
            |name| Ok(if name == "fiber.fsmatpk" { memory_pack(delta, false) }
                else { pack(name, 7800.0, 200_000.0, None) }),
        ).unwrap();
        let (mut soft, mut stiff) = (make(0.6e9), make(3e9));
        assert_eq!(soft.spec.bytes, stiff.spec.bytes);
        let model = soft.performance.bank().source().fiber_ribbon().unwrap().model();
        let memory = model.material_relaxation().unwrap();
        assert!(model.strain().is_some() && memory.has_finite_strain());
        // Three original cells, two material arms, each with three strain
        // histories (mean + two axial harmonics) and two bending histories.
        assert_eq!(memory.memory_states(), 30);
        assert_eq!(soft.performance.bank().system().material_memory().len(), 12);
        let initial: Vec<_> = soft.performance.bank().system().nonlinear_material_memory()
            .into_iter().map(<[f64]>::to_vec).collect();
        assert_eq!(initial.iter().map(Vec::len).sum::<usize>(), 18);
        let mut output = Vec::new();
        let mut pressure_change = 0.0;
        let mut contact_change = 0.0;
        for _ in 0..129 {
            let (mut a, mut b) = ([0.0], [0.0]);
            soft.performance.block(&mut a).unwrap();
            stiff.performance.block(&mut b).unwrap();
            output.push(a[0]);
            pressure_change += (a[0] - b[0]).abs();
            let a = soft.performance.bank().system().last_frame().unwrap();
            let b = stiff.performance.bank().system().last_frame().unwrap();
            contact_change += a.ports.iter().zip(&b.ports)
                .map(|(a, b)| (a.traction_n - b.traction_n).abs()
                    + (a.normal_force_n - b.normal_force_n).abs()).sum::<f64>();
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
        }
        assert!(pressure_change > 1e-12 && contact_change > 1e-12);
        assert!(soft.performance.bank().system().nonlinear_material_memory().iter()
            .zip(&initial).any(|(after, before)| after.iter().zip(before)
                .any(|(a, b)| (a - b).abs() > 1e-16)));
        let mut chunked = make(0.6e9);
        assert_eq!(output, pressure(&mut chunked, 37));
        assert_eq!(soft.performance.bank().system().nonlinear_material_memory(),
            chunked.performance.bank().system().nonlinear_material_memory());
        assert_eq!(soft.performance.bank().system().material_memory(),
            chunked.performance.bank().system().material_memory());
        let (wav, summary) = make(0.6e9).write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let (expected, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&output, 48000, 0.01).unwrap();
        assert_eq!(wav.into_inner(), expected);
        assert_eq!(summary.samples, 129);
        assert_eq!(summary.clipped_samples, clips as u64);
    });
}

#[test]
fn full_fiber_strain_memory_requires_its_geometric_mean_budget_and_original_sources() {
    let gate = CancelGate::new_clock_free();
    let text = strain_memory_text();
    let spec = input::Spec::parse(text.as_bytes()).unwrap();
    let (model, _) = compile_fiber_source(&spec, memory_pack(0.6e9, false), &gate).unwrap();
    assert_eq!(model.material_relaxation().unwrap().memory_states(), 30);
    for (text, omit_time) in [
        (text.replace("bow_fiber_memory 30", "bow_fiber_memory 29"), false),
        (text, true),
    ] {
        let spec = input::Spec::parse(text.as_bytes()).unwrap();
        assert!(compile_fiber_source(&spec, memory_pack(0.6e9, omit_time), &gate).is_err());
    }
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
    spec.bow_fibers.as_mut().unwrap().prestress = StringPrestress::FixedTension(tension);
    spec.compile_with_profiles(
        cx, gate, None, Some(crate::stops::Profile::parse(profile.as_bytes()).unwrap()),
        |name| Ok(if name == "fiber.fsmatpk" {
            fiber_pack(1300.0, false)
        } else {
            pack(name, 7800.0, 200_000.0, None)
        }),
    ).unwrap()
}

fn pressing_fiber_text() -> String {
    include_str!("../fiber-contact.material-strings")
        .replace("48000 4096", "48000 512")
        .replace("forces retain", concat!(
            "bow_fiber_motion 1 1000 0.001 1048576\n",
            "bow_fiber_normal_program 0 0.1 20 0.01 1048576 4\n",
            "fiber_normal_speed 64 -0.02\nfiber_normal_speed 128 0\n",
            "fiber_normal_speed 192 0.02\nfiber_normal_speed 256 0\n",
            "forces retain"))
}

#[test]
fn normal_support_pressing_changes_solved_load_and_material_pressure_through_pcm() {
    scoped(|cx, gate| {
        let text = pressing_fiber_text();
        let profile = include_str!("../fiber-contact.stops").replace("-0.000001", "0.000005");
        let build = |density| input::Spec::parse(text.as_bytes()).unwrap().compile_with_profiles(
            cx, gate, None, Some(crate::stops::Profile::parse(profile.as_bytes()).unwrap()),
            |name| Ok(if name == "fiber.fsmatpk" { fiber_pack(density, false) }
                else { pack(name, 7800.0, 200_000.0, None) }),
        ).unwrap();
        let mut run = build(1300.0);
        assert_eq!(run.performance.fiber_normal_support_program().unwrap().knots.len(), 4);
        let mut output = Vec::new();
        let (mut loaded, mut released, mut support_work) = (false, false, 0.0);
        for sample in 1..=512 {
            let mut pressure = [0.0];
            run.performance.block(&mut pressure).unwrap();
            output.push(pressure[0]);
            let normal_motion = run.performance.fiber_support_frame().unwrap().normal.unwrap();
            support_work += normal_motion.contact_translation_work_j.abs();
            let system = run.performance.bank().system();
            for (bow, contact) in system.last_frame().unwrap().ports.iter()
                .zip(system.last_obstacle_frames().unwrap()) {
                assert_eq!(bow.normal_force_n, contact.normal_force_n);
                loaded |= bow.normal_force_n > 0.0;
                released |= sample > 256 && bow.normal_force_n == 0.0;
                if bow.normal_force_n == 0.0 { assert_eq!(bow.traction_n, 0.0); }
            }
        }
        assert!(loaded && released && support_work > 0.0);
        assert!(output.iter().any(|p| p.abs() > 1e-12));
        assert_eq!(output, pressure(&mut build(1300.0), 37));
        assert!(output.iter().zip(pressure(&mut build(5200.0), 128))
            .any(|(a, b)| (a-b).abs() > 1e-12));
        let mut wav = build(1300.0);
        let (bytes, summary) = wav.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let (expected, clips) = fs_couple::pcm_wav::encode_pcm16_wav(
            &output, wav.spec.rate, wav.spec.full_scale,
        ).unwrap();
        assert_eq!(bytes.into_inner(), expected);
        assert_eq!(summary.clipped_samples, clips as u64);
        assert_eq!(run.performance.fiber_support_frame(), wav.performance.fiber_support_frame());
    });
}

#[test]
fn normal_support_input_requires_an_explicit_frame_and_finite_motion() {
    let text = pressing_fiber_text();
    for (from, to) in [
        ("bow_fiber_motion 1 1000 0.001 1048576\n", ""),
        ("bow_fiber_normal_program 0 0.1 20", "bow_fiber_normal_program 0 0.1 1"),
        ("fiber_normal_speed 128 0", "fiber_normal_speed 64 0"),
        ("fiber_normal_speed 256 0", "fiber_normal_speed 513 0"),
        ("fiber_normal_speed 192 0.02", "fiber_normal_speed 192 0.2"),
    ] {
        assert!(input::Spec::parse(text.replace(from, to).as_bytes()).is_err(), "{from}");
    }
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

fn support_spec(support: &str) -> input::Spec {
    input::Spec::parse(include_str!("../fiber-ribbon.material-strings")
        .replace("48000 4096", "48000 64")
        .replace("180 40 0.4", &format!("180 {support} 0.4")).as_bytes()).unwrap()
}

fn compile_fiber_source(spec: &input::Spec, source: NormalizedPack, gate: &CancelGate)
    -> Result<(fs_couple::fiber_ribbon::TensionedFiberRibbon, Catalog)>
{
    let mut catalog = Catalog::new();
    let model = spec.bow_fibers.as_ref().unwrap().compile(
        spec, &mut catalog, &mut |_| Ok(source.clone().into()), gate)?;
    Ok((model, catalog))
}

#[test]
fn fixed_fiber_extension_uses_material_stiffness_and_physical_count_not_ribbon_cells() {
    let gate = CancelGate::new_clock_free();
    let base = support_spec("fixed-extension 0.594 0.05");
    for (count, cells, young) in [(180, 3, 3e9), (360, 3, 3e9), (180, 5, 3e9), (180, 3, 6e9)] {
        let mut spec = base.clone();
        spec.bow_fibers.as_mut().unwrap().fiber_count = count;
        spec.bow_ribbon.as_mut().unwrap().1 = cells;
        let (model, catalog) = compile_fiber_source(&spec, fiber_elastic_pack(1300.0, young, false), &gate).unwrap();
        let area = core::f64::consts::PI * 0.00004_f64.powi(2);
        let strain = (0.6 - 0.594) / 0.594;
        let expected = count as f64 * young * area * strain;
        assert!((model.geometry().tension_n - expected).abs() < 1e-12 * expected);
        let mass = count as f64 * 1300.0 * area * 0.6;
        assert!((model.mass_kg() - mass).abs() < 1e-12 * mass);
        let k = core::f64::consts::PI / 0.6;
        let axial = (young / 1300.0).sqrt() * k;
        let normal = (young / 1300.0 * (strain * k*k + 0.00004_f64.powi(2) / 4.0 * k.powi(4))).sqrt();
        assert!((model.cell_model().modes()[0].angular_frequency_rad_s - axial).abs() < 1e-10);
        assert!((model.cell_model().modes()[2].angular_frequency_rad_s - normal).abs() < 1e-10);
        assert_eq!(model.axial_material().elastic().layers()[0].material.card_identity(),
            catalog["fiber.fsmatpk"].card.content_hash());
        assert_eq!(spec.bow_fibers.as_ref().unwrap().prestress,
            StringPrestress::FixedExtension { stress_free_length_m: 0.594, linear_strain_limit: 0.05 });
    }
    // The original numerical field remains an independent, maintained TOTAL force.
    for young in [3e9, 6e9] {
        let (model, _) = compile_fiber_source(&support_spec("40"),
            fiber_elastic_pack(1300.0, young, false), &gate).unwrap();
        assert_eq!(model.geometry().tension_n, 40.0);
    }
}

fn thermal_fiber_pack(factor: f64) -> NormalizedPack {
    use fs_material::state_point::{LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
        INVERSE_TEMPERATURE_DIMS};
    let mut claims = fiber_pack(1300.0, false).claims().clone();
    let provenance = Provenance {
        source: "manufactured piecewise fiber expansion; not measured horsehair".into(),
        license: "CC0-1.0".into(), artifact: Some(fs_blake3::ContentHash([27; 32])),
    };
    let observation = claims.register_observation(ObservationDataset {
        specimen: "manufactured fiber".into(), method: "analytic piecewise alpha(T)".into(),
        artifact: fs_blake3::ContentHash([27; 32]), caveats: "not measured material data".into(),
        provenance: provenance.clone(),
    }).unwrap();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ALPHA, INVERSE_TEMPERATURE_DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: fs_qty::Temperature::DIMS,
            knots: [(280.0,1e-5),(290.0,1e-5),(295.0,3e-5),(300.0,1e-5),(320.0,1e-5)]
                .map(|(t, alpha)| (t, factor * alpha)).to_vec(), dims: INVERSE_TEMPERATURE_DIMS },
        validity: ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,200_000.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![observation], provenance,
    }).unwrap();
    let pack = NormalizedPack::new("thermal-fiber", "manufactured-test-v1",
        fs_blake3::ContentHash([27;32]), "manufactured redistribution permitted",
        claims, vec![], vec![]).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

#[test]
fn source_fiber_expansion_changes_contact_load_and_sound_in_the_original_joint_owner() {
    scoped(|cx, gate| {
        let text = include_str!("../fiber-contact.material-strings")
            .replace("48000 4096", "48000 64")
            .replace("state 293.15", "state 300")
            .replace("180 40 0.4", "180 thermal-extension 0.594 0.05 290 0.4");
        let make = |factor| input::Spec::parse(text.as_bytes()).unwrap().compile_with_profiles(
            cx, gate, None, Some(crate::stops::Profile::parse(
                include_str!("../fiber-contact.stops").as_bytes()).unwrap()),
            |name| Ok(if name == "fiber.fsmatpk" { thermal_fiber_pack(factor) }
                else { pack(name, 7800.0, 200_000.0, None) })).unwrap();
        let mut low = make(1.0);
        let mut high = make(20.0);
        assert_eq!(low.spec.bytes, high.spec.bytes);
        assert_ne!(low.catalog["fiber.fsmatpk"].pack.content_hash(),
            high.catalog["fiber.fsmatpk"].pack.content_hash());
        for (run, factor) in [(&low, 1.0), (&high, 20.0)] {
            let model = run.performance.bank().source().fiber_ribbon().unwrap().model();
            let ea = 3e9 * core::f64::consts::PI * 0.00004_f64.powi(2);
            // The crossed alpha knot contributes 0.0002; endpoint alpha*dT
            // would give 0.0001 and therefore the wrong installed preload.
            let expected = 180.0 * ea * ((0.6 - 0.594) / 0.594 - factor * 0.0002);
            assert!((model.geometry().tension_n - expected).abs() < 1e-11);
            assert_eq!(run.spec.bow_fibers.as_ref().unwrap().thermal_reference_kelvin, Some(290.0));
            assert_eq!(model.axial_material().elastic().layers()[0].material.card_identity(),
                run.catalog["fiber.fsmatpk"].card.content_hash());
        }
        let first = low.performance.bank().source().fiber_ribbon().unwrap().component_start();
        for (a, b) in low.performance.bank().system().components()[..first].iter()
            .zip(&high.performance.bank().system().components()[..first]) {
            assert_eq!(a.modes(), b.modes(), "only the fiber expansion source changed");
        }
        let (mut load_change, mut pressure_change) = (0.0, 0.0);
        for _ in 0..64 {
            let (mut a, mut b) = ([0.0], [0.0]);
            low.performance.block(&mut a).unwrap();
            high.performance.block(&mut b).unwrap();
            pressure_change += (a[0] - b[0]).abs();
            let a = low.performance.bank().system().last_frame().unwrap();
            let b = high.performance.bank().system().last_frame().unwrap();
            load_change += a.ports.iter().zip(&b.ports)
                .map(|(a, b)| (a.normal_force_n - b.normal_force_n).abs()).sum::<f64>();
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
        }
        assert!(load_change > 1e-10 && pressure_change > 1e-12);
    });
}

#[test]
fn fiber_supports_refuse_missing_expansion_slack_strain_and_uninitialized_memory() {
    let raw = include_str!("../fiber-ribbon.material-strings");
    for support in ["fixed-extension 0.6 0.05", "fixed-extension 0.61 0.05",
        "fixed-extension 0.594 0.001", "fixed-extension 0.594 1",
        "thermal-extension 0.594 0.05 0", "thermal-extension NaN 0.05 290"] {
        assert!(input::Spec::parse(raw.replace("180 40 0.4", &format!("180 {support} 0.4")).as_bytes()).is_err(), "{support}");
    }
    let fixed = raw.replace("180 40 0.4", "180 fixed-extension 0.594 0.05 0.4")
        .replace("forces retain", "bow_fiber_memory 24 young_modulus 0\nforces retain");
    assert!(input::Spec::parse(fixed.as_bytes()).is_err(),
        "a maintained-force spectrum does not initialize fixed-extension preload history");
    let gate = CancelGate::new_clock_free();
    let mut thermal = support_spec("thermal-extension 0.594 0.05 290");
    thermal.temperature = 300.0;
    assert!(compile_fiber_source(&thermal, fiber_pack(1300.0, false), &gate).is_err());
    assert!(compile_fiber_source(&thermal, thermal_fiber_pack(100.0), &gate).is_err(),
        "thermal slack cannot become a clamped positive tension");
    assert!(compile_fiber_source(&thermal, thermal_fiber_pack(1000.0), &gate).is_err(),
        "the shared owner must admit thermal and elastic strain separately");
    thermal.bow_fibers.as_mut().unwrap().thermal_reference_kelvin = Some(270.0);
    assert!(compile_fiber_source(&thermal, thermal_fiber_pack(1.0), &gate).is_err());
}

fn quasistatic_text() -> String {
    include_str!("../compliant-ribbon.material-strings").replace("48000 4096", "48000 256")
}

fn quasistatic_pack(young: f64, viscosity: Option<f64>) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("young_modulus_static", Pressure::DIMS, Some(young)),
        ("density", Density::DIMS, Some(1100.0)),
        ("kelvin_voigt_young_viscosity_pw98_model", DynViscosity::DIMS, viscosity),
    ] {
        let Some(value) = value else { continue; };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with("frequency", 0.0, 0.0).with("relative_humidity_known", 0.0, 0.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "synthetic quasistatic fiber integration test".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    NormalizedPack::new("quasistatic-test", "fixture", fs_blake3::ContentHash([13; 32]),
        "manufactured freely redistributable data", claims, vec![], vec![]).unwrap()
}

fn compile_quasistatic_source(s: &input::Spec) -> fs_couple::bowed_string::runtime::reciprocal::bank::multiple::compliant_ribbon::HairRibbonProperties {
    let gate = CancelGate::new_clock_free();
    s.bow_compliance.as_ref().unwrap().compile(s, &mut Catalog::new(),
        &mut |_| Ok(quasistatic_pack(7e9, Some(1e6)).into()), &gate).unwrap()
}

#[test]
fn quasistatic_fiber_compliance_scales_with_material_geometry_not_quadrature() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::compliant_ribbon::CompliantBowRibbon;
    let s = input::Spec::parse(quasistatic_text().as_bytes()).unwrap();
    let ribbon = |s: &input::Spec| {
        let (width, cells) = s.bow_ribbon.unwrap();
        CompliantBowRibbon::new(StringBowLocation { string: 0, station_fraction: 0.2 },
            0.5, width, cells, compile_quasistatic_source(s)).unwrap()
    };
    let baseline = ribbon(&s);
    let (k, d) = (baseline.lumped_tangential_stiffness_n_m(), baseline.lumped_tangential_damping_kg_s());
    for (width, cells) in [(0.005, 1), (0.02, 3)] {
        let mut changed = s.clone();
        changed.bow_ribbon = Some((width, cells));
        let actual = ribbon(&changed);
        assert_eq!(actual.lumped_tangential_stiffness_n_m(), k);
        assert_eq!(actual.lumped_tangential_damping_kg_s(), d);
        assert!((actual.cell_tangential_stiffness_n_m() * cells as f64 - k).abs() < 1e-10);
    }
    let mut doubled = s.clone();
    doubled.bow_compliance.as_mut().unwrap().fiber_count *= 2;
    assert_eq!(ribbon(&doubled).lumped_tangential_stiffness_n_m(), 2.0 * k);
    assert_eq!(ribbon(&doubled).lumped_tangential_damping_kg_s(), 2.0 * d);
    let mut longer = s.clone();
    longer.bow_compliance.as_mut().unwrap().length_m *= 2.0;
    assert_eq!(ribbon(&longer).lumped_tangential_stiffness_n_m(), 0.5 * k);
    assert_eq!(ribbon(&longer).lumped_tangential_damping_kg_s(), 0.5 * d);
}

#[test]
fn quasistatic_material_compliance_reaches_shared_pressure_and_exact_wav_path() {
    scoped(|cx, gate| {
        let make = |young, viscosity| {
            input::Spec::parse(quasistatic_text().as_bytes()).unwrap().compile(cx, gate, |name| {
                Ok(if name.starts_with("bow-hair-") { quasistatic_pack(young, Some(viscosity)) }
                else { pack(name, 7800.0, 200_000.0, None) })
            }).unwrap()
        };
        let mut soft = make(7e9, 1e6);
        let mut stiff = make(28e9, 1e6);
        let mut viscous = make(7e9, 4e6);
        assert_eq!(soft.spec.bytes, stiff.spec.bytes);
        assert_eq!(soft.spec.bytes, viscous.spec.bytes);
        let baseline = soft.performance.bank().source().compliant_ribbon().unwrap();
        let changed = stiff.performance.bank().source().compliant_ribbon().unwrap();
        assert_eq!(changed.lumped_tangential_stiffness_n_m(), 4.0 * baseline.lumped_tangential_stiffness_n_m());
        assert_eq!(changed.lumped_tangential_damping_kg_s(), baseline.lumped_tangential_damping_kg_s());
        let expected = pressure(&mut soft, 37);
        let stiffness_pressure = pressure(&mut stiff, 128);
        let viscosity_pressure = pressure(&mut viscous, 128);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        assert!(expected.iter().zip(stiffness_pressure).any(|(a,b)| (a-b).abs() > 1e-12));
        assert!(expected.iter().zip(viscosity_pressure).any(|(a,b)| (a-b).abs() > 1e-12));
        let frame = soft.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.ports.iter().any(|p| p.series_drive_dissipation_j > 0.0));
        let (encoded, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48000, 0.01).unwrap();
        let (output, summary) = make(7e9, 1e6).write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(output.into_inner(), encoded);
        assert_eq!(summary.clipped_samples, clips as u64);
    });
}

#[test]
fn quasistatic_fibers_refuse_implicit_loss_or_conflicting_material_query() {
    let text = quasistatic_text();
    for (from, to) in [
        ("quasistatic-continuation", "measured-broadband"),
        ("bow_ribbon 0.01 3\n", ""),
        ("16 3 0.8", "4 3 0.8"),
        ("0.65 50 0.5", "0.65 0 0.5"),
        ("0.65 50 0.5", "0.65 50 1"),
        ("relative_humidity_known 0", "frequency 0"),
        ("density viscosity", "young_modulus_static viscosity"),
        ("forces retain", "bow_fibers fiber.fsmatpk 0.00004 0.6 180 40 0.4 2 2\nforces retain"),
    ] {
        assert!(input::Spec::parse(text.replace(from, to).as_bytes()).is_err(), "{to}");
    }
    let gate = CancelGate::new_clock_free();
    let original = input::Spec::parse(text.as_bytes()).unwrap();
    let resolve = |s: &input::Spec| s.bow_compliance.as_ref().unwrap().compile(s,
        &mut Catalog::new(), &mut |_| Ok(quasistatic_pack(7e9, None).into()), &gate);
    assert!(resolve(&original).is_err(), "missing viscosity became zero damping");
    let lossless = input::Spec::parse(text.replace(
        "viscosity kelvin_voigt_young_viscosity_pw98_model", "lossless-reference").as_bytes()).unwrap();
    assert_eq!(resolve(&lossless).unwrap().distributed_damping_kg_s_m, 0.0);
    let mut outside = lossless.clone();
    outside.bow_compliance.as_mut().unwrap().axes[0].1 = 1000.0;
    assert!(resolve(&outside).is_err(), "source frequency validity was silently widened");
    let mut conflict = original.clone();
    conflict.material_states.insert(conflict.bow_compliance.as_ref().unwrap().pack.clone(),
        vec![("relative_humidity_known".into(), 1.0)]);
    assert!(conflict.network_dimensions(false).is_err());
}

#[test]
fn quasistatic_hair_and_omitted_string_bending_compose_without_retuning() {
    scoped(|cx, gate| {
        let text = quasistatic_text().replace("cavity neck",
            "residual-flexibility exact-static-green\nbow-residual-flexibility exact-static-green\ncavity neck");
        let make = |young| input::Spec::parse(text.as_bytes()).unwrap().compile(cx, gate, |name| {
            Ok(if name.starts_with("bow-hair-") { quasistatic_pack(young, Some(1e6)) }
            else { pack(name, 7800.0, 200_000.0, None) })
        }).unwrap();
        let mut soft = make(7e9);
        let mut stiff = make(28e9);
        let source = soft.performance.bank().source();
        assert!(source.bank().residual_flexibility());
        assert!(source.compliant_ribbon().is_some());
        let tail = soft.performance.bank().system().friction_residual_compliance().unwrap();
        assert_eq!(tail.len(), 9);
        assert!(tail[0] > 0.0 && tail[1] != 0.0);
        assert_eq!(Some(tail), stiff.performance.bank().system().friction_residual_compliance());
        for (a, b) in soft.performance.bank().system().components().iter()
            .zip(stiff.performance.bank().system().components()) {
            assert_eq!(a.modes(), b.modes(), "hair data must not retune string/body modes");
        }
        let expected = pressure(&mut soft, 37);
        let changed = pressure(&mut stiff, 128);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        assert!(expected.iter().zip(changed).any(|(a, b)| (a-b).abs() > 1e-12));
        let frame = soft.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.ports.iter().any(|p| p.series_drive_dissipation_j > 0.0));
        assert!(frame.friction_residual_energy_j > 0.0);
        let (encoded, _) = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48000, 0.01).unwrap();
        let (output, _) = make(7e9).write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(output.into_inner(), encoded);
    });
}
