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
