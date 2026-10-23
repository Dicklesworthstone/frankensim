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
