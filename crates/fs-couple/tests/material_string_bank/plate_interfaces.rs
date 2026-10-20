//! Sourced joint data reaches the original shared string/body pressure solve.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::MaterialBridgeInterface;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::state_point::{
    ResolvedLinearInterface, interface_stiffness_per_area_dims, interface_viscosity_per_area_dims,
    resolve_interface_state_point,
};

fn interface(
    pair: [MaterialStateId; 2],
    k: f64,
    c: f64,
    t: f64,
    upper: f64,
) -> ResolvedLinearInterface {
    interface_with_band(pair, k, c, t, (0.0, upper))
}

fn interface_with_band(
    pair: [MaterialStateId; 2],
    k: f64,
    c: f64,
    t: f64,
    band: (f64, f64),
) -> ResolvedLinearInterface {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        ("joint-k", interface_stiffness_per_area_dims(), k),
        ("joint-c", interface_viscosity_per_area_dims(), c),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 310.0)
                    .with("omega", band.0, band.1),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source: "synthetic string/body interface".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
                observations: vec![],
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::NonNegative)
                .unwrap(),
        );
    }
    let [first, second] = pair;
    let surface = |material| SurfaceSpec {
        material,
        texture_frame: "declared contact axes".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface(first),
        surface(second),
        SystemContext {
            medium: "synthetic bilateral joint".into(),
            third_body: None,
            environment: "frozen state".into(),
            history: "source preparation retained".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let source = resolve_interface_state_point(
        &card,
        &QueryPoint::new()
            .with("T", t)
            .unwrap()
            .with("omega", band.0.max(1.0))
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedLinearInterface::new(source, "joint-k", "joint-c", 1e-4).unwrap()
}

fn sourced(
    cx: &Cx<'_>,
    gate: &CancelGate,
    k: f64,
    c: f64,
) -> (ReciprocalStringBank, ResolvedLinearInterface) {
    let mut body = source_plate(7800.0, 2);
    let pair = [
        specimen(7800.0, 200e9).material().material().clone(),
        body.specimen.material().material().clone(),
    ];
    let interface = interface(pair, k, c, 293.15, 100_000.0);
    body.bridges[0] = SharedPlateBridge::from_interface(0.8, 0.0, &interface).unwrap();
    let bank = ReciprocalStringBank::from_material_plate(
        cx,
        sources(),
        gas(101_325.0),
        body,
        law(),
        options(0),
        gate,
    )
    .unwrap();
    (bank, interface)
}
fn bind(
    bank: ReciprocalStringBank,
    interface: ResolvedLinearInterface,
    gate: &CancelGate,
) -> ReciprocalStringBank {
    bank.with_material_bridge_interfaces(
        vec![MaterialBridgeInterface {
            string: 0,
            interface,
        }],
        gate,
    )
    .unwrap()
}

#[test]
fn source_binding_preserves_live_state_and_material_data_changes_actual_pressure() {
    with_cx(|cx, gate| {
        let (mut raw, original) = sourced(cx, gate, 3e7, 1e3);
        let (bound, _) = sourced(cx, gate, 3e7, 1e3);
        let mut bound = bind(bound, original.clone(), gate);
        let (changed, changed_source) = sourced(cx, gate, 1.2e8, 4e3);
        let mut changed = bind(changed, changed_source, gate);
        assert_eq!(
            raw.cards()[0].fundamental_hz(),
            changed.cards()[0].fundamental_hz()
        );
        let zero = vec![0.0; raw.system().mode_count()];
        let mut pressure_delta = 0.0_f64;
        for _ in 0..128 {
            let a = *raw.step(&zero, drive()).unwrap();
            assert_eq!(a, *bound.step(&zero, drive()).unwrap());
            let b = *changed.step(&zero, drive()).unwrap();
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
            pressure_delta += (a.observer_pressure_pa - b.observer_pressure_pa).abs();
        }
        assert!(
            pressure_delta > 1e-12,
            "source bridge data must change actual microphone pressure"
        );
        let states: Vec<_> = raw
            .system()
            .components()
            .iter()
            .map(|c| c.states().to_vec())
            .collect();
        let frame = *raw.system().last_frame().unwrap();
        let bound_live = bind(raw, original.clone(), gate);
        assert_eq!(bound_live.system().samples_rendered(), 128);
        assert_eq!(bound_live.system().last_frame(), Some(&frame));
        assert_eq!(
            bound_live
                .system()
                .components()
                .iter()
                .map(|c| c.states().to_vec())
                .collect::<Vec<_>>(),
            states
        );
        assert_eq!(
            bound_live.source().material_bridge_interfaces()[0].interface,
            original
        );
    });
}

#[test]
fn exact_original_interface_survives_scheduled_performance() {
    with_cx(|cx, gate| {
        let (bank, original) = sourced(cx, gate, 3e7, 1e3);
        let bank = bind(bank, original.clone(), gate);
        assert_eq!(
            bank.source().bridge_surface_materials(0).unwrap(),
            *original.resolved().surface_materials()
        );
        let schedule = GestureSchedule::try_new(
            1000,
            vec![GestureTrack {
                id: "bow".into(),
                target: GestureTarget::BowStroke { string: 0 },
                initial: GestureValue::Bow {
                    velocity_m_per_s: 0.01,
                    normal_force_n: 0.1,
                    station: 0.2,
                },
                events: vec![],
            }],
        )
        .unwrap();
        let mut performance = bank
            .into_scheduled(
                &schedule,
                "bow",
                FrictionGestureConfig {
                    sample_rate_hz: 48_000,
                    samples: 32,
                    max_block: 32,
                    max_work: 1000,
                    max_events: 16,
                    station_fraction: 0.2,
                },
            )
            .unwrap();
        let mut pressure = [0.0; 32];
        performance.renderer_mut().block(&mut pressure).unwrap();
        assert!(pressure.iter().any(|p| *p != 0.0));
        assert_eq!(
            performance.source().material_bridge_interfaces()[0].interface,
            original
        );
    });
}

#[test]
fn cavity_neck_is_a_physical_coordinate_in_the_conservative_interface_band() {
    with_cx(|cx, gate| {
        for lower in [0.0, 1.0] {
            let mut body = source_plate(7800.0, 2);
            let original = interface_with_band(
                [
                    specimen(7800.0, 200e9).material().material().clone(),
                    body.specimen.material().material().clone(),
                ],
                3e7,
                1e3,
                293.15,
                (lower, 100_000.0),
            );
            body.bridges[0] = SharedPlateBridge::from_interface(0.8, 0.0, &original).unwrap();
            let bank = ReciprocalStringBank::from_material_plate_with_cavity(
                cx,
                sources(),
                gas(101_325.0),
                body,
                cavity(0.003, gas(101_325.0)),
                law(),
                options(0),
                gate,
            )
            .unwrap();
            assert_eq!(
                bank.system().components().last().unwrap().modes()[0].angular_frequency_rad_s,
                0.0
            );
            let result = bank.with_material_bridge_interfaces(
                vec![MaterialBridgeInterface {
                    string: 0,
                    interface: original,
                }],
                gate,
            );
            assert_eq!(result.is_ok(), lower == 0.0);
        }
    });
}

#[test]
fn mismatched_joint_source_surfaces_temperature_band_and_rebinding_refuse() {
    with_cx(|cx, gate| {
        for defect in 0..8 {
            let (bank, original) = sourced(cx, gate, 3e7, 1e3);
            let pair = bank.source().bridge_surface_materials(0).unwrap();
            let mut binding = MaterialBridgeInterface {
                string: 0,
                interface: original.clone(),
            };
            match defect {
                0 => binding.string = 2,
                1 => binding.interface = interface(pair, 3.1e7, 1e3, 293.15, 100_000.0),
                2 => {
                    let mut wrong = pair;
                    wrong[1].revision += 1;
                    binding.interface = interface(wrong, 3e7, 1e3, 293.15, 100_000.0);
                }
                3 => binding.interface = interface(pair, 3e7, 1e3, 300.0, 100_000.0),
                4 => binding.interface = interface(pair, 3e7, 1e3, 293.15, 10.0),
                5 => {
                    assert!(bank.with_material_bridge_interfaces(vec![], gate).is_err());
                    continue;
                }
                6 => {
                    assert!(
                        bank.with_material_bridge_interfaces(vec![binding.clone(), binding], gate)
                            .is_err()
                    );
                    continue;
                }
                _ => {
                    let free_upper = bank
                        .system()
                        .components()
                        .iter()
                        .flat_map(|c| c.modes())
                        .map(|m| m.angular_frequency_rad_s)
                        .fold(0.0_f64, f64::max);
                    let coupled_upper = bank.system().linear_frequency_upper_bound_rad_s();
                    assert!(coupled_upper > free_upper);
                    // Every individual mode fits; the physical joint stiffness
                    // still raises the complete-network applicability bound.
                    binding.interface =
                        interface(pair, 3e7, 1e3, 293.15, 0.5 * (free_upper + coupled_upper));
                }
            }
            assert!(
                bank.with_material_bridge_interfaces(vec![binding], gate)
                    .is_err(),
                "defect {defect}"
            );
        }
        let (bank, original) = sourced(cx, gate, 3e7, 1e3);
        assert!(
            bind(bank, original.clone(), gate)
                .with_material_bridge_interfaces(
                    vec![MaterialBridgeInterface {
                        string: 0,
                        interface: original.clone()
                    }],
                    gate
                )
                .is_err()
        );
        for (station, rest) in [
            (0.0, 0.0),
            (1.0, 0.0),
            (f64::NAN, 0.0),
            (0.8, f64::INFINITY),
        ] {
            assert!(SharedPlateBridge::from_interface(station, rest, &original).is_err());
        }
        let raw = super::super::build(cx, gate, 7800.0);
        assert!(raw.source().bridge_surface_materials(0).is_err());
        assert!(
            raw.with_material_bridge_interfaces(
                vec![MaterialBridgeInterface {
                    string: 0,
                    interface: original
                }],
                gate
            )
            .is_err()
        );
    });
}
