//! One manufactured two-string/two-layer fixture exercises actual native offsets.
//! No measured string heating or instrument validation is asserted.
use super::*;

const INITIAL: f64 = 293.15;
const COUNT: usize = 24;
const CORE: f64 = 0.0002;
const OUTER: f64 = 0.0003;
const LOSS: &str = "bending-loss prony 2 inactive_modulus inactive_time bending_delta bending_time";

fn volume() -> f64 {
    core::f64::consts::PI * (OUTER * OUTER - CORE * CORE) * 0.5
}

fn text(rates: bool) -> String {
    let option = if rates { " source-temperature" } else { "" };
    format!(
        "frankensim-contact-thermal-v4\nboundary adiabatic\n\
        solver 512 24 1e-9 1e-12 1e-12 4e-12 1e-13 1e-10 1 1e-8\n\
        cells 1\ncell 30 sleeve.fsmatpk 280 {INITIAL} 320 0 0 0 0.5 0 0 {:.17e} 0.001\n\
        dry-contacts 0\nfilms 0\nfilm-surfaces 0\nseries 0\n\
        material-heat-budget 256 67108864 100000000000\nfiber-heat 0\n\
        string-heat 1\nstring-body 1 1 30{option}\nlinks 0\n",
        volume()
    )
}

fn pack(name: &str, cp: f64, rates: bool) -> NormalizedPack {
    let base = material(name, cp, Some(ETA));
    let mut claims = base.claims().clone();
    for (name, dims, value) in [
        (ETA, DynViscosity::DIMS, 0.0),
        (
            fs_couple::string_specimen::EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            200e9,
        ),
        ("inactive_modulus", Pressure::DIMS, 0.0),
        ("inactive_time", fs_qty::Time::DIMS, 0.02),
        ("bending_delta", Pressure::DIMS, 80e9),
        ("bending_time", fs_qty::Time::DIMS, 0.0005),
    ] {
        let mut property = claim(
            name,
            dims,
            value,
            ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0)
                .with("omega", 0.0, 400_000.0),
        );
        if name == "bending_time" && rates {
            property.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(280.0, 0.0008), (320.0, 0.0002)],
                dims,
            };
            property.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(property).unwrap();
    }
    NormalizedPack::new(
        name,
        "manufactured-string-heat",
        fs_blake3::ContentHash([57; 32]),
        "manufactured CC0 integration data",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

fn spec(profile: &Profile, explicit: bool) -> input::Spec {
    let text = include_str!("../relaxing-layered-strings.material-strings")
        .replace(
            "bending-loss prony 1 relaxing_bending_modulus bending_relaxation_time",
            LOSS,
        )
        .replace(
            "layer relaxing-wire.fsmatpk 0.0003",
            "layer sleeve.fsmatpk 0.0003",
        )
        .replace("port 0 0.2 surface 0.002", "port 1 0.2 surface 0.1");
    let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
    spec.samples = COUNT as u64;
    spec.block = 7;
    for wire in &mut spec.strings {
        wire.bending_modes = 2;
        wire.twist_modes = 0;
    }
    if explicit {
        spec.strings[1].bending_loss = input::BendingLoss::KelvinVoigt;
        spec.bending_memory = Some(
            bending_memory::Input::new(
                16,
                vec![bending_memory::Selection {
                    string: 1,
                    layers: vec![
                        bending_memory::LayerProperties {
                            equilibrium_modulus: "young_modulus".into(),
                            branches: vec![
                                ("inactive_modulus".into(), "inactive_time".into()),
                                ("bending_delta".into(), "bending_time".into())
                            ]
                        };
                        2
                    ],
                }],
                &spec.strings,
            )
            .unwrap(),
        );
    }
    profile.configure(&mut spec);
    spec
}

fn ready(
    cx: &Cx<'_>,
    gate: &CancelGate,
    profile: &Profile,
    cp: f64,
    rates: bool,
    explicit: bool,
) -> Ready {
    spec(profile, explicit)
        .compile(cx, gate, |name| {
            Ok(if name == "body.fsmatpk" {
                material(name, 20.0, None)
            } else {
                pack(
                    name,
                    if name == "sleeve.fsmatpk" { cp } else { 20.0 },
                    rates,
                )
            })
        })
        .unwrap()
}

fn torsional_ready(
    cx: &Cx<'_>,
    gate: &CancelGate,
    profile: &Profile,
    cp: f64,
    explicit: bool,
) -> Ready {
    let mut specification = spec(profile, explicit);
    for wire in &mut specification.strings {
        wire.twist_modes = 1;
        wire.torsion_memory = Some(vec![
            torsion_memory::LayerProperties {
                equilibrium_modulus: SHEAR.into(),
                branches: vec![
                    ("inactive_modulus".into(), "inactive_time".into()),
                    ("shear_delta".into(), "shear_time".into())
                ],
            };
            2
        ]);
    }
    specification
        .compile(cx, gate, |name| {
            if name == "body.fsmatpk" {
                return Ok(material(name, 20.0, None));
            }
            let base = pack(name, if name == "sleeve.fsmatpk" { cp } else { 20.0 }, true);
            let mut claims = ClaimSet::new();
            for id in base.claims().observation_ids() {
                claims
                    .register_observation(base.claims().observation(id).unwrap().clone())
                    .unwrap();
            }
            for (_, property) in base.claims().claims_ordered() {
                if property.key.name() != SHEAR_ETA {
                    claims.insert_claim(property.clone()).unwrap();
                }
            }
            for (key, dims, value) in [
                (SHEAR_ETA, DynViscosity::DIMS, 0.0),
                ("shear_delta", Pressure::DIMS, 30e9),
                ("shear_time", fs_qty::Time::DIMS, 0.0005),
            ] {
                let mut property = claim(
                    key,
                    dims,
                    value,
                    ValidityDomain::unconstrained()
                        .with("T", 280.0, 320.0)
                        .with("omega", 0.0, 400_000.0),
                );
                if key == "shear_time" {
                    // Independent, oppositely directed source mobility. This is
                    // manufactured integration data, not a measured metal law.
                    property.value = PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![(280.0, 0.0002), (320.0, 0.0008)],
                        dims,
                    };
                    property.interpolation = InterpolationPolicy::LinearInside;
                }
                claims.insert_claim(property).unwrap();
            }
            Ok(NormalizedPack::new(
                name,
                "manufactured-bending-torsion-heat",
                fs_blake3::ContentHash([58; 32]),
                "manufactured CC0 integration data",
                claims,
                vec![],
                vec![],
            )
            .unwrap())
        })
        .unwrap()
}

#[test]
fn v4_string_heat_keeps_second_string_outer_layer_planes_and_source_arm_order() {
    scoped(|cx, gate| {
        let profile = Profile::parse(text(true).as_bytes()).unwrap();
        for (explicit, torsion) in [(false, false), (true, false), (false, true), (true, true)] {
            let ready = if torsion {
                torsional_ready(cx, gate, &profile, 0.001, explicit)
            } else {
                ready(cx, gate, &profile, 0.001, true, explicit)
            };
            let prepared = profile.prepare(&ready, gate).unwrap();
            let network = prepared.network(gate).unwrap();
            let (routes, rates) = prepared.material_routes(&ready, &network, gate).unwrap();
            assert_eq!(routes.len(), 1);
            let (bending, twist) = if torsion && explicit {
                (16..20, Some(11..12))
            } else {
                (12..16, torsion.then_some(19..20))
            };
            let mut groups = vec![bending.clone()];
            let mut times: Vec<_> = bending
                .map(|branch| (branch, "bending_time".to_owned()))
                .collect();
            if let Some(twist) = twist {
                groups.push(twist.clone());
                times.extend(twist.map(|branch| (branch, "shear_time".to_owned())));
            }
            groups.sort_unstable_by_key(|range| range.start);
            times.sort_unstable_by_key(|(branch, _)| *branch);
            assert_eq!(routes[0].channel.as_ref().unwrap().groups(), groups);
            assert_eq!(rates[0].times, times);
            assert!((routes[0].volume_m3 - volume()).abs() < 1e-14 * volume());
            assert_eq!(
                routes[0].source.card_identity(),
                prepared.sources[0].card.content_hash()
            );
            // Four channels are the two planes times two sine modes of the
            // positive arm, after an explicitly zero source arm was skipped.
            let branches = ready
                .performance
                .bank()
                .system()
                .material_relaxation_branches();
            assert_eq!(branches.len(), if torsion { 20 } else { 16 });
            for (branch, name) in &rates[0].times {
                assert_eq!(
                    branches[*branch].relaxation_time_s.to_bits(),
                    routes[0]
                        .source
                        .property(name)
                        .unwrap()
                        .value_si()
                        .to_bits()
                );
            }
        }
    });
}

#[test]
fn v4_string_heat_routes_bending_and_torsion_through_one_temperature_and_pcm() {
    scoped(|cx, gate| {
        let profile = Profile::parse(text(true).as_bytes()).unwrap();
        let mut outcomes = Vec::new();
        for cp in [0.001, 0.004] {
            let mut direct = torsional_ready(cx, gate, &profile, cp, false);
            let original_times: Vec<_> = direct
                .performance
                .bank()
                .system()
                .material_relaxation_branches()
                .iter()
                .map(|branch| branch.relaxation_time_s)
                .collect();
            let prepared = profile.prepare(&direct, gate).unwrap();
            let network = prepared.network(gate).unwrap();
            let (routes, rates) = prepared.material_routes(&direct, &network, gate).unwrap();
            let mut renderer = direct
                .performance
                .with_material_contact_thermal_network(
                    network,
                    vec![],
                    vec![],
                    routes,
                    profile.controls,
                    gate,
                    |endpoint| prepared.resolve(endpoint),
                )
                .unwrap()
                .with_material_relaxation_temperature(rates)
                .unwrap();
            let mut pressure = [0.0; COUNT];
            let mut heat = 0.0;
            for sample in &mut pressure[..COUNT / 2] {
                renderer.block(std::slice::from_mut(sample)).unwrap();
                heat += renderer.last_frame().unwrap().routed_material_heat_j;
            }
            let (network, routes, rates) = renderer.into_thermal_with_material_rates();
            assert_eq!(
                routes[0].source.state_coordinate("T"),
                Some(network.temperature_k(CellId(30)).unwrap())
            );
            let mut renderer = direct
                .performance
                .with_material_contact_thermal_network(
                    network,
                    vec![],
                    vec![],
                    routes,
                    profile.controls,
                    gate,
                    |endpoint| prepared.resolve(endpoint),
                )
                .unwrap()
                .with_material_relaxation_temperature(rates)
                .unwrap();
            for sample in &mut pressure[COUNT / 2..] {
                renderer.block(std::slice::from_mut(sample)).unwrap();
                heat += renderer.last_frame().unwrap().routed_material_heat_j;
                assert!(
                    renderer.last_frame().unwrap().total_energy_residual_j.abs()
                        <= profile.controls.total_energy_tolerance_j
                );
            }
            let temperature = renderer.thermal().temperature_k(CellId(30)).unwrap();
            let cell = &renderer.thermal().cells()[&CellId(30)];
            assert!(heat > 0.0 && temperature > INITIAL);
            assert!((cell.mass_kg * cell.specific_h_j_kg - heat).abs() < 4e-12);
            let branches = renderer
                .performance()
                .bank()
                .system()
                .material_relaxation_branches();
            for (index, branch) in branches.iter().enumerate() {
                if !(12..16).contains(&index) && index != 19 {
                    assert_eq!(branch.relaxation_time_s, original_times[index]);
                }
            }
            let time = branches[12].relaxation_time_s;
            assert!(time < original_times[12]);
            let analytical_time = 0.0008 + (temperature - 280.0) * (0.0002 - 0.0008) / 40.0;
            assert!((time - analytical_time).abs() < 1e-16);
            for branch in &branches[12..16] {
                assert_eq!(branch.relaxation_time_s, time);
            }
            let shear_time = branches[19].relaxation_time_s;
            assert!(shear_time > original_times[19]);
            assert!((shear_time - (0.0002 + (temperature - 280.0) * 0.0006 / 40.0)).abs() < 1e-16);
            drop(renderer);
            if outcomes.is_empty() {
                let mut streamed = torsional_ready(cx, gate, &profile, cp, false);
                let prepared = profile.prepare(&streamed, gate).unwrap();
                let mut trace = Vec::new();
                let (wav, summary) = prepared
                    .write(
                        &mut streamed,
                        Cursor::new(Vec::new()),
                        Some(&mut trace),
                        gate,
                    )
                    .unwrap();
                assert_eq!(
                    wav.into_inner(),
                    fs_couple::pcm_wav::encode_pcm16_wav(&pressure, 48_000, 0.01)
                        .unwrap()
                        .0
                );
                assert_eq!(summary.samples, COUNT as u64);
                assert_eq!(String::from_utf8(trace).unwrap().lines().count(), COUNT + 1);
            }
            outcomes.push((pressure, temperature, time));
        }
        assert!(outcomes[0].1 > outcomes[1].1);
        assert!(outcomes[0].2 < outcomes[1].2);
        assert_eq!(outcomes[0].0[0], outcomes[1].0[0]);
        assert_ne!(outcomes[0].0, outcomes[1].0);
    });
}

#[test]
fn v4_string_heat_deposits_every_native_bending_and_torsional_loss_once() {
    scoped(|cx, gate| {
        let core_volume = core::f64::consts::PI * CORE * CORE * 0.5;
        let one = text(true);
        let begin = one.find("cells 1\n").unwrap();
        let end = one.find("dry-contacts").unwrap();
        let mut profile_text = one.clone();
        let mut cells = "cells 4\n".to_owned();
        for (id, pack, volume) in [
            (10, "relaxing-wire.fsmatpk", core_volume),
            (11, "sleeve.fsmatpk", volume()),
            (20, "relaxing-wire.fsmatpk", core_volume),
            (30, "sleeve.fsmatpk", volume()),
        ] {
            cells +=
                &format!("cell {id} {pack} 280 {INITIAL} 320 0 0 0 0.5 0 0 {volume:.17e} 0.001\n");
        }
        profile_text.replace_range(begin..end, &cells);
        profile_text = profile_text.replace("string-heat 1\nstring-body 1 1 30 source-temperature",
            "string-heat 4\nstring-body 0 0 10 source-temperature\nstring-body 0 1 11 source-temperature\nstring-body 1 0 20 source-temperature\nstring-body 1 1 30 source-temperature");
        let profile = Profile::parse(profile_text.as_bytes()).unwrap();
        let mut ready = torsional_ready(cx, gate, &profile, 0.001, false);
        let prepared = profile.prepare(&ready, gate).unwrap();
        let network = prepared.network(gate).unwrap();
        let (routes, rates) = prepared.material_routes(&ready, &network, gate).unwrap();
        let mut renderer = ready
            .performance
            .with_material_contact_thermal_network(
                network,
                vec![],
                vec![],
                routes,
                profile.controls,
                gate,
                |endpoint| prepared.resolve(endpoint),
            )
            .unwrap()
            .with_material_relaxation_temperature(rates)
            .unwrap();
        for _ in 0..3 {
            renderer.block(&mut [0.0]).unwrap();
            let original = renderer
                .performance()
                .bank()
                .system()
                .coupled_frame()
                .unwrap();
            let heat = renderer.last_frame().unwrap().routed_material_heat_j;
            assert!(heat > 0.0);
            assert!(
                (heat - original.material_dissipation_j).abs() < 4e-12,
                "all annular cells must receive exactly the original complete material loss"
            );
        }
    });
}

#[test]
fn v4_string_heat_refuses_incomplete_material_routes_before_advancing() {
    for (from, to) in [
        (
            "string-heat 1\nstring-body 1 1 30",
            "string-heat 2\nstring-body 1 1 30\nstring-body 1 1 30",
        ),
        ("string-body 1 1 30", "string-body 1 1 31"),
        (
            "string-body 1 1 30",
            "string-body 1 1 30 fitted-temperature",
        ),
        ("fiber-heat 0", "fiber-heat 1\nfiber-body 0 0 30"),
        ("films 0", "films 1\nfilm 0 30"),
    ] {
        assert!(Profile::parse(text(false).replace(from, to).as_bytes()).is_err());
    }
    scoped(|cx, gate| {
        let fixed = Profile::parse(text(false).as_bytes()).unwrap();
        let original = ready(cx, gate, &fixed, 0.001, true, false);
        assert!(fixed.prepare(&original, gate).is_err());
        assert_eq!(original.performance.samples_rendered(), 0);
        for (from, to) in [
            ("string-body 1 1 30", "string-body 2 1 30"),
            ("string-body 1 1 30", "string-body 1 2 30"),
            ("cell 30 sleeve.fsmatpk", "cell 30 relaxing-wire.fsmatpk"),
        ] {
            let profile = Profile::parse(text(true).replace(from, to).as_bytes()).unwrap();
            let original = ready(cx, gate, &profile, 0.001, true, false);
            assert!(profile.prepare(&original, gate).is_err());
            assert_eq!(original.performance.samples_rendered(), 0);
        }
        let mut profile = Profile::parse(text(true).as_bytes()).unwrap();
        profile.cells[0].geometry.volume_m3 *= 2.0;
        let original = ready(cx, gate, &profile, 0.001, true, false);
        assert!(profile.prepare(&original, gate).is_err());
        let profile = Profile::parse(text(false).as_bytes()).unwrap();
        let mut twist = spec(&profile, false);
        twist.strings[1].twist_modes = 1;
        let original = twist
            .compile(cx, gate, |name| Ok(pack(name, 20.0, false)))
            .unwrap();
        assert!(profile.prepare(&original, gate).is_err());

        let profile = Profile::parse(text(true).as_bytes()).unwrap();
        let original = torsional_ready(cx, gate, &profile, 0.001, false);
        let prepared = profile.prepare(&original, gate).unwrap();
        let network = prepared.network(gate).unwrap();
        let (routes, rates) = prepared.material_routes(&original, &network, gate).unwrap();
        let system = original.performance.bank().system();
        let before = system.material_memory().to_vec();
        for omission in [0, 4] {
            let mut incomplete = rates.clone();
            incomplete[0].times.remove(omission);
            assert!(
                system
                    .admit_material_rate_thermal_feedback(&network, &routes, &incomplete, gate)
                    .is_err()
            );
        }
        let mut unrelated = rates.clone();
        unrelated[0].times[4].0 = 18;
        assert!(
            system
                .admit_material_rate_thermal_feedback(&network, &routes, &unrelated, gate)
                .is_err()
        );
        assert_eq!(system.material_memory(), before);
        assert_eq!(original.performance.samples_rendered(), 0);
    });
}
