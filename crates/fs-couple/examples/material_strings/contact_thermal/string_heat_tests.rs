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
    viscous_ready(cx, gate, profile, cp, explicit, [0.0; 2], false, false)
}

#[allow(clippy::too_many_arguments)]
fn viscous_ready(
    cx: &Cx<'_>,
    gate: &CancelGate,
    profile: &Profile,
    cp: f64,
    explicit: bool,
    viscosity: [f64; 2],
    only_kelvin: bool,
    varying_viscosity: bool,
) -> Ready {
    viscous_ready_with_motion(cx, gate, profile, cp, explicit, viscosity,
        only_kelvin, varying_viscosity, false, false)
}

#[allow(clippy::too_many_arguments)]
fn viscous_ready_with_motion(
    cx: &Cx<'_>, gate: &CancelGate, profile: &Profile, cp: f64, explicit: bool,
    viscosity: [f64; 2], only_kelvin: bool, varying_viscosity: bool, moving: bool,
    solved_normal: bool,
) -> Ready {
    let mut specification = spec(profile, explicit);
    if moving {
        let motion = "bow_station_motion 2 0.00003 1048576 2\n\
            bow_station_knot 6 0.20012\nbow_station_knot 12 0.2\n";
        specification.bow_station_motion = input::Spec::parse(
            include_str!("../relaxing-layered-strings.material-strings")
                .replace("forces preload", &format!("{motion}forces preload")).as_bytes())
            .unwrap().bow_station_motion;
    }
    if moving || solved_normal {
        specification.schedule = GestureSchedule::try_new(48_000, vec![GestureTrack {
            id: "bow".into(), target: GestureTarget::BowStroke { string: 1 },
            initial: GestureValue::Bow { velocity_m_per_s: 0.01,
                normal_force_n: if solved_normal { 0.0 } else { 0.001 }, station: 0.2 },
            events: vec![],
        }]).unwrap();
    }
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
        if only_kelvin {
            wire.bending_loss = input::BendingLoss::KelvinVoigt;
            wire.torsion_memory = None;
        }
    }
    if only_kelvin {
        specification.bending_memory = None;
    }
    let normal = solved_normal.then(|| crate::stops::Profile::parse(concat!(
        "frankensim-string-stops-v1\n",
        "stop bow_surface 1 0.2 normal 0 -1 fixed -0.000001 5000 1 0 10 0.001\n",
        "bow-contact 0 bow_surface\n",
        "motion 0.001 0.1 3\n",
        "knot 3 -0.0000001\nknot 6 -0.0000002\nknot 12 0\n",
    ).as_bytes()).unwrap());
    specification
        .compile_with_profiles(cx, gate, None, normal, |name| {
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
                if ![ETA, SHEAR_ETA].contains(&property.key.name()) {
                    claims.insert_claim(property.clone()).unwrap();
                }
            }
            for (key, dims, value) in [
                (ETA, DynViscosity::DIMS, viscosity[0]),
                (SHEAR_ETA, DynViscosity::DIMS, viscosity[1]),
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
                if varying_viscosity && [ETA, SHEAR_ETA].contains(&key) {
                    property.value = PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: if key == ETA {
                            vec![(280.0, value), (320.0, value * 2.0)]
                        } else {
                            vec![(280.0, value * 2.0), (320.0, value)]
                        },
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
        // The original torsional viscosity now has its own native heat rows.
        let prepared = profile.prepare(&original, gate).unwrap();
        let network = prepared.network(gate).unwrap();
        let (routes, _) = prepared.material_routes(&original, &network, gate).unwrap();
        assert_eq!(
            routes[0]
                .channel
                .as_ref()
                .unwrap()
                .viscosity_sources()
                .len(),
            1
        );

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

#[test]
fn v4_string_heat_uses_source_viscosities_with_one_layer_inventory_and_native_history() {
    scoped(|cx, gate| {
        for only_kelvin in [false, true] {
            let profile = Profile::parse(text(!only_kelvin).as_bytes()).unwrap();
            let mut observations = Vec::new();
            for scale in [1.0, 3.0] {
                let viscosity = [2e9 * scale, 2e5 * scale];
                let mut ready = viscous_ready(
                    cx,
                    gate,
                    &profile,
                    0.001,
                    true,
                    viscosity,
                    only_kelvin,
                    false,
                );
                let prepared = profile.prepare(&ready, gate).unwrap();
                let network = prepared.network(gate).unwrap();
                let mass = network.cells()[&CellId(30)].mass_kg;
                let (routes, rates) = prepared.material_routes(&ready, &network, gate).unwrap();
                let channel = routes[0].channel.as_ref().unwrap();
                assert_eq!(channel.groups().is_empty(), only_kelvin);
                let [source] = channel.viscosity_sources() else {
                    panic!("one original annulus");
                };
                assert_eq!((source.component(), source.region()), (1, 1));
                // Independent annular geometry and native mass-normalized
                // velocity powers, including both planes and the twist mode.
                let inertia = core::f64::consts::PI * (OUTER.powi(4) - CORE.powi(4)) / 4.0;
                let mu = 7800.0 * core::f64::consts::PI * OUTER.powi(2);
                let rotary = 7800.0 * core::f64::consts::PI * OUTER.powi(4) / 2.0;
                let k = core::f64::consts::PI / 0.5;
                let first = viscosity[0] * inertia * k.powi(4) / mu;
                let twist = viscosity[1] * 2.0 * inertia * k.powi(2) / rotary;
                let expected = [
                    (5, first),
                    (6, first * 16.0),
                    (7, first),
                    (8, first * 16.0),
                    (9, twist),
                ];
                assert_eq!(source.rates_per_s().len(), expected.len());
                for (&(mode, rate), (index, value)) in source.rates_per_s().iter().zip(expected) {
                    assert_eq!(mode, index);
                    assert!((rate - value).abs() <= 2e-13 * value);
                }
                let memory = ready.performance.bank().system().material_memory().len();
                assert_eq!(memory, if only_kelvin { 0 } else { 20 });
                let renderer = ready
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
                    .unwrap();
                let mut renderer = if only_kelvin {
                    renderer
                } else {
                    renderer
                        .with_material_relaxation_temperature(rates)
                        .unwrap()
                };
                let mut pressure = [0.0; 6];
                let mut heat = 0.0;
                for sample in &mut pressure {
                    renderer.block(std::slice::from_mut(sample)).unwrap();
                    let frame = renderer.last_frame().unwrap();
                    heat += frame.routed_material_heat_j;
                    assert!(
                        frame.total_energy_residual_j.abs()
                            <= profile.controls.total_energy_tolerance_j
                    );
                }
                let temperature = renderer.thermal().temperature_k(CellId(30)).unwrap();
                let cell = &renderer.thermal().cells()[&CellId(30)];
                assert!(heat > 0.0 && temperature > INITIAL);
                assert!((cell.mass_kg * cell.specific_h_j_kg - heat).abs() < 4e-12);
                assert_eq!(
                    renderer
                        .performance()
                        .bank()
                        .system()
                        .material_memory()
                        .len(),
                    memory
                );
                observations.push((pressure, temperature, mass));
            }
            assert_eq!(observations[0].2, observations[1].2);
            assert_ne!(observations[0].0, observations[1].0);
            assert_ne!(observations[0].1, observations[1].1);
            let fixed = Profile::parse(text(false).as_bytes()).unwrap();
            let varying = viscous_ready(
                cx,
                gate,
                &fixed,
                0.001,
                true,
                [2e9, 2e5],
                only_kelvin,
                true,
            );
            assert!(
                fixed.prepare(&varying, gate).is_err(),
                "temperature-dependent viscosity requires explicit source-temperature continuation"
            );
            assert_eq!(varying.performance.samples_rendered(), 0);
        }
    });
}

#[test]
fn v4_source_viscosities_and_times_continue_a_moving_bow_and_warm_annulus() {
    scoped(|cx, gate| {
        let profile = Profile::parse(text(true).as_bytes()).unwrap();
        for only_kelvin in [false, true] {
            let mut outcomes = Vec::new();
            for cp in [0.001, 0.004] {
                let mut ready = viscous_ready_with_motion(cx, gate, &profile, cp, true,
                    [2e9, 2e5], only_kelvin, true, true, false);
                assert!(ready.spec.material_body_thermal && !ready.spec.dry_contact_thermal);
                let initial_map = ready.performance.bank().system().port(0).unwrap().left.clone();
                let original_modes: Vec<_> = ready.performance.bank().system().components()
                    .iter().flat_map(|m| m.modes()).copied().collect();
                let original_times: Vec<_> = ready.performance.bank().system().material_relaxation_branches()
                    .iter().map(|b| b.relaxation_time_s).collect();
                let prepared = profile.prepare(&ready, gate).unwrap();
                let network = prepared.network(gate).unwrap();
                let (routes, rates) = prepared.material_routes(&ready, &network, gate).unwrap();
                assert_eq!(rates[0].times.is_empty(), only_kelvin);
                let original_rows = routes[0].channel.as_ref().unwrap().viscosity_sources()[0].rates_per_s().to_vec();
                let mass = network.cells()[&CellId(30)].mass_kg;
                let mut renderer = ready.performance.with_material_contact_thermal_network(
                    network, vec![], vec![], routes, profile.controls, gate, |endpoint| prepared.resolve(endpoint))
                    .unwrap().with_material_relaxation_temperature(rates).unwrap();
                let mut pressure = [0.0; 6];
                let mut heat = 0.0;
                for sample in &mut pressure[..3] {
                    renderer.block(std::slice::from_mut(sample)).unwrap();
                    heat += renderer.last_frame().unwrap().routed_material_heat_j;
                }
                let (network, routes, rates) = renderer.into_thermal_with_material_rates();
                let temperature = network.temperature_k(CellId(30)).unwrap();
                assert!(temperature > INITIAL);
                assert_eq!(routes[0].source.state_coordinate("T"), Some(temperature));
                let eta_b = routes[0].source.property(ETA).unwrap().value_si();
                let eta_s = routes[0].source.property(SHEAR_ETA).unwrap().value_si();
                assert!((eta_b / 2e9 - (1.0 + (temperature - 280.0) / 40.0)).abs() < 1e-14);
                assert!((eta_s / 2e5 - (2.0 - (temperature - 280.0) / 40.0)).abs() < 1e-14);
                let rows = routes[0].channel.as_ref().unwrap().viscosity_sources()[0].rates_per_s();
                assert_eq!(rows.iter().map(|r| r.0).collect::<Vec<_>>(), vec![5, 6, 7, 8, 9]);
                let original_b = 2e9 * (1.0 + (INITIAL - 280.0) / 40.0);
                let original_s = 2e5 * (2.0 - (INITIAL - 280.0) / 40.0);
                let system = ready.performance.bank().system();
                let map = system.port(0).unwrap().left.clone();
                assert_ne!(map, initial_map);
                let stations = ready.performance.bank().bow_station_fractions().unwrap().to_vec();
                assert!((stations[0] - 0.20006).abs() < 1e-15);
                for (index, (mode, original)) in system.components().iter().flat_map(|m| m.modes())
                    .zip(&original_modes).enumerate() {
                    assert_eq!(mode.angular_frequency_rad_s, original.angular_frequency_rad_s);
                    assert_eq!(mode.pressure_per_modal_velocity, original.pressure_per_modal_velocity);
                    if let Some(at) = rows.iter().position(|row| row.0 == index) {
                        let ratio = if index == 9 { eta_s / original_s } else { eta_b / original_b };
                        assert!((rows[at].1 - original_rows[at].1 * ratio).abs() < 2e-13 * rows[at].1);
                        // The other annulus and any gas resistance remain in
                        // the original native coordinate with their old rates.
                        let before = 2.0 * original.damping_ratio * original.angular_frequency_rad_s - original_rows[at].1;
                        let after = 2.0 * mode.damping_ratio * mode.angular_frequency_rad_s - rows[at].1;
                        assert!((before-after).abs() < 2e-12 * (1.0+before.abs()));
                    } else { assert_eq!(mode, original); }
                }
                for (index, branch) in system.material_relaxation_branches().iter().enumerate() {
                    let expected = rates[0].times.iter().find(|(i, _)| *i == index)
                        .map_or(original_times[index], |(_, name)| routes[0].source.property(name).unwrap().value_si());
                    assert_eq!(branch.relaxation_time_s, expected);
                }
                let state: Vec<_> = system.components().iter().flat_map(|m| m.states()).copied().collect();
                let memory = system.material_memory().to_vec();
                let clock = ready.performance.samples_rendered();
                let mut foreign = rates.clone();
                if only_kelvin { foreign[0].times.push((0, "invented-time".into())); }
                else { foreign[0].times[0].1 = ETA.into(); }
                assert!(system.admit_material_rate_thermal_feedback(&network, &routes, &foreign, gate).is_err());
                let refuse = std::cell::Cell::new(false);
                let refusal_observed = std::cell::Cell::new(false);
                let mut renderer = ready.performance.with_material_contact_thermal_network(
                    network, vec![], vec![], routes, profile.controls, gate, |endpoint| {
                        if refuse.get() {
                            refusal_observed.set(true);
                            return Err("late moving-body thermal source refusal".into());
                        }
                        prepared.resolve(endpoint)
                    })
                    .unwrap().with_material_relaxation_temperature(rates).unwrap();
                assert_eq!(renderer.performance().samples_rendered(), clock);
                assert_eq!(renderer.performance().bank().system().material_memory(), memory);
                assert_eq!(renderer.performance().bank().system().components().iter().flat_map(|m| m.states()).copied().collect::<Vec<_>>(), state);
                for sample in &mut pressure[3..] {
                    renderer.block(std::slice::from_mut(sample)).unwrap();
                    let frame = renderer.last_frame().unwrap();
                    heat += frame.routed_material_heat_j;
                    assert!(frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j);
                }
                let cell = &renderer.thermal().cells()[&CellId(30)];
                assert_eq!(cell.mass_kg, mass);
                assert!((cell.mass_kg * cell.specific_h_j_kg - heat).abs() < 4e-12);
                assert_eq!(renderer.performance().bank().system().material_memory().len(), if only_kelvin { 0 } else { 20 });
                assert_eq!(renderer.performance().samples_rendered(), 6);
                assert_eq!(renderer.performance().bank().bow_station_fractions().unwrap(), &[0.20012]);
                // Finish the successful warm continuation before the injected
                // refusal: a failed block deliberately poisons this renderer.
                let system = renderer.performance().bank().system();
                let state: Vec<_> = system.components().iter().flat_map(|m| m.states()).copied().collect();
                let memory = system.material_memory().to_vec();
                let map = system.port(0).unwrap().left.clone();
                let stations = renderer.performance().bank().bow_station_fractions().unwrap().to_vec();
                let clock = renderer.performance().samples_rendered();
                let modes: Vec<_> = renderer.performance().bank().system().components().iter().flat_map(|m| m.modes()).copied().collect();
                let times: Vec<_> = renderer.performance().bank().system().material_relaxation_branches().iter().map(|b| b.relaxation_time_s).collect();
                let thermal_time = renderer.thermal().time_s();
                let enthalpy = renderer.thermal().cells()[&CellId(30)].specific_h_j_kg;
                refuse.set(true);
                let mut untouched = [123.0];
                assert!(renderer.block(&mut untouched).is_err());
                assert!(refusal_observed.get(), "the failure must reach the staged thermal source");
                assert_eq!(untouched, [123.0]);
                assert_eq!(renderer.performance().samples_rendered(), clock);
                assert_eq!(renderer.thermal().time_s(), thermal_time);
                assert_eq!(renderer.thermal().cells()[&CellId(30)].specific_h_j_kg, enthalpy);
                assert_eq!(renderer.performance().bank().bow_station_fractions().unwrap(), stations);
                assert_eq!(renderer.performance().bank().system().port(0).unwrap().left, map);
                assert_eq!(renderer.performance().bank().system().material_memory(), memory);
                assert_eq!(renderer.performance().bank().system().components().iter().flat_map(|m| m.states()).copied().collect::<Vec<_>>(), state);
                assert_eq!(renderer.performance().bank().system().components().iter().flat_map(|m| m.modes()).copied().collect::<Vec<_>>(), modes);
                assert_eq!(renderer.performance().bank().system().material_relaxation_branches().iter().map(|b| b.relaxation_time_s).collect::<Vec<_>>(), times);
                outcomes.push((pressure, renderer.thermal().temperature_k(CellId(30)).unwrap()));
            }
            assert!(outcomes[0].1 > outcomes[1].1);
            assert_eq!(outcomes[0].0[0], outcomes[1].0[0]);
            assert_ne!(outcomes[0].0, outcomes[1].0);
        }
    });
}

#[test]
fn v4_source_mobilities_follow_solved_normal_motion_and_keep_one_energy_clock() {
    scoped(|cx, gate| {
        let profile = Profile::parse(text(true).as_bytes()).unwrap();
        for moving in [false, true] {
            // The same parsed annular source selects eta_b(T), eta_s(T) and
            // disjoint bending/torsional clocks. The stop profile supplies an
            // authored conservative normal spring, never a measured bow law.
            let mut ready = viscous_ready_with_motion(cx, gate, &profile, 0.001, true,
                [2e9, 2e5], false, true, moving, true);
            let system = ready.performance.bank().system();
            assert_eq!(system.contact_normal_load(0), Some(0));
            let initial_normal = system.obstacle_normal_column(0).unwrap().to_vec();
            let original_modes: Vec<_> = system.components().iter().flat_map(|m| m.modes()).copied().collect();
            let mut balance = system.total_energy_j().unwrap();
            let prepared = profile.prepare(&ready, gate).unwrap();
            let mut network = prepared.network(gate).unwrap();
            let (mut routes, mut rates) = prepared.material_routes(&ready, &network, gate).unwrap();
            let mut received = 0.0;
            let mut normal_work = 0.0;
            let mut loaded = false;
            // A declared 20 mW heater changes the finite source temperature
            // enough to expose stale normal responses above root tolerances.
            // Its external work stays separate from native dissipative heat.
            let heater = |endpoint: &DepositionSimulation<'_>| -> Result<StepCoefficients<()>> {
                let mut coefficients = prepared.resolve(endpoint)?;
                coefficients.sources_w.insert(CellId(30), 0.02);
                Ok(coefficients)
            };
            // Recover and rebind the actual finite bodies at a warm sample
            // boundary, keeping the original modal, memory and contact state.
            for _ in 0..2 {
                let system = ready.performance.bank().system();
                let warm: Vec<_> = system.components().iter().flat_map(|m| m.states()).copied().collect();
                let memory = system.material_memory().to_vec();
                let normal = system.last_obstacle_frames().map(|frames| frames.to_vec());
                let mut renderer = ready.performance.with_material_contact_thermal_network(
                    network, vec![], vec![], routes, profile.controls, gate,
                    heater)
                    .unwrap().with_material_relaxation_temperature(rates).unwrap();
                assert_eq!(renderer.performance().bank().system().components().iter()
                    .flat_map(|m| m.states()).copied().collect::<Vec<_>>(), warm);
                assert_eq!(renderer.performance().bank().system().material_memory(), memory);
                assert_eq!(renderer.performance().bank().system().last_obstacle_frames(), normal.as_deref());
                for _ in 0..3 {
                    let system = renderer.performance().bank().system();
                    let before: Vec<_> = system.components().iter().flat_map(|m| m.states())
                        .map(|s| s.displacement_m_sqrt_kg).collect();
                    let column = system.obstacle_normal_column(0).unwrap().to_vec();
                    let offset = system.obstacle_offsets_m().unwrap()[0];
                    let old_x = column.iter().zip(&before).map(|(n, q)| n*q).sum::<f64>() - offset;
                    renderer.block(&mut [0.0]).unwrap();
                    let system = renderer.performance().bank().system();
                    let after: Vec<_> = system.components().iter().flat_map(|m| m.states())
                        .map(|s| s.displacement_m_sqrt_kg).collect();
                    let end = system.obstacle_normal_column(0).unwrap();
                    let contact = system.last_obstacle_frames().unwrap()[0];
                    let new_x = end.iter().zip(&after).map(|(n, q)| n*q).sum::<f64>() - contact.support_offset_m;
                    let (p0, p1) = ((old_x + 1e-6).max(0.0), (new_x + 1e-6).max(0.0));
                    // Independent secant of U = k max(x-gap,0)^2 / 2.
                    let reaction = if p0 > 0.0 && p1 > 0.0 { 2500.0 * (p0+p1) }
                        else if new_x == old_x { 5000.0 * p0 }
                        else { 2500.0 * (p1*p1-p0*p0) / (new_x-old_x) };
                    assert!((contact.normal_force_n - reaction).abs() <= contact.force_tolerance_n);
                    let transport = end.iter().zip(&column).zip(before.iter().zip(&after))
                        .map(|((b, a), (q0, q1))| 0.5*(b-a)*(q0+q1)).sum::<f64>();
                    let work = contact.normal_force_n * (transport - (contact.support_offset_m-offset));
                    assert!((contact.actuator_work_j-work).abs() < 1e-13);
                    let frame = renderer.last_frame().unwrap();
                    let mechanical = &frame.mechanics;
                    assert_eq!(mechanical.ports[0].normal_force_n, contact.normal_force_n);
                    loaded |= contact.normal_force_n > 0.0 && mechanical.ports[0].traction_n != 0.0;
                    normal_work += contact.actuator_work_j.abs();
                    received += frame.routed_material_heat_j;
                    let external_heat = 0.02 * system.sample_period_s();
                    assert!((frame.external_heat_j-external_heat).abs() < 1e-13);
                    balance += mechanical.external_work_j + mechanical.drive_work_j + work + external_heat
                        - mechanical.network_dissipation_j + frame.routed_material_heat_j
                        - mechanical.friction_dissipation_j - mechanical.declared_slip_work_j
                        - mechanical.obstacle_dissipation_j - mechanical.obstacle_residual_loss_j
                        - mechanical.friction_residual_loss_j;
                    let cell = &renderer.thermal().cells()[&CellId(30)];
                    let heat = cell.mass_kg * cell.specific_h_j_kg;
                    assert!((heat-received-0.02*renderer.thermal().time_s()).abs() < 4e-12);
                    assert!((system.total_energy_j().unwrap()+heat-balance).abs() < 3e-11);
                    assert!(frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j);
                }
                (network, routes, rates) = renderer.into_thermal_with_material_rates();
            }
            assert!(loaded && normal_work > 0.0 && received > 0.0);
            let temperature = network.temperature_k(CellId(30)).unwrap();
            assert!(temperature > INITIAL);
            assert_eq!(routes[0].source.state_coordinate("T"), Some(temperature));
            assert!((routes[0].source.property(ETA).unwrap().value_si()/2e9
                - (1.0 + (temperature-280.0)/40.0)).abs() < 1e-14);
            let system = ready.performance.bank().system();
            assert!(system.components().iter().flat_map(|m| m.modes()).zip(&original_modes)
                .any(|(a, b)| a.damping_ratio != b.damping_ratio));
            for (index, name) in &rates[0].times {
                assert_eq!(system.material_relaxation_branches()[*index].relaxation_time_s,
                    routes[0].source.property(name).unwrap().value_si());
            }
            let normal_map = system.obstacle_normal_column(0).unwrap().to_vec();
            assert_eq!(normal_map != initial_normal, moving);
            let tangent = system.port(0).unwrap().left.clone();
            let normal = system.last_obstacle_frames().unwrap().to_vec();
            let offsets = system.obstacle_offsets_m().unwrap().to_vec();
            let state: Vec<_> = system.components().iter().flat_map(|m| m.states()).copied().collect();
            let modes: Vec<_> = system.components().iter().flat_map(|m| m.modes()).copied().collect();
            let memory = system.material_memory().to_vec();
            let times: Vec<_> = system.material_relaxation_branches().iter().map(|b| b.relaxation_time_s).collect();
            let source = routes[0].source.clone();
            let clock = network.time_s();
            let enthalpy = network.cells()[&CellId(30)].specific_h_j_kg;
            let observed = std::cell::Cell::new(false);
            let mut renderer = ready.performance.with_material_contact_thermal_network(
                network, vec![], vec![], routes, profile.controls, gate, |endpoint| {
                    observed.set(true);
                    if endpoint.time_s() >= clock { return Err("late solved-normal thermal refusal".into()); }
                    prepared.resolve(endpoint)
                }).unwrap().with_material_relaxation_temperature(rates).unwrap();
            let mut untouched = [123.0];
            assert!(renderer.block(&mut untouched).is_err());
            assert!(observed.get());
            assert_eq!(untouched, [123.0]);
            let system = renderer.performance().bank().system();
            assert_eq!(renderer.performance().samples_rendered(), 6);
            assert_eq!(system.obstacle_normal_column(0).unwrap(), normal_map);
            assert_eq!(system.port(0).unwrap().left, tangent);
            assert_eq!(system.last_obstacle_frames().unwrap(), normal);
            assert_eq!(system.obstacle_offsets_m().unwrap(), offsets);
            assert_eq!(system.components().iter().flat_map(|m| m.states()).copied().collect::<Vec<_>>(), state);
            assert_eq!(system.components().iter().flat_map(|m| m.modes()).copied().collect::<Vec<_>>(), modes);
            assert_eq!(system.material_memory(), memory);
            assert_eq!(system.material_relaxation_branches().iter().map(|b| b.relaxation_time_s).collect::<Vec<_>>(), times);
            let (network, routes, _) = renderer.into_thermal_with_material_rates();
            assert_eq!(network.time_s(), clock);
            assert_eq!(network.cells()[&CellId(30)].specific_h_j_kg, enthalpy);
            assert_eq!(routes[0].source, source);
        }
    });
}
