//! Manufactured native input/heat/PCM integration, not measured bow or rosin data.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, PropertyClaim, PropertyValue, Provenance, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_scenario::gesture::{GestureTrack, GestureValue};
use std::io::Cursor;

#[path = "film_surface_tests.rs"]
mod film_surface;

#[path = "material_heat_tests.rs"]
mod material_heat;

#[path = "string_heat_tests.rs"]
mod string_heat;

#[path = "caloric_reference_tests.rs"]
mod caloric_reference;

fn scoped(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 0,
                    kernel_id: 0,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn claim(name: &str, dims: Dims, value: f64, validity: ValidityDomain) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "MANUFACTURED finite-cell dry-heat integration".into(),
            license: "CC0-1.0".into(),
            artifact: Some(fs_blake3::ContentHash([47; 32])),
        },
    }
}

/// The raw observation every compiled fixture pack must retain: a manufactured record, not measured data.
fn observed(claims: &mut ClaimSet, specimen: &str, artifact: [u8; 32]) -> fs_matdb::ObservationId {
    let artifact = fs_blake3::ContentHash(artifact);
    claims
        .register_observation(fs_matdb::ObservationDataset {
            specimen: specimen.into(),
            method: "MANUFACTURED finite-cell dry-heat integration".into(),
            artifact,
            caveats: "manufactured fixture; not measured bow, rosin, hair or string data".into(),
            provenance: Provenance {
                source: "MANUFACTURED finite-cell dry-heat integration".into(),
                license: "CC0-1.0".into(),
                artifact: Some(artifact),
            },
        })
        .unwrap()
}

fn material(name: &str, cp: f64, omit: Option<&str>) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let observation = observed(&mut claims, name, [47; 32]);
    let fiber = name == "fiber.fsmatpk";
    for (key, dims, value) in [
        (
            "density",
            Density::DIMS,
            if fiber { 1300.0 } else { 7800.0 },
        ),
        (
            "young_modulus",
            Pressure::DIMS,
            if fiber { 3e9 } else { 200e9 },
        ),
        ("poisson_ratio", Dims::NONE, 0.3),
        (ETA, DynViscosity::DIMS, 50.0),
        (
            fs_material::section::axial::YOUNG_VISCOSITY_PROPERTY,
            DynViscosity::DIMS,
            50.0,
        ),
        (SHEAR, Pressure::DIMS, 80e9),
        (SHEAR_ETA, DynViscosity::DIMS, 1000.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, cp),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.2),
    ] {
        if Some(key) != omit {
            let mut item = claim(
                key,
                dims,
                value,
                ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, 400_000.0),
            );
            item.observations = vec![observation];
            claims.insert_claim(item).unwrap();
        }
    }
    NormalizedPack::new(
        name,
        "manufactured-v1",
        fs_blake3::ContentHash([47; 32]),
        "manufactured freely redistributable integration data",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

fn interface(receiver: &str) -> NormalizedInterfacePack {
    let mut claims = ClaimSet::new();
    let observation = observed(&mut claims, "manufactured dry contact", [48; 32]);
    for (key, value) in [("mu_static", 0.8), ("mu_kinetic", 0.4)] {
        let mut item = claim(
            key,
            Dims::NONE,
            value,
            ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0)
                .with("p", 0.0, 1e6)
                .with("v", 0.0, 10.0),
        );
        item.observations = vec![observation];
        claims.insert_claim(item).unwrap();
    }
    let claims = NormalizedPack::new(
        "manufactured-dry-contact",
        "manufactured-v1",
        fs_blake3::ContentHash([48; 32]),
        "manufactured CC0 fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    let surface = |name| SurfaceSpec {
        material: MaterialPack::from(material(name, 20.0, None))
            .card()
            .unwrap()
            .id()
            .clone(),
        texture_frame: "manufactured tangent direction".into(),
    };
    NormalizedInterfacePack::new(
        surface("wire.fsmatpk"),
        surface(receiver),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "declared gas".into(),
            history: "manufactured unchanged dry surfaces".into(),
        },
        claims,
    )
    .unwrap()
}

fn profile_text(fiber: bool) -> String {
    let receiver = if fiber {
        "fiber.fsmatpk"
    } else {
        "body.fsmatpk"
    };
    let mapping = if fiber { "fiber -0.238" } else { "fixed" };
    let mut text = String::from(
        "frankensim-contact-thermal-v2\nboundary adiabatic\n\
        solver 512 24 1e-9 1e-12 1e-12 4e-12 1e-13 1e-10 8 1e-8\ncells 8\n",
    );
    for (side, pack) in ["wire.fsmatpk", receiver].into_iter().enumerate() {
        for cell in 0..4 {
            let start = cell as f64 * 0.002;
            let end = (cell + 1) as f64 * 0.002;
            text.push_str(&format!(
                "cell {} {pack} 280 293.15 320 {start} 0 {side} {end} 0 {side} 4e-10 4.048e-5\n",
                10 + side * 10 + cell
            ));
        }
    }
    text.push_str(&format!("dry-contacts 1\ndry-contact 0 -0.098 {mapping}\nfootprint 0.001 0.003 0.01 0.001 0.003 0.01\n"));
    for side in 0..2 {
        text.push_str("surface 4\n");
        for cell in 0..4 {
            text.push_str(&format!(
                "strip {} {} {} 0.01 0.00001\n",
                10 + side * 10 + cell,
                cell as f64 * 0.002,
                (cell + 1) as f64 * 0.002
            ));
        }
    }
    text.push_str("films 0\nseries 0\nlinks 0\n");
    text
}

fn input(fiber: bool, block: usize, rate: u32) -> input::Spec {
    let original = if fiber {
        include_str!("../fiber-ribbon.material-strings")
    } else {
        include_str!("../pluck.material-strings")
    };
    let mut text = original.split_once("\nschedule\n").unwrap().0.to_owned();
    text = text
        .replace(
            "audio 48000 4096 0.01 128",
            &format!("audio {rate} 128 0.01 {block}"),
        )
        .replace("bow_ribbon 0.01 3", "bow_ribbon 0.01 1");
    let law = "friction interfaces 1\nbow-interface bow dry.fsintpk 293.15 0.00002 50 0 coulomb mu_static mu_kinetic";
    text = text
        .replace("friction stribeck 0.8 0.4 0.1 0", law)
        .replace("friction coulomb 0.8 0.4", law);
    let motion = if fiber {
        "bow_fiber_program 0 0.1 100 0.01 1048576 3\nfiber_speed 32 0.01\nfiber_speed 64 -0.01\nfiber_speed 96 0\n"
    } else {
        "bow_station_motion 10 0.0002 1048576 2\nbow_station_knot 32 0.208\nbow_station_knot 64 0.2\n"
    };
    text = text.replace("forces ", &format!("{motion}forces "));
    let schedule = GestureSchedule::try_new(
        48_000,
        vec![GestureTrack {
            id: "bow".into(),
            target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow {
                velocity_m_per_s: if fiber { 0.0 } else { 0.03 },
                normal_force_n: 0.001,
                station: 0.2,
            },
            events: vec![],
        }],
    )
    .unwrap();
    text.push_str("\nschedule\n");
    text.push_str(std::str::from_utf8(&schedule.to_canonical_bytes()).unwrap());
    input::Spec::parse(text.as_bytes()).unwrap()
}

fn ready(
    cx: &Cx<'_>,
    gate: &CancelGate,
    profile: &Profile,
    fiber: bool,
    cp: f64,
    block: usize,
    rate: u32,
    omit: Option<&str>,
) -> Ready {
    let mut source = input(fiber, block, rate);
    profile.configure(&mut source);
    source
        .friction
        .load(gate, |_| {
            Ok(interface(if fiber {
                "fiber.fsmatpk"
            } else {
                "body.fsmatpk"
            }))
        })
        .unwrap();
    source
        .compile(cx, gate, |name| Ok(material(name, cp, omit)))
        .unwrap()
}

fn conduction_profile(record: &str) -> Profile {
    // Two extra finite wire cells share one internal link, separate from the
    // moving dry footprint. Geometry and all source data are manufactured.
    let cells = "cell 30 wire.fsmatpk 280 310 320 0 0 2 0.001 0 2 1e-9 1e-4\n\
        cell 31 wire.fsmatpk 280 290 320 0.001 0 2 0.002 0 2 1e-9 1e-4\n";
    let text = profile_text(false)
        .replace("cells 8", "cells 10")
        .replace("dry-contacts 1", &format!("{cells}dry-contacts 1"))
        .replace("links 0", &format!("links 1\n{record}"));
    Profile::parse(text.as_bytes()).unwrap()
}

const INTEGRATED_LINK: &str = "integrated-link 9 30 31 1e-6 1e-5 1e-5 96 1e-12 1e-12";
const ENDPOINT_LINK: &str = "link 9 30 31 1e-6 1e-5 1e-5";

fn conductivity_pack(name: &str, values: &[(f64, f64)]) -> NormalizedPack {
    let base = material(name, 20.0, Some("thermal_conductivity"));
    let mut claims = base.claims().clone();
    let observation = observed(&mut claims, "manufactured nonlinear conduction", [49; 32]);
    let mut item = claim(
        "thermal_conductivity",
        THERMAL_CONDUCTIVITY_DIMS,
        0.0,
        ValidityDomain::unconstrained()
            .with("T", values[0].0, values.last().unwrap().0)
            .with("omega", 0.0, 400_000.0),
    );
    item.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: fs_qty::Temperature::DIMS,
        knots: values.to_vec(),
        dims: THERMAL_CONDUCTIVITY_DIMS,
    };
    item.interpolation = InterpolationPolicy::LinearInside;
    item.observations = vec![observation];
    claims.insert_claim(item).unwrap();
    NormalizedPack::new(
        name,
        "manufactured-v1",
        fs_blake3::ContentHash([49; 32]),
        "manufactured freely redistributable conduction fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

fn conduction_ready(
    cx: &Cx<'_>,
    gate: &CancelGate,
    profile: &Profile,
    values: &[(f64, f64)],
) -> Ready {
    let mut source = input(false, 37, 48_000);
    profile.configure(&mut source);
    source
        .friction
        .load(gate, |_| Ok(interface("body.fsmatpk")))
        .unwrap();
    source
        .compile(cx, gate, |name| {
            Ok(if name == "wire.fsmatpk" {
                conductivity_pack(name, values)
            } else {
                material(name, 20.0, None)
            })
        })
        .unwrap()
}

#[test]
fn native_integrated_link_uses_source_integrals_in_the_original_enthalpy_step() {
    scoped(|cx, gate| {
        let profile = conduction_profile(INTEGRATED_LINK);
        let values = [(280.0, 0.1), (300.0, 0.2), (320.0, 0.3)];
        let source = conduction_ready(cx, gate, &profile, &values);
        let prepared = profile.prepare(&source, gate).unwrap();
        let original = prepared.network(gate).unwrap();
        let chart = prepared.sources[8]
            .integrated_conductivity
            .as_ref()
            .unwrap();
        assert_eq!(
            chart.curve.knots(),
            &[[280.0, 0.1], [300.0, 0.2], [320.0, 0.3]]
        );
        assert_eq!(chart.source.knots.len(), 3);
        assert!(chart.source.knots.iter().all(|(_, answer)| {
            answer.receipt.selected
                == prepared.sources[8]
                    .state
                    .property("thermal_conductivity")
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
        }));

        // Identical affine k(T), equal capacities: Ta+Tb=600 K, hence the
        // integral over both half-cells gives G=A*k(300)/(da+db)=0.01 W/K.
        // The old endpoint harmonic mean would give only 0.009375 W/K.
        let coefficient = prepared.resolve(&original).unwrap();
        assert!(coefficient.sources_w.is_empty());
        let g = 0.01;
        assert!((coefficient.conductances_w_k[&9] - g).abs() < 1e-12);
        let endpoint = series_conductance_w_k(1e-6, [1e-5; 2], [0.25, 0.15]).unwrap();
        assert!((endpoint - g).abs() > 0.0006);
        let dt = prepared.sample_period_s;
        let capacity = 7800.0 * 1e-9 * 20.0;
        let difference = 20.0 / (1.0 + 2.0 * g * dt / capacity);
        let exchange = capacity * 0.5 * (20.0 - difference);
        let mut increments = vec![0.0; 10];
        increments[8] = -exchange;
        increments[9] = exchange;
        let step =
            fs_thermal::enthalpy::backward_euler::BackwardEulerStep::new(&original, dt, || Ok(()))
                .unwrap();
        let trial = step
            .evaluate(&increments, |state| prepared.resolve(state), || Ok(()))
            .unwrap();
        assert!(trial.maximum_cell_residual_j() < 1e-12);
        let (next, (), audit) = trial.accept(1e-12, 1e-12, || Ok(())).unwrap();
        assert!(
            (next.temperature_k(CellId(30)).unwrap() - (300.0 + 0.5 * difference)).abs() < 1e-9
        );
        assert!(
            (next.temperature_k(CellId(31)).unwrap() - (300.0 - 0.5 * difference)).abs() < 1e-9
        );
        assert_eq!(audit.external_heat_j, 0.0);
        assert!(audit.residual_j.abs() < 1e-12);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(next.time_s(), dt);
        assert_eq!(next.links()[&9], original.links()[&9]);
        assert_eq!(source.performance.samples_rendered(), 0);
        let fresh = prepared.resolve(&next).unwrap().conductances_w_k[&9]
            * (next.temperature_k(CellId(30)).unwrap() - next.temperature_k(CellId(31)).unwrap());
        assert!(fresh < g * 20.0);
        assert!(
            prepared
                .conductance(&profile.links[0], [279.0, 300.0])
                .is_err()
        );
        assert!(
            prepared
                .conductance(&profile.links[0], [300.0, 321.0])
                .is_err()
        );
    });
}

#[test]
fn integrated_links_retain_constant_k_and_legacy_zero_insulation() {
    scoped(|cx, gate| {
        for record in [ENDPOINT_LINK, INTEGRATED_LINK] {
            let profile = conduction_profile(record);
            let source = conduction_ready(cx, gate, &profile, &[(280.0, 0.2), (320.0, 0.2)]);
            let prepared = profile.prepare(&source, gate).unwrap();
            let (g, residual) = prepared
                .conductance(&profile.links[0], [310.0, 290.0])
                .unwrap();
            assert!((g - 0.01).abs() < 1e-14);
            assert!(residual < 1e-12);
            // The independent dry footprint remains conducting at 293.15 K;
            // both extra link endpoints have explicit source k=0 instead.
            let insulating = conduction_ready(
                cx,
                gate,
                &profile,
                &[
                    (280.0, 0.0),
                    (290.0, 0.0),
                    (293.15, 0.2),
                    (300.0, 0.2),
                    (310.0, 0.0),
                    (320.0, 0.0),
                ],
            );
            if record == ENDPOINT_LINK {
                let prepared = profile.prepare(&insulating, gate).unwrap();
                assert_eq!(
                    prepared
                        .conductance(&profile.links[0], [310.0, 290.0])
                        .unwrap(),
                    (0.0, 0.0)
                );
            } else {
                assert!(profile.prepare(&insulating, gate).is_err());
            }
        }
    });
}

#[test]
fn integrated_links_refuse_missing_source_paths_heat_error_and_cancellation() {
    scoped(|cx, gate| {
        let profile = conduction_profile(INTEGRATED_LINK);
        let source = conduction_ready(
            cx,
            gate,
            &profile,
            &[(290.0, 0.15), (300.0, 0.2), (310.0, 0.25)],
        );
        assert!(
            profile.prepare(&source, gate).is_err(),
            "endpoint support is insufficient for a full source path"
        );
        assert_eq!(source.performance.samples_rendered(), 0);

        let loose = conduction_profile("integrated-link 9 30 31 1e-6 1e-5 1e-5 1 1 0");
        let values = [(280.0, 0.1), (300.0, 0.2), (320.0, 0.3)];
        let source = conduction_ready(cx, gate, &loose, &values);
        let error = match loose.prepare(&source, gate) {
            Ok(_) => panic!("loose layer tolerance bypassed the finite-cell heat budget"),
            Err(error) => error,
        };
        assert!(error.to_string().contains("layer residual"), "{error}");
        assert_eq!(source.performance.samples_rendered(), 0);

        let source = conduction_ready(cx, gate, &profile, &values);
        let prepared = profile.prepare(&source, gate).unwrap();
        let network = prepared.network(gate).unwrap();
        gate.request();
        assert!(prepared.resolve(&network).is_err());
        assert_eq!(network.time_s(), 0.0);
        assert_eq!(source.performance.samples_rendered(), 0);
    });
}

#[test]
fn native_dry_heat_follows_moving_and_revisited_material_cells_and_matches_pcm() {
    scoped(|cx, gate| {
        let profile = Profile::parse(profile_text(false).as_bytes()).unwrap();
        let mut direct = ready(cx, gate, &profile, false, 20.0, 37, 48_000, None);
        let prepared = profile.prepare(&direct, gate).unwrap();
        let routes = prepared.dry_routes(&direct).unwrap();
        let network = prepared.network(gate).unwrap();
        let mut coupled = direct
            .performance
            .with_dry_contact_thermal_network(
                network,
                routes,
                vec![],
                profile.controls,
                gate,
                |endpoint| prepared.resolve(endpoint),
            )
            .unwrap();
        let mut expected = vec![0.0; 128];
        coupled.block(&mut expected[..32]).unwrap();
        let reached = coupled.last_frame().unwrap().dry_contacts[0]
            .footprint
            .intervals_m[0];
        assert!((reached[0] - (0.005 - 0.5 * 0.004 / 32.0)).abs() < 1e-14);
        let retained = coupled.thermal().cells()[&CellId(12)].specific_h_j_kg;
        assert!(retained > 0.0);
        coupled.block(&mut expected[32..64]).unwrap();
        let returned = coupled.last_frame().unwrap().dry_contacts[0]
            .footprint
            .intervals_m[0];
        assert!((returned[0] - (0.001 + 0.5 * 0.004 / 32.0)).abs() < 1e-14);
        let abandoned = coupled.thermal().cells()[&CellId(12)].specific_h_j_kg;
        assert!(abandoned > 0.0);
        coupled.block(&mut expected[64..96]).unwrap();
        coupled.block(&mut expected[96..]).unwrap();
        let final_heat = coupled
            .thermal()
            .cells()
            .values()
            .map(|cell| cell.mass_kg * cell.specific_h_j_kg)
            .sum::<f64>();
        assert!(final_heat > 0.0);
        assert_eq!(
            coupled.thermal().cells()[&CellId(12)].specific_h_j_kg,
            abandoned
        );
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        let temperatures: Vec<_> = coupled
            .thermal()
            .cells()
            .keys()
            .map(|&id| coupled.thermal().temperature_k(id).unwrap())
            .collect();
        assert!(temperatures.iter().all(|t| *t >= 293.15));
        let final_frame = coupled.last_frame().unwrap().mechanics.clone();
        drop(coupled);
        let encoded = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48_000, 0.01).unwrap();
        let mut streamed = ready(cx, gate, &profile, false, 20.0, 128, 48_000, None);
        let stream = profile.prepare(&streamed, gate).unwrap();
        let mut trace = Vec::new();
        let (wav, summary) = stream
            .write(
                &mut streamed,
                Cursor::new(Vec::new()),
                Some(&mut trace),
                gate,
            )
            .unwrap();
        assert_eq!(wav.into_inner(), encoded.0);
        assert_eq!(summary.samples, 128);
        assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 129);
        assert_eq!(
            streamed.performance.bank().system().last_frame(),
            Some(&final_frame)
        );
    });
}

#[test]
fn native_fiber_acceleration_and_material_heat_capacity_reach_the_same_dry_owner() {
    scoped(|cx, gate| {
        let profile = Profile::parse(profile_text(true).as_bytes()).unwrap();
        let mut temperatures = Vec::new();
        for cp in [20.0, 2000.0] {
            let mut source = ready(cx, gate, &profile, true, cp, 37, 48_000, None);
            let prepared = profile.prepare(&source, gate).unwrap();
            let routes = prepared.dry_routes(&source).unwrap();
            let mut renderer = source
                .performance
                .with_dry_contact_thermal_network(
                    prepared.network(gate).unwrap(),
                    routes,
                    vec![],
                    profile.controls,
                    gate,
                    |endpoint| prepared.resolve(endpoint),
                )
                .unwrap();
            let mut output = vec![0.0; 128];
            renderer.block(&mut output[..32]).unwrap();
            let moved = renderer.last_frame().unwrap().dry_contacts[0]
                .footprint
                .intervals_m[1];
            assert_ne!(moved, [0.001, 0.003]);
            for chunk in output[32..].chunks_mut(32) {
                renderer.block(chunk).unwrap();
            }
            let heat: f64 = renderer
                .thermal()
                .cells()
                .values()
                .map(|cell| cell.mass_kg * cell.specific_h_j_kg)
                .sum();
            assert!(heat > 0.0);
            temperatures.push(renderer.thermal().temperature_k(CellId(20)).unwrap());
            let frame = renderer.performance().fiber_support_frame().unwrap();
            assert_eq!(frame.velocity_m_s, 0.0);
            assert_eq!(frame.sample, 128);
            assert!(output.iter().any(|p| p.abs() > 1e-12));
            drop(renderer);
            let mut streamed = ready(cx, gate, &profile, true, cp, 128, 48_000, None);
            let stream = profile.prepare(&streamed, gate).unwrap();
            let (wav, _) = stream
                .write(&mut streamed, Cursor::new(Vec::new()), None, gate)
                .unwrap();
            assert_eq!(
                wav.into_inner(),
                fs_couple::pcm_wav::encode_pcm16_wav(&output, 48_000, 0.01)
                    .unwrap()
                    .0
            );
        }
        assert!(temperatures[0] > temperatures[1]);
        assert!(temperatures[1] > 293.15);
    });
}

#[test]
fn native_dry_profiles_refuse_wrong_participants_missing_thermal_data_and_false_geometry() {
    scoped(|cx, gate| {
        for (from, to) in [
            ("cell 20 body.fsmatpk", "cell 20 wire.fsmatpk"),
            ("4e-10 4.048e-5", "8e-10 4.048e-5"),
            ("footprint 0.001 0.003 0.01", "footprint 0.001 0.003 0.02"),
            ("dry-contact 0 -0.098 fixed", "dry-contact 0 -0.097 fixed"),
            (
                "dry-contact 0 -0.098 fixed",
                "dry-contact 0 -0.098 fiber -0.238",
            ),
            ("links 0", "links 1\nlink 1 10 20 0.00002 0.00001 0.00001"),
        ] {
            let profile = Profile::parse(profile_text(false).replace(from, to).as_bytes()).unwrap();
            let source = ready(cx, gate, &profile, false, 20.0, 128, 48_000, None);
            assert!(
                profile.prepare(&source, gate).is_err(),
                "accepted {from} -> {to}"
            );
            assert_eq!(source.performance.samples_rendered(), 0);
        }
        let profile = Profile::parse(profile_text(false).as_bytes()).unwrap();
        for absent in ["specific_heat_capacity", "thermal_conductivity"] {
            let source = ready(cx, gate, &profile, false, 20.0, 128, 48_000, Some(absent));
            assert!(profile.prepare(&source, gate).is_err());
        }
        for (from, to) in [
            ("surface 4", "surface 0"),
            ("strip 11", "strip 10"),
            ("strip 10 0 0.002", "strip 10 0.002 0"),
            ("films 0", "films 1\nfilm 0 10"),
        ] {
            assert!(Profile::parse(profile_text(false).replace(from, to).as_bytes()).is_err());
        }
    });
}

fn temperature_interface(slope: f64) -> NormalizedInterfacePack {
    let original = interface("body.fsmatpk");
    let mut claims = ClaimSet::new();
    // The copied claims keep their observation, so its record comes across first.
    let source = original.card().claims();
    for id in source.observation_ids() {
        claims
            .register_observation(source.observation(id).unwrap().clone())
            .unwrap();
    }
    // The sourced Coulomb owner refuses a query whose nominal pressure is not a typed static pressure, and `interface()`
    // declares T, p and v bare for the non-thermal paths. The same ranges are declared here with the quantities the
    // resolver builds its query in (as tests/bow_interfaces.rs held_environment::source does); T stays bare, which the
    // owner admits as the legacy kelvin axis.
    let typed = || {
        fs_evidence::ValidityDomain::unconstrained()
            .with("T", 280.0, 320.0)
            .with_quantity(
                "p",
                fs_qty::QuantitySpec::semantic(fs_qty::semantic::SemanticType::new(
                    fs_qty::semantic::QuantityKind::Pressure,
                    fs_qty::semantic::ValueForm::Static,
                )),
                0.0,
                1e6,
            )
            .with_quantity(
                "v",
                fs_qty::QuantitySpec::dimensional(fs_qty::Velocity::DIMS),
                0.0,
                10.0,
            )
    };
    for (_, item) in original.card().claims().claims_ordered() {
        let mut item = item.clone();
        item.validity = typed();
        if item.key.name() == "mu_kinetic" {
            item.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: [280.0, 320.0]
                    .map(|t| (t, 0.4 - slope * (t - 293.15)))
                    .to_vec(),
                dims: Dims::NONE,
            };
            item.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(item).unwrap();
    }
    let pack = NormalizedPack::new(
        "manufactured-thermal-dry-contact",
        "manufactured-v1",
        fs_blake3::ContentHash([49; 32]),
        "manufactured CC0 fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedInterfacePack::new(
        original.card().surface_a().clone(),
        original.card().surface_b().clone(),
        original.card().context().clone(),
        pack,
    )
    .unwrap()
}

fn temperature_ready(
    cx: &Cx<'_>,
    gate: &CancelGate,
    profile: &Profile,
    cp: f64,
    slope: f64,
    block: usize,
) -> Ready {
    let mut spec = input(false, block, 48_000);
    profile.configure(&mut spec);
    spec.friction
        .load(gate, |_| Ok(temperature_interface(slope)))
        .unwrap();
    spec.compile(cx, gate, |name| Ok(material(name, cp, None)))
        .unwrap()
}

#[test]
fn native_temperature_source_changes_traction_and_pressure_through_finite_heat() {
    scoped(|cx, gate| {
        let profile = Profile::parse(profile_text(false).as_bytes()).unwrap();
        let mut runs = Vec::new();
        for (cp, slope) in [(0.02, 0.005), (20.0, 0.005), (0.02, 0.0)] {
            let mut source = temperature_ready(cx, gate, &profile, cp, slope, 37);
            let prepared = profile.prepare(&source, gate).unwrap();
            let routes = prepared.dry_routes(&source).unwrap();
            let mut render = source
                .performance
                .with_dry_contact_thermal_network(
                    prepared.network(gate).unwrap(),
                    routes,
                    vec![],
                    profile.controls,
                    gate,
                    |endpoint| prepared.resolve(endpoint),
                )
                .unwrap();
            let mut pressure = vec![0.0; 128];
            for block in pressure.chunks_mut(37) {
                render.block(block).unwrap();
            }
            let frame = render.last_frame().unwrap();
            let temperature = frame.dry_contacts[0].interface_temperature_k.unwrap();
            let binding = render.performance().bank().bow_interface(0).unwrap();
            let receipt = binding.resolved_temperature_source().unwrap();
            assert_eq!(receipt.query().axes()["T"], temperature);
            assert!(
                (receipt.property("mu_kinetic").unwrap().value_si()
                    - (0.4 - slope * (temperature - 293.15)))
                    .abs()
                    < 1e-13
            );
            assert_eq!(
                receipt
                    .property("mu_kinetic")
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                binding
                    .temperature_source()
                    .unwrap()
                    .initial_source()
                    .property("mu_kinetic")
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert!(
                frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j
            );
            assert!(temperature > 293.15);
            runs.push((pressure, temperature));
        }
        assert!(runs[0].1 > runs[1].1);
        assert!(
            runs[0]
                .0
                .iter()
                .zip(&runs[1].0)
                .any(|(a, b)| a.to_bits() != b.to_bits())
        );
        assert!(
            runs[0]
                .0
                .iter()
                .zip(&runs[2].0)
                .any(|(a, b)| a.to_bits() != b.to_bits())
        );
        let mut source = temperature_ready(cx, gate, &profile, 0.02, 0.005, 128);
        let prepared = profile.prepare(&source, gate).unwrap();
        let mut trace = Vec::new();
        let (wav, summary) = prepared
            .write(&mut source, Cursor::new(Vec::new()), Some(&mut trace), gate)
            .unwrap();
        assert_eq!(summary.samples, 128);
        assert_eq!(
            wav.into_inner(),
            fs_couple::pcm_wav::encode_pcm16_wav(&runs[0].0, 48_000, 0.01)
                .unwrap()
                .0
        );
        assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 129);
    });
}

#[test]
fn native_temperature_source_requires_complete_heat_routes_and_survives_ribbon_partition() {
    scoped(|cx, gate| {
        let mut spec = input(false, 37, 48_000);
        spec.friction
            .load(gate, |_| Ok(temperature_interface(0.005)))
            .unwrap();
        let mut reads = 0;
        let error = spec
            .compile(cx, gate, |name| {
                reads += 1;
                Ok(material(name, 20.0, None))
            })
            .err()
            .unwrap();
        assert!(error.to_string().contains("finite-cell dry-contact"));
        assert_eq!(reads, 0);

        let profile = Profile::parse(profile_text(false).as_bytes()).unwrap();
        let mut spec = input(false, 37, 48_000);
        profile.configure(&mut spec);
        spec.bow_ribbon = Some((0.001, 3));
        spec.friction
            .load(gate, |_| Ok(temperature_interface(0.005)))
            .unwrap();
        let source = spec
            .compile(cx, gate, |name| Ok(material(name, 20.0, None)))
            .unwrap();
        let first = source.performance.bank().bow_interface(0).unwrap();
        for port in 0..3 {
            let binding = source.performance.bank().bow_interface(port).unwrap();
            assert_eq!(
                binding.temperature_source().unwrap().identity(),
                first.temperature_source().unwrap().identity()
            );
            assert!((binding.nominal_area_m2() - 0.00002 / 3.0).abs() < 1e-20);
        }
        // This original one-port profile cannot silently leave two evolving laws unheated.
        let error = profile.prepare(&source, gate).err().unwrap();
        assert!(
            error
                .to_string()
                .contains("every temperature-dependent dry port")
        );
        assert_eq!(source.performance.samples_rendered(), 0);
    });
}
