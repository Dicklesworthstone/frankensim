//! Manufactured native input/heat/PCM integration, not measured bow or rosin data.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, PropertyClaim, PropertyValue, Provenance, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_scenario::gesture::{GestureTrack, GestureValue};
use std::io::Cursor;

fn scoped(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        fs_exec::StreamKey { seed: 0, kernel_id: 0, tile: 0, iteration: 0 },
        fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic), &gate));
}

fn claim(name: &str, dims: Dims, value: f64, validity: ValidityDomain) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
        validity, uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
        provenance: Provenance { source: "MANUFACTURED finite-cell dry-heat integration".into(),
            license: "CC0-1.0".into(), artifact: Some(fs_blake3::ContentHash([47; 32])) },
    }
}

fn material(name: &str, cp: f64, omit: Option<&str>) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let fiber = name == "fiber.fsmatpk";
    for (key, dims, value) in [
        ("density", Density::DIMS, if fiber { 1300.0 } else { 7800.0 }),
        ("young_modulus", Pressure::DIMS, if fiber { 3e9 } else { 200e9 }),
        ("poisson_ratio", Dims::NONE, 0.3),
        (ETA, DynViscosity::DIMS, 50.0),
        (fs_material::section::axial::YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS, 50.0),
        (SHEAR, Pressure::DIMS, 80e9),
        (SHEAR_ETA, DynViscosity::DIMS, 1000.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, cp),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.2),
    ] {
        if Some(key) != omit {
            claims.insert_claim(claim(key, dims, value, ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0).with("omega", 0.0, 400_000.0))).unwrap();
        }
    }
    NormalizedPack::new(name, "manufactured-v1", fs_blake3::ContentHash([47; 32]),
        "manufactured freely redistributable integration data", claims, vec![], vec![]).unwrap()
}

fn interface(receiver: &str) -> NormalizedInterfacePack {
    let mut claims = ClaimSet::new();
    for (key, value) in [("mu_static", 0.8), ("mu_kinetic", 0.4)] {
        claims.insert_claim(claim(key, Dims::NONE, value, ValidityDomain::unconstrained()
            .with("T", 280.0, 320.0).with("p", 0.0, 1e6).with("v", 0.0, 10.0))).unwrap();
    }
    let claims = NormalizedPack::new("manufactured-dry-contact", "manufactured-v1",
        fs_blake3::ContentHash([48; 32]), "manufactured CC0 fixture", claims, vec![], vec![]).unwrap();
    let surface = |name| SurfaceSpec {
        material: MaterialPack::from(material(name, 20.0, None)).card().unwrap().id().clone(),
        texture_frame: "manufactured tangent direction".into(),
    };
    NormalizedInterfacePack::new(surface("wire.fsmatpk"), surface(receiver), SystemContext {
        medium: "dry".into(), third_body: None, environment: "declared gas".into(),
        history: "manufactured unchanged dry surfaces".into(),
    }, claims).unwrap()
}

fn profile_text(fiber: bool) -> String {
    let receiver = if fiber { "fiber.fsmatpk" } else { "body.fsmatpk" };
    let mapping = if fiber { "fiber -0.238" } else { "fixed" };
    let mut text = String::from("frankensim-contact-thermal-v2\nboundary adiabatic\n\
        solver 512 24 1e-9 1e-12 1e-12 4e-12 1e-13 1e-10 8 1e-8\ncells 8\n");
    for (side, pack) in ["wire.fsmatpk", receiver].into_iter().enumerate() {
        for cell in 0..4 {
            let start = cell as f64 * 0.002;
            let end = (cell + 1) as f64 * 0.002;
            text.push_str(&format!("cell {} {pack} 280 293.15 320 {start} 0 {side} {end} 0 {side} 4e-10 4.048e-5\n", 10 + side * 10 + cell));
        }
    }
    text.push_str(&format!("dry-contacts 1\ndry-contact 0 -0.098 {mapping}\nfootprint 0.001 0.003 0.01 0.001 0.003 0.01\n"));
    for side in 0..2 {
        text.push_str("surface 4\n");
        for cell in 0..4 {
            text.push_str(&format!("strip {} {} {} 0.01 0.00001\n", 10 + side * 10 + cell,
                cell as f64 * 0.002, (cell + 1) as f64 * 0.002));
        }
    }
    text.push_str("films 0\nseries 0\nlinks 0\n");
    text
}

fn input(fiber: bool, block: usize, rate: u32) -> input::Spec {
    let original = if fiber { include_str!("../fiber-ribbon.material-strings") }
        else { include_str!("../pluck.material-strings") };
    let mut text = original.split_once("\nschedule\n").unwrap().0.to_owned();
    text = text.replace("audio 48000 4096 0.01 128", &format!("audio {rate} 128 0.01 {block}"))
        .replace("bow_ribbon 0.01 3", "bow_ribbon 0.01 1");
    let law = "friction interfaces 1\nbow-interface bow dry.fsintpk 293.15 0.00002 50 0 coulomb mu_static mu_kinetic";
    text = text.replace("friction stribeck 0.8 0.4 0.1 0", law)
        .replace("friction coulomb 0.8 0.4", law);
    let motion = if fiber {
        "bow_fiber_program 0 0.1 100 0.01 1048576 3\nfiber_speed 32 0.01\nfiber_speed 64 -0.01\nfiber_speed 96 0\n"
    } else {
        "bow_station_motion 10 0.0002 1048576 2\nbow_station_knot 32 0.208\nbow_station_knot 64 0.2\n"
    };
    text = text.replace("forces ", &format!("{motion}forces "));
    let schedule = GestureSchedule::try_new(48_000, vec![GestureTrack {
        id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
        initial: GestureValue::Bow { velocity_m_per_s: if fiber { 0.0 } else { 0.03 },
            normal_force_n: 0.001, station: 0.2 }, events: vec![],
    }]).unwrap();
    text.push_str("\nschedule\n");
    text.push_str(std::str::from_utf8(&schedule.to_canonical_bytes()).unwrap());
    input::Spec::parse(text.as_bytes()).unwrap()
}

fn ready(cx: &Cx<'_>, gate: &CancelGate, profile: &Profile, fiber: bool,
    cp: f64, block: usize, rate: u32, omit: Option<&str>) -> Ready {
    let mut source = input(fiber, block, rate);
    profile.configure(&mut source);
    source.friction.load(gate, |_| Ok(interface(if fiber { "fiber.fsmatpk" } else { "body.fsmatpk" }))).unwrap();
    source.compile(cx, gate, |name| Ok(material(name, cp, omit))).unwrap()
}

#[test]
fn native_dry_heat_follows_moving_and_revisited_material_cells_and_matches_pcm() {
    scoped(|cx, gate| {
        let profile = Profile::parse(profile_text(false).as_bytes()).unwrap();
        let mut direct = ready(cx, gate, &profile, false, 20.0, 37, 48_000, None);
        let prepared = profile.prepare(&direct, gate).unwrap();
        let routes = prepared.dry_routes(&direct).unwrap();
        let network = prepared.network(gate).unwrap();
        let mut coupled = direct.performance.with_dry_contact_thermal_network(network, routes, vec![],
            profile.controls, gate, |endpoint| prepared.resolve(endpoint)).unwrap();
        let mut expected = vec![0.0; 128];
        coupled.block(&mut expected[..32]).unwrap();
        let reached = coupled.last_frame().unwrap().dry_contacts[0].footprint.intervals_m[0];
        assert!((reached[0] - (0.005 - 0.5 * 0.004 / 32.0)).abs() < 1e-14);
        let retained = coupled.thermal().cells()[&CellId(12)].specific_h_j_kg;
        assert!(retained > 0.0);
        coupled.block(&mut expected[32..64]).unwrap();
        let returned = coupled.last_frame().unwrap().dry_contacts[0].footprint.intervals_m[0];
        assert!((returned[0] - (0.001 + 0.5 * 0.004 / 32.0)).abs() < 1e-14);
        let abandoned = coupled.thermal().cells()[&CellId(12)].specific_h_j_kg;
        assert!(abandoned > 0.0);
        coupled.block(&mut expected[64..96]).unwrap();
        coupled.block(&mut expected[96..]).unwrap();
        let final_heat = coupled.thermal().cells().values()
            .map(|cell| cell.mass_kg * cell.specific_h_j_kg).sum::<f64>();
        assert!(final_heat > 0.0);
        assert_eq!(coupled.thermal().cells()[&CellId(12)].specific_h_j_kg, abandoned);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        let temperatures: Vec<_> = coupled.thermal().cells().keys()
            .map(|&id| coupled.thermal().temperature_k(id).unwrap()).collect();
        assert!(temperatures.iter().all(|t| *t >= 293.15));
        let final_frame = coupled.last_frame().unwrap().mechanics.clone();
        drop(coupled);
        let encoded = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48_000, 0.01).unwrap();
        let mut streamed = ready(cx, gate, &profile, false, 20.0, 128, 48_000, None);
        let stream = profile.prepare(&streamed, gate).unwrap();
        let mut trace = Vec::new();
        let (wav, summary) = stream.write(&mut streamed, Cursor::new(Vec::new()), Some(&mut trace), gate).unwrap();
        assert_eq!(wav.into_inner(), encoded.0);
        assert_eq!(summary.samples, 128);
        assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 129);
        assert_eq!(streamed.performance.bank().system().last_frame(), Some(&final_frame));
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
            let mut renderer = source.performance.with_dry_contact_thermal_network(
                prepared.network(gate).unwrap(), routes, vec![], profile.controls, gate,
                |endpoint| prepared.resolve(endpoint)).unwrap();
            let mut output = vec![0.0; 128];
            renderer.block(&mut output[..32]).unwrap();
            let moved = renderer.last_frame().unwrap().dry_contacts[0].footprint.intervals_m[1];
            assert_ne!(moved, [0.001, 0.003]);
            for chunk in output[32..].chunks_mut(32) { renderer.block(chunk).unwrap(); }
            let heat: f64 = renderer.thermal().cells().values().map(|cell|
                cell.mass_kg * cell.specific_h_j_kg).sum();
            assert!(heat > 0.0);
            temperatures.push(renderer.thermal().temperature_k(CellId(20)).unwrap());
            let frame = renderer.performance().fiber_support_frame().unwrap();
            assert_eq!(frame.velocity_m_s, 0.0);
            assert_eq!(frame.sample, 128);
            assert!(output.iter().any(|p| p.abs() > 1e-12));
            drop(renderer);
            let mut streamed = ready(cx, gate, &profile, true, cp, 128, 48_000, None);
            let stream = profile.prepare(&streamed, gate).unwrap();
            let (wav, _) = stream.write(&mut streamed, Cursor::new(Vec::new()), None, gate).unwrap();
            assert_eq!(wav.into_inner(), fs_couple::pcm_wav::encode_pcm16_wav(&output, 48_000, 0.01).unwrap().0);
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
            ("dry-contact 0 -0.098 fixed", "dry-contact 0 -0.098 fiber -0.238"),
            ("links 0", "links 1\nlink 1 10 20 0.00002 0.00001 0.00001"),
        ] {
            let profile = Profile::parse(profile_text(false).replace(from, to).as_bytes()).unwrap();
            let source = ready(cx, gate, &profile, false, 20.0, 128, 48_000, None);
            assert!(profile.prepare(&source, gate).is_err(), "accepted {from} -> {to}");
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
    for (_, item) in original.card().claims().claims_ordered() {
        let mut item = item.clone();
        if item.key.name() == "mu_kinetic" {
            item.value = PropertyValue::Curve { abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: [280.0, 320.0].map(|t| (t, 0.4 - slope * (t - 293.15))).to_vec(),
                dims: Dims::NONE };
            item.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(item).unwrap();
    }
    let pack = NormalizedPack::new("manufactured-thermal-dry-contact", "manufactured-v1",
        fs_blake3::ContentHash([49; 32]), "manufactured CC0 fixture", claims, vec![], vec![]).unwrap();
    NormalizedInterfacePack::new(original.card().surface_a().clone(),
        original.card().surface_b().clone(), original.card().context().clone(), pack).unwrap()
}

fn temperature_ready(cx: &Cx<'_>, gate: &CancelGate, profile: &Profile,
    cp: f64, slope: f64, block: usize) -> Ready {
    let mut spec = input(false, block, 48_000);
    profile.configure(&mut spec);
    spec.friction.load(gate, |_| Ok(temperature_interface(slope))).unwrap();
    spec.compile(cx, gate, |name| Ok(material(name, cp, None))).unwrap()
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
            let mut render = source.performance.with_dry_contact_thermal_network(
                prepared.network(gate).unwrap(), routes, vec![], profile.controls, gate,
                |endpoint| prepared.resolve(endpoint)).unwrap();
            let mut pressure = vec![0.0; 128];
            for block in pressure.chunks_mut(37) { render.block(block).unwrap(); }
            let frame = render.last_frame().unwrap();
            let temperature = frame.dry_contacts[0].interface_temperature_k.unwrap();
            let binding = render.performance().bank().bow_interface(0).unwrap();
            let receipt = binding.resolved_temperature_source().unwrap();
            assert_eq!(receipt.query().axes()["T"], temperature);
            assert!((receipt.property("mu_kinetic").unwrap().value_si()
                - (0.4 - slope * (temperature - 293.15))).abs() < 1e-13);
            assert_eq!(receipt.property("mu_kinetic").unwrap().answer().receipt.selected,
                binding.temperature_source().unwrap().initial_source()
                    .property("mu_kinetic").unwrap().answer().receipt.selected);
            assert!(frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j);
            assert!(temperature > 293.15);
            runs.push((pressure, temperature));
        }
        assert!(runs[0].1 > runs[1].1);
        assert!(runs[0].0.iter().zip(&runs[1].0).any(|(a, b)| a.to_bits() != b.to_bits()));
        assert!(runs[0].0.iter().zip(&runs[2].0).any(|(a, b)| a.to_bits() != b.to_bits()));
        let mut source = temperature_ready(cx, gate, &profile, 0.02, 0.005, 128);
        let prepared = profile.prepare(&source, gate).unwrap();
        let mut trace = Vec::new();
        let (wav, summary) = prepared.write(&mut source, Cursor::new(Vec::new()),
            Some(&mut trace), gate).unwrap();
        assert_eq!(summary.samples, 128);
        assert_eq!(wav.into_inner(), fs_couple::pcm_wav::encode_pcm16_wav(
            &runs[0].0, 48_000, 0.01).unwrap().0);
        assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 129);
    });
}

#[test]
fn native_temperature_source_requires_complete_heat_routes_and_survives_ribbon_partition() {
    scoped(|cx, gate| {
        let mut spec = input(false, 37, 48_000);
        spec.friction.load(gate, |_| Ok(temperature_interface(0.005))).unwrap();
        let mut reads = 0;
        let error = spec.compile(cx, gate, |name| {
            reads += 1; Ok(material(name, 20.0, None))
        }).err().unwrap();
        assert!(error.to_string().contains("finite-cell dry-contact"));
        assert_eq!(reads, 0);

        let profile = Profile::parse(profile_text(false).as_bytes()).unwrap();
        let mut spec = input(false, 37, 48_000);
        profile.configure(&mut spec);
        spec.bow_ribbon = Some((0.001, 3));
        spec.friction.load(gate, |_| Ok(temperature_interface(0.005))).unwrap();
        let source = spec.compile(cx, gate, |name| Ok(material(name, 20.0, None))).unwrap();
        let first = source.performance.bank().bow_interface(0).unwrap();
        for port in 0..3 {
            let binding = source.performance.bank().bow_interface(port).unwrap();
            assert_eq!(binding.temperature_source().unwrap().identity(),
                first.temperature_source().unwrap().identity());
            assert!((binding.nominal_area_m2() - 0.00002 / 3.0).abs() < 1e-20);
        }
        // This original one-port profile cannot silently leave two evolving laws unheated.
        let error = profile.prepare(&source, gate).err().unwrap();
        assert!(error.to_string().contains("every temperature-dependent dry port"));
        assert_eq!(source.performance.samples_rendered(), 0);
    });
}
