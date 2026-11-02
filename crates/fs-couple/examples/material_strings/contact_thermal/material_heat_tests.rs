//! Manufactured source/body/pressure integration; no identified instrument data.
use super::*;

fn profile() -> String {
    let volume = core::f64::consts::PI * (40e-6_f64).powi(2) * 0.6 * 180.0;
    format!("frankensim-contact-thermal-v4\nboundary adiabatic\n\
        solver 512 24 1e-9 1e-12 1e-12 4e-12 1e-13 1e-10 1 1e-8\n\
        cells 1\ncell 30 fiber.fsmatpk 280 293.15 320 0 0 0 0.6 0 0 {volume:.17e} 0.027\n\
        dry-contacts 0\nfilms 0\nfilm-surfaces 0\nseries 0\n\
        material-heat-budget 256 67108864 100000000000\nfiber-heat 1\nfiber-body 0 0 30\nlinks 0\n")
}

fn fiber_pack(cp: f64) -> NormalizedPack {
    fiber_pack_with_temperature(cp, false)
}

fn fiber_pack_with_temperature(cp: f64, temperature_rates: bool) -> NormalizedPack {
    let original = material("fiber.fsmatpk", cp, None);
    let mut claims = original.claims().clone();
    for (name, dims, value) in [
        ("fiber_delta", Pressure::DIMS, 1e9),
        ("fiber_tau", fs_qty::Time::DIMS, 0.0005),
    ] {
        let mut property = claim(name, dims, value, ValidityDomain::unconstrained()
            .with("T", 280.0, 320.0).with("omega", 0.0, 400_000.0));
        if name == "fiber_tau" && temperature_rates {
            // Reuse the manufactured tau(T) law in material_multi_bow's
            // finite fiber-heat fixture. This is not empirical horsehair data.
            property.value = PropertyValue::Curve { abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(280.0, 0.0006), (340.0, 0.0003)], dims };
            property.validity = ValidityDomain::unconstrained()
                .with("T", 280.0, 340.0).with("omega", 0.0, 400_000.0);
            property.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(property).unwrap();
    }
    NormalizedPack::new("manufactured-hereditary-thermal-fiber", "fixture",
        fs_blake3::ContentHash([53; 32]), "manufactured CC0 integration data",
        claims, vec![], vec![]).unwrap()
}

fn source(cx: &Cx<'_>, gate: &CancelGate, profile: &Profile,
    full_strain: bool, cp: f64, block: usize) -> Ready {
    source_with_pack(cx, gate, profile, full_strain, cp, block, fiber_pack(cp))
}

fn source_with_pack(cx: &Cx<'_>, gate: &CancelGate, profile: &Profile,
    full_strain: bool, other_cp: f64, block: usize, pack: NormalizedPack) -> Ready {
    let strain = if full_strain {
        "bow_fiber_strain 0.2 0.01 128 1e8 1e-10 1e-9\n"
    } else { "" };
    let text = include_str!("../fiber-ribbon.material-strings")
        .replace("audio 48000 4096 0.01 128", &format!("audio 48000 16 0.01 {block}"))
        .replace("bow_ribbon 0.01 3", "bow_ribbon 0.01 1")
        .replace("forces retain", &format!(
            "bow_fiber_memory 16 young_modulus 1 fiber_delta fiber_tau\n{strain}forces retain"));
    let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
    profile.configure(&mut spec);
    spec.compile(cx, gate, |name| Ok(if name == "fiber.fsmatpk" { pack.clone() }
        else { material(name, other_cp, None) })).unwrap()
}

#[test]
fn v4_material_route_rejects_duplicate_bodies_aliases_missing_cells_and_unbounded_work() {
    let original = profile();
    let parsed = Profile::parse(original.as_bytes()).unwrap();
    assert_eq!(parsed.material_heat.as_ref().unwrap().fibers.len(), 1);
    assert!(!parsed.material_heat.as_ref().unwrap().fibers[0].source_temperature);
    let selected = Profile::parse(original.replace("fiber-body 0 0 30",
        "fiber-body 0 0 30 source-temperature").as_bytes()).unwrap();
    assert!(selected.material_heat.as_ref().unwrap().fibers[0].source_temperature);
    for (from, to) in [
        ("fiber-heat 1\nfiber-body 0 0 30", "fiber-heat 2\nfiber-body 0 0 30\nfiber-body 0 0 30"),
        ("fiber-body 0 0 30", "fiber-body 0 0 31"),
        ("fiber-body 0 0 30", "fiber-body 0 0 30 fitted-temperature"),
        ("fiber-body 0 0 30", "fiber-body 0 0 30 source-temperature extra"),
        ("films 0", "films 1\nfilm 0 30"),
        ("material-heat-budget 256", "material-heat-budget 257"),
        ("67108864 100000000000", "0 100000000000"),
        ("67108864 100000000000", "67108864 0"),
    ] {
        assert!(Profile::parse(original.replace(from, to).as_bytes()).is_err(), "{from} -> {to}");
    }
    // V3 cannot silently consume new records as an unrecognized optional tail.
    assert!(Profile::parse(original.replace("-v4", "-v3").as_bytes()).is_err());
    // Existing versions still need a declared physical heat route.
    assert!(Profile::parse(original.replace("fiber-heat 1\nfiber-body 0 0 30", "fiber-heat 0")
        .as_bytes()).is_err());
}

#[test]
fn v4_fiber_body_heat_reaches_finite_mass_and_the_existing_pcm_route() {
    scoped(|cx, gate| {
        let profile = Profile::parse(profile().as_bytes()).unwrap();
        for full_strain in [false, true] {
            let mut direct = source(cx, gate, &profile, full_strain, 20.0, 7);
            let prepared = profile.prepare(&direct, gate).unwrap();
            let network = prepared.network(gate).unwrap();
            let (routes, rates) = prepared.material_routes(&direct, &network, gate).unwrap();
            assert!(rates.is_empty());
            let volume = core::f64::consts::PI * (40e-6_f64).powi(2) * 0.6 * 180.0;
            assert!((routes[0].volume_m3 - volume).abs() < 1e-14 * volume);
            let mut renderer = direct.performance.with_material_contact_thermal_network(
                network, vec![], vec![], routes, profile.controls, gate,
                |endpoint| prepared.resolve(endpoint)).unwrap();
            let mut pressure = [0.0; 16];
            let mut loss = 0.0;
            for sample in &mut pressure {
                renderer.block(std::slice::from_mut(sample)).unwrap();
                let frame = renderer.last_frame().unwrap();
                loss += frame.routed_material_heat_j;
                assert!(frame.routed_material_heat_j >= 0.0);
                assert!(frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j);
            }
            let cell = &renderer.thermal().cells()[&CellId(30)];
            assert!(loss > 0.0);
            assert!((cell.mass_kg * cell.specific_h_j_kg - loss).abs() < 2e-12);
            assert!(renderer.thermal().temperature_k(CellId(30)).unwrap() > 293.15);
            let accepted = renderer.last_frame().unwrap().mechanics.clone();
            drop(renderer);

            let mut streamed = source(cx, gate, &profile, full_strain, 20.0, 16);
            let prepared = profile.prepare(&streamed, gate).unwrap();
            let mut trace = Vec::new();
            let (wav, summary) = prepared.write(&mut streamed, Cursor::new(Vec::new()),
                Some(&mut trace), gate).unwrap();
            assert_eq!(wav.into_inner(),
                fs_couple::pcm_wav::encode_pcm16_wav(&pressure, 48_000, 0.01).unwrap().0);
            assert_eq!(summary.samples, 16);
            assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 17);
            assert_eq!(streamed.performance.bank().system().last_frame(), Some(&accepted));
        }
    });
}

#[test]
fn v4_rejects_wrong_source_volume_missing_region_and_budget_before_advancing() {
    scoped(|cx, gate| {
        for (from, to) in [
            ("cell 30 fiber.fsmatpk", "cell 30 wire.fsmatpk"),
            ("fiber-body 0 0 30", "fiber-body 1 0 30"),
            ("fiber-body 0 0 30", "fiber-body 0 1 30"),
            ("material-heat-budget 256", "material-heat-budget 1"),
        ] {
            let profile = Profile::parse(profile().replace(from, to).as_bytes()).unwrap();
            let ready = source(cx, gate, &profile, false, 20.0, 16);
            assert!(profile.prepare(&ready, gate).is_err(), "{from} -> {to}");
            assert_eq!(ready.performance.samples_rendered(), 0);
        }
        let mut wrong_volume = Profile::parse(profile().as_bytes()).unwrap();
        wrong_volume.cells[0].geometry.volume_m3 *= 2.0;
        let ready = source(cx, gate, &wrong_volume, true, 20.0, 16);
        assert!(wrong_volume.prepare(&ready, gate).is_err());
        assert_eq!(ready.performance.samples_rendered(), 0);
    });
}

#[test]
fn v4_source_temperature_rates_require_explicit_opt_in_and_linear_material_history() {
    scoped(|cx, gate| {
        let fixed = Profile::parse(profile().as_bytes()).unwrap();
        let ready = source_with_pack(cx, gate, &fixed, false, 20.0, 16,
            fiber_pack_with_temperature(20.0, true));
        assert!(fixed.prepare(&ready, gate).is_err());
        assert_eq!(ready.performance.samples_rendered(), 0);

        let selected = Profile::parse(profile().replace("fiber-body 0 0 30",
            "fiber-body 0 0 30 source-temperature").as_bytes()).unwrap();
        let nonlinear = source_with_pack(cx, gate, &selected, true, 20.0, 16,
            fiber_pack_with_temperature(20.0, true));
        assert!(selected.prepare(&nonlinear, gate).is_err());
        assert_eq!(nonlinear.performance.samples_rendered(), 0);
    });
}

#[test]
fn v4_source_temperature_fiber_rates_reach_original_pressure_and_pcm() {
    scoped(|cx, gate| {
        let profile = Profile::parse(profile().replace("fiber-body 0 0 30",
            "fiber-body 0 0 30 source-temperature").as_bytes()).unwrap();
        let mut observed = Vec::new();
        // Change ONLY the fiber material's caloric source. Geometry, mechanical
        // coefficients at the initial temperature, other cards and drive agree.
        for cp in [20.0, 40.0] {
            let pack = fiber_pack_with_temperature(cp, true);
            let mut direct = source_with_pack(cx, gate, &profile, false, 20.0, 7, pack.clone());
            let initial_time = direct.performance.bank().system()
                .material_relaxation_branches()[0].relaxation_time_s;
            let prepared = profile.prepare(&direct, gate).unwrap();
            assert_eq!(direct.performance.samples_rendered(), 0);
            let network = prepared.network(gate).unwrap();
            let (routes, rates) = prepared.material_routes(&direct, &network, gate).unwrap();
            assert_eq!(rates.len(), 1);
            assert!(rates[0].times.iter().all(|(_, name)| name == "fiber_tau"));
            let mut renderer = direct.performance.with_material_contact_thermal_network(
                network, vec![], vec![], routes, profile.controls, gate,
                |endpoint| prepared.resolve(endpoint)).unwrap()
                .with_material_relaxation_temperature(rates).unwrap();
            let mut pressure = [0.0; 16];
            let mut loss = 0.0;
            for sample in &mut pressure[..8] {
                renderer.block(std::slice::from_mut(sample)).unwrap();
                let frame = renderer.last_frame().unwrap();
                loss += frame.routed_material_heat_j;
                assert!(frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j);
            }
            // Resume the SAME performance with its accepted finite mass, updated
            // source snapshots and heat integrals. Recreating initial sources
            // here would either refuse or silently reset physical memory rates.
            let (network, routes, rates) = renderer.into_thermal_with_material_rates();
            assert_eq!(routes[0].source.state_coordinate("T"), Some(network.temperature_k(CellId(30)).unwrap()));
            let mut renderer = direct.performance.with_material_contact_thermal_network(
                network, vec![], vec![], routes, profile.controls, gate,
                |endpoint| prepared.resolve(endpoint)).unwrap()
                .with_material_relaxation_temperature(rates).unwrap();
            for sample in &mut pressure[8..] {
                renderer.block(std::slice::from_mut(sample)).unwrap();
                let frame = renderer.last_frame().unwrap();
                loss += frame.routed_material_heat_j;
                assert!(frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j);
            }
            let temperature = renderer.thermal().temperature_k(CellId(30)).unwrap();
            let time = renderer.performance().bank().system()
                .material_relaxation_branches()[0].relaxation_time_s;
            assert!(temperature > 293.15 && time < initial_time);
            let cell = &renderer.thermal().cells()[&CellId(30)];
            assert!(loss > 0.0);
            assert!((cell.mass_kg * cell.specific_h_j_kg - loss).abs() < 2e-12);
            let accepted = renderer.last_frame().unwrap().mechanics.clone();
            drop(renderer);

            // The declarative path, tracing and PCM writer consume that same
            // rate-enabled owner; no second simulation or thermal audio effect.
            let mut streamed = source_with_pack(cx, gate, &profile, false, 20.0, 16, pack);
            let prepared = profile.prepare(&streamed, gate).unwrap();
            let mut trace = Vec::new();
            let (wav, summary) = prepared.write(&mut streamed, Cursor::new(Vec::new()),
                Some(&mut trace), gate).unwrap();
            assert_eq!(wav.into_inner(),
                fs_couple::pcm_wav::encode_pcm16_wav(&pressure, 48_000, 0.01).unwrap().0);
            assert_eq!(summary.samples, 16);
            assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 17);
            assert_eq!(streamed.performance.bank().system().last_frame(), Some(&accepted));
            observed.push((pressure, temperature, time));
        }
        assert_eq!(observed[0].0[0], observed[1].0[0]); // same initial physical operator
        assert!(observed[0].1 > observed[1].1); // larger finite capacity warms less
        assert!(observed[0].2 < observed[1].2); // the ORIGINAL source tau(T) responds
        let difference: f64 = observed[0].0.iter().zip(&observed[1].0).map(|(a, b)| (a-b).abs()).sum();
        assert!(difference > 0.0, "finite-body caloric source must reach original physical pressure");
    });
}
