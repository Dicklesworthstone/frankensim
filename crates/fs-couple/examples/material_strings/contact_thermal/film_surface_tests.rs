//! Manufactured retained-film/sliding-surface integration, not measured rosin.
use super::*;

const FILM: CellId = CellId(10);

fn film_material(name: &str, cp: f64) -> NormalizedPack {
    if name != "manufactured-shear-layer.fsmatpk" {
        return material(name, 20.0, None);
    }
    let mut claims = ClaimSet::new();
    for (key, dims, cold, hot) in [
        ("density", Density::DIMS, 1000.0, 1000.0),
        ("shear_modulus", Pressure::DIMS, 8e6, 8e6),
        ("shear_yield_stress", Pressure::DIMS, 200.0, 100.0),
        ("plastic_shear_viscosity", DynViscosity::DIMS, 6.0, 3.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, cp, cp),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.2, 0.2),
    ] {
        let mut source = claim(key, dims, cold, ValidityDomain::unconstrained()
            .with("T", 280.0, 320.0).with("omega", 0.0, 400_000.0));
        if cold != hot {
            source.value = PropertyValue::Curve { abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS, knots: vec![(280.0, cold), (320.0, hot)], dims };
            source.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(source).unwrap();
    }
    NormalizedPack::new(name, "manufactured-film-heat", fs_blake3::ContentHash([51; 32]),
        "manufactured CC0 integration data", claims, vec![], vec![]).unwrap()
}

fn film_profile(narrow_surface: bool) -> String {
    let mut text = String::from("frankensim-contact-thermal-v3\nboundary adiabatic\n\
        solver 512 24 1e-9 1e-12 1e-12 4e-12 1e-13 1e-10 8 1e-8\ncells 5\n\
        cell 10 manufactured-shear-layer.fsmatpk 280 293.15 320 0.00025 0 0 0.00035 0 0 1e-10 4.02e-6\n");
    for i in 0..4 {
        text.push_str(&format!("cell {} body.fsmatpk 280 293.15 320 {} 0 1 {} 0 1 2e-11 2.404e-6\n",
            20 + i, i as f64 * 0.0001, (i + 1) as f64 * 0.0001));
    }
    text.push_str("dry-contacts 0\nfilms 1\nfilm 0 10\nfilm-surfaces 1\n\
        film-surface 0 +1\nfilm-footprint 0.00025 0.00035 0.01 0.00005\n");
    text.push_str(if narrow_surface { "surface 2\n" } else { "surface 4\n" });
    for i in (if narrow_surface { 2 } else { 0 })..4 {
        text.push_str(&format!("strip {} {} {} 0.01 0.00001\n", 20 + i,
            i as f64 * 0.0001, (i + 1) as f64 * 0.0001));
    }
    text.push_str("series 0\nlinks 0\n");
    text
}

fn film_ready(cx: &Cx<'_>, gate: &CancelGate, profile: &Profile, cp: f64, block: usize) -> Ready {
    let mut text = include_str!("../shear-layer-bow.material-strings")
        .split_once("\nschedule\n").unwrap().0
        .replace("audio 48000 4096 0.01 128", &format!("audio 48000 128 0.01 {block}"))
        .replace("bow_ribbon 0.01 3", "bow_ribbon 0.01 1");
    let schedule = GestureSchedule::try_new(48_000, vec![GestureTrack {
        id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
        initial: GestureValue::Bow { velocity_m_per_s: 0.1, normal_force_n: 0.01, station: 0.2 },
        events: vec![],
    }]).unwrap();
    text.push_str("\nschedule\n");
    text.push_str(std::str::from_utf8(&schedule.to_canonical_bytes()).unwrap());
    let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
    profile.configure(&mut spec);
    spec.compile(cx, gate, |name| Ok(film_material(name, cp))).unwrap()
}

fn film_bindings(prepared: &Prepared) -> Vec<LayerThermalCell<'_>> {
    prepared.profile.films.iter().map(|&(port, cell)| {
        let index = prepared.profile.cells.iter().position(|c| c.id == cell).unwrap();
        LayerThermalCell { port, cell, material: &prepared.sources[index].card,
            caloric: &prepared.sources[index].curve, temperature_axis: "T" }
    }).collect()
}

#[test]
fn native_retained_film_heats_moving_receivers_and_replays_the_same_pressure() {
    scoped(|cx, gate| {
        let profile = Profile::parse(film_profile(false).as_bytes()).unwrap();
        let mut runs = Vec::new();
        for cp in [20.0, 2000.0] {
            let mut ready = film_ready(cx, gate, &profile, cp, 37);
            let prepared = profile.prepare(&ready, gate).unwrap();
            let network = prepared.network(gate).unwrap();
            let routes = prepared.film_routes(&ready, &network).unwrap();
            let mut renderer = ready.performance.with_contact_thermal_network(network,
                film_bindings(&prepared), vec![], profile.controls, gate,
                |endpoint| prepared.resolve(endpoint)).unwrap()
                .with_trial_heat_source(|endpoint, ports| prepared.resolve_with_films(endpoint, ports, &routes)).unwrap();
            let mut pressure = vec![0.0; 128];
            renderer.block(&mut pressure[..32]).unwrap();
            let original_receiver_heat = renderer.thermal().cells()[&CellId(23)].specific_h_j_kg;
            assert!(original_receiver_heat > 0.0);
            for block in pressure[32..].chunks_mut(32) { renderer.block(block).unwrap(); }
            let frame = renderer.last_frame().unwrap();
            let layer = frame.mechanics.ports[0].shear_layer.as_ref().unwrap();
            let slip = layer.elastic_slip_m + layer.plastic_slip_m;
            assert!(slip > 0.0001 && slip < 0.00025);
            assert!(renderer.thermal().cells()[&CellId(20)].specific_h_j_kg > 0.0);
            assert!(renderer.thermal().cells()[&CellId(23)].specific_h_j_kg > 0.0);
            let exchange = prepared.resolve_with_films(renderer.thermal(), &frame.mechanics.ports, &routes).unwrap();
            assert_eq!(exchange.sources_w.get(&CellId(23)).copied().unwrap_or(0.0), 0.0);
            assert!(frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j);
            assert!(pressure.iter().any(|p| p.abs() > 1e-12));
            let temperature = renderer.thermal().temperature_k(FILM).unwrap();
            assert!(temperature > 293.15);
            let final_frame = frame.mechanics.clone();
            drop(renderer);
            let mut streamed = film_ready(cx, gate, &profile, cp, 128);
            let stream = profile.prepare(&streamed, gate).unwrap();
            let mut trace = Vec::new();
            let (wav, summary) = stream.write(&mut streamed, Cursor::new(Vec::new()), Some(&mut trace), gate).unwrap();
            assert_eq!(wav.into_inner(), fs_couple::pcm_wav::encode_pcm16_wav(&pressure, 48_000, 0.01).unwrap().0);
            assert_eq!(summary.samples, 128);
            assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 129);
            assert_eq!(streamed.performance.bank().system().last_frame(), Some(&final_frame));
            runs.push((pressure, temperature));
        }
        assert!(runs[0].1 > runs[1].1);
        assert!(runs[0].0.iter().zip(&runs[1].0).any(|(a, b)| a.to_bits() != b.to_bits()));
    });
}

#[test]
fn native_retained_film_refuses_late_surface_exit_without_publishing_the_step() {
    scoped(|cx, gate| {
        let profile = Profile::parse(film_profile(true).as_bytes()).unwrap();
        let mut ready = film_ready(cx, gate, &profile, 20.0, 37);
        let prepared = profile.prepare(&ready, gate).unwrap();
        let network = prepared.network(gate).unwrap();
        let routes = prepared.film_routes(&ready, &network).unwrap();
        let mut renderer = ready.performance.with_contact_thermal_network(network,
            film_bindings(&prepared), vec![], profile.controls, gate,
            |endpoint| prepared.resolve(endpoint)).unwrap()
            .with_trial_heat_source(|endpoint, ports| prepared.resolve_with_films(endpoint, ports, &routes)).unwrap();
        for _ in 0..128 {
            let count = renderer.samples_rendered();
            let cells = renderer.thermal().cells().clone();
            let frame = renderer.last_frame().map(|f| f.mechanics.clone());
            let mut output = [12345.0];
            if renderer.block(&mut output).is_err() {
                assert!(count > 0);
                assert_eq!(renderer.samples_rendered(), count);
                assert_eq!(renderer.thermal().cells(), &cells);
                assert_eq!(renderer.last_frame().map(|f| &f.mechanics), frame.as_ref());
                assert_eq!(output, [12345.0]);
                assert!(renderer.block(&mut output).is_err());
                assert_eq!(renderer.thermal().cells(), &cells);
                return;
            }
        }
        panic!("sliding film unexpectedly stayed inside the truncated receiving surface");
    });
}
