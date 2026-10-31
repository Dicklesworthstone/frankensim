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

fn moving_film_spec(profile: &Profile, block: usize) -> input::Spec {
    let mut text = include_str!("../shear-layer-bow.material-strings")
        .split_once("\nschedule\n").unwrap().0
        .replace("audio 48000 4096 0.01 128", &format!("audio 48000 128 0.01 {block}"))
        .replace("bow_ribbon 0.01 3", "bow_ribbon 0.01 1")
        .replace("forces ", "bow_fibers fiber.fsmatpk 4e-5 0.6 180 40 0.4 2 2\n\
            bow_fiber_program 0 0.1 1000 0.01 1048576 3\n\
            fiber_speed 32 0.05\nfiber_speed 64 -0.05\nfiber_speed 96 0\nforces ");
    let schedule = GestureSchedule::try_new(48_000, vec![GestureTrack {
        id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
        initial: GestureValue::Bow { velocity_m_per_s: 0.0, normal_force_n: 0.01, station: 0.2 },
        events: vec![],
    }]).unwrap();
    text.push_str("\nschedule\n");
    text.push_str(std::str::from_utf8(&schedule.to_canonical_bytes()).unwrap());
    let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
    profile.configure(&mut spec);
    spec
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

#[test]
fn native_accelerating_fibers_keep_string_film_heat_history_and_replay_pressure() {
    scoped(|cx, gate| {
        let profile = Profile::parse(film_profile(false).replace("body.fsmatpk", "fiber.fsmatpk").as_bytes()).unwrap();
        let mut runs = Vec::new();
        for cp in [20.0, 2000.0] {
            let mut ready = moving_film_spec(&profile, 37).compile(cx, gate,
                |name| Ok(film_material(name, cp))).unwrap();
            let prepared = profile.prepare(&ready, gate).unwrap();
            let network = prepared.network(gate).unwrap();
            let routes = prepared.film_routes(&ready, &network).unwrap();
            let mut renderer = ready.performance.with_contact_thermal_network(network,
                film_bindings(&prepared), vec![], profile.controls, gate,
                |endpoint| prepared.resolve(endpoint)).unwrap()
                .with_trial_heat_source(|endpoint, ports| prepared.resolve_with_films(endpoint, ports, &routes)).unwrap();
            let initial_heat = renderer.thermal().energy_j();
            let mut pressure = vec![0.0; 128];
            renderer.block(&mut pressure[..32]).unwrap();
            let support = renderer.performance().fiber_support_frame().unwrap();
            assert_eq!(support.sample, 32);
            assert_eq!(support.velocity_m_s, 0.05);
            assert!((support.translation_m - 0.5 * 0.05 * 32.0 / 48_000.0).abs() < 1e-14);
            let first_layer = renderer.last_frame().unwrap().mechanics.ports[0]
                .shear_layer.as_ref().unwrap();
            assert!(first_layer.plastic_slip_m.abs() > 0.0);
            assert!(renderer.thermal().temperature_k(FILM).unwrap() > 293.15);
            for block in pressure[32..].chunks_mut(32) { renderer.block(block).unwrap(); }
            let support = *renderer.performance().fiber_support_frame().unwrap();
            assert_eq!(support.sample, 128);
            assert_eq!(support.velocity_m_s, 0.0);
            assert!(support.translation_m.abs() < 1e-14);
            assert_eq!(support.laboratory_kinetic_correction_j, 0.0);
            assert!(renderer.thermal().energy_j() > initial_heat);
            assert!([CellId(22), CellId(23)].into_iter()
                .any(|cell| renderer.thermal().temperature_k(cell).unwrap() > 293.15));
            let temperature = renderer.thermal().temperature_k(FILM).unwrap();
            assert!(temperature > 293.15);
            assert!(pressure.iter().any(|p| p.abs() > 1e-12));
            let frame = renderer.last_frame().unwrap();
            assert!(frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j);
            assert!(frame.mechanics.ports[0].shear_layer.is_some());
            let final_frame = frame.mechanics.clone();
            drop(renderer);

            let mut streamed = moving_film_spec(&profile, 128).compile(cx, gate,
                |name| Ok(film_material(name, cp))).unwrap();
            let stream = profile.prepare(&streamed, gate).unwrap();
            let mut trace = Vec::new();
            let (wav, summary) = stream.write(&mut streamed, Cursor::new(Vec::new()), Some(&mut trace), gate).unwrap();
            assert_eq!(wav.into_inner(), fs_couple::pcm_wav::encode_pcm16_wav(&pressure, 48_000, 0.01).unwrap().0);
            assert_eq!(summary.samples, 128);
            assert_eq!(String::from_utf8(trace).unwrap().lines().count(), 129);
            assert_eq!(streamed.performance.bank().system().last_frame(), Some(&final_frame));
            assert_eq!(streamed.performance.fiber_support_frame(), Some(&support));
            runs.push((pressure, temperature));
        }
        assert!(runs[0].1 > runs[1].1);
        assert!(runs[0].0.iter().zip(&runs[1].0).any(|(a, b)| a.to_bits() != b.to_bits()));
    });
}

#[test]
fn native_moving_film_requires_complete_patch_ownership_and_actual_fiber_material() {
    scoped(|cx, gate| {
        let complete = film_profile(false);
        let incomplete = format!("{}film-surfaces 0\nseries 0\nlinks 0\n",
            complete.split_once("film-surfaces 1").unwrap().0);
        let incomplete = Profile::parse(incomplete.as_bytes()).unwrap();
        assert!(moving_film_spec(&incomplete, 128).compile(cx, gate,
            |name| Ok(film_material(name, 20.0))).is_err());

        // The body pack exists in the same assembled instrument but is not the
        // material of the fibers actually touching this attached film.
        let wrong_surface = Profile::parse(complete.as_bytes()).unwrap();
        let ready = moving_film_spec(&wrong_surface, 128).compile(cx, gate,
            |name| Ok(film_material(name, 20.0))).unwrap();
        assert!(wrong_surface.prepare(&ready, gate).is_err());
        assert_eq!(ready.performance.samples_rendered(), 0);
    });
}

// Two different source laws act on the SAME string/body/pressure trajectory.
// Surface cells use distinct physical volumes; the explicit link carries heat
// between the film and the dry string patch in the same nonlinear candidate.
fn mixed_profile(narrow_surface: bool) -> String {
    let film = film_profile(narrow_surface);
    let (cells, routes) = film.split_once("dry-contacts 0\n").unwrap();
    let mut text = cells.replace("cells 5", "cells 13");
    for (side, pack) in ["wire.fsmatpk", "body.fsmatpk"].into_iter().enumerate() {
        for i in 0..4 {
            text.push_str(&format!(
                "cell {} {pack} 280 293.15 320 {} 1 {side} {} 1 {side} 4e-10 4.048e-5\n",
                30 + side * 10 + i, i as f64 * 0.002, (i + 1) as f64 * 0.002));
        }
    }
    text.push_str("dry-contacts 1\ndry-contact 1 -0.148 fixed\n\
        footprint 0.001 0.003 0.01 0.001 0.003 0.01\n");
    for side in 0..2 {
        text.push_str("surface 4\n");
        for i in 0..4 {
            text.push_str(&format!("strip {} {} {} 0.01 0.00001\n", 30 + side * 10 + i,
                i as f64 * 0.002, (i + 1) as f64 * 0.002));
        }
    }
    text.push_str(&routes.replace("links 0", "links 1\nlink 9 10 30 5e-8 0.00005 0.00001"));
    text
}

fn mixed_spec_text(block: usize) -> String {
    let mut text = include_str!("../shear-layer-bow.material-strings")
        .split_once("\nschedule\n").unwrap().0
        .replace("audio 48000 4096 0.01 128", &format!("audio 48000 128 0.01 {block}"))
        .replace("bow_ribbon 0.01 3\n", "")
        .replace("friction layers 1\nshear-layer bow", "friction contacts 2\n\
            bow-interface dry dry.fsintpk 293.15 0.00002 50 0 coulomb mu_static mu_kinetic\n\
            shear-layer film");
    // Deliberately reverse declaration order relative to physical contact order.
    let schedule = GestureSchedule::try_new(48_000, vec![
        GestureTrack { id: "film".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow { velocity_m_per_s: 0.1, normal_force_n: 0.01, station: 0.2 },
            events: vec![] },
        GestureTrack { id: "dry".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow { velocity_m_per_s: 0.03, normal_force_n: 0.001, station: 0.3 },
            events: vec![] },
    ]).unwrap();
    text.push_str("\nschedule\n");
    text.push_str(std::str::from_utf8(&schedule.to_canonical_bytes()).unwrap());
    text
}

fn mixed_ready(cx: &Cx<'_>, gate: &CancelGate, profile: &Profile, cp: f64, block: usize) -> Ready {
    let mut spec = input::Spec::parse(mixed_spec_text(block).as_bytes()).unwrap();
    profile.configure(&mut spec);
    spec.friction.load(gate, |_| Ok(temperature_interface(0.005))).unwrap();
    spec.compile(cx, gate, |name| Ok(film_material(name, cp))).unwrap()
}

#[test]
fn mixed_contacts_require_complete_disjoint_track_and_heat_ownership() {
    let text = mixed_spec_text(37);
    input::Spec::parse(text.as_bytes()).unwrap();
    for (from, to) in [
        ("bow-interface dry ", "bow-interface film "),
        ("bow-interface dry ", "bow-interface absent "),
        ("friction contacts 2", "friction interfaces 2"),
        ("friction contacts 2", "friction layers 2"),
    ] {
        assert!(input::Spec::parse(text.replace(from, to).as_bytes()).is_err(), "accepted {to}");
    }
    let heat = mixed_profile(false);
    Profile::parse(heat.as_bytes()).unwrap();
    assert!(Profile::parse(heat.replace("dry-contact 1", "dry-contact 0").as_bytes()).is_err());
}

#[test]
fn native_mixed_sources_share_heat_history_and_replay_the_same_pressure() {
    scoped(|cx, gate| {
        let profile = Profile::parse(mixed_profile(false).as_bytes()).unwrap();
        let mut responses = Vec::new();
        for cp in [20.0, 2000.0] {
            let mut ready = mixed_ready(cx, gate, &profile, cp, 37);
            assert!(ready.performance.bank().bow_shear_layer(0).is_some());
            assert!(ready.performance.bank().bow_interface(0).is_none());
            assert!(ready.performance.bank().bow_shear_layer(1).is_none());
            assert!(ready.performance.bank().bow_interface(1).unwrap().temperature_source().is_some());
            let prepared = profile.prepare(&ready, gate).unwrap();
            let network = prepared.network(gate).unwrap();
            let routes = prepared.film_routes(&ready, &network).unwrap();
            let dry = prepared.dry_routes(&ready).unwrap();
            let mut renderer = ready.performance.with_combined_contact_thermal_network(network,
                film_bindings(&prepared), vec![], dry, vec![], profile.controls, gate,
                |endpoint| prepared.resolve(endpoint)).unwrap()
                .with_trial_heat_source(|endpoint, ports| prepared.resolve_with_films(endpoint, ports, &routes)).unwrap();
            let initial_heat = renderer.thermal().energy_j();
            let mut deposited = 0.0;
            let mut dry_heat = 0.0;
            let mut film_heat = 0.0;
            let mut pressure = vec![0.0; 128];
            for sample in pressure.chunks_mut(1) {
                renderer.block(sample).unwrap();
                let frame = renderer.last_frame().unwrap();
                assert!(frame.mechanics.ports[0].shear_layer.is_some());
                assert!(frame.mechanics.ports[1].shear_layer.is_none());
                assert_eq!(frame.dry_contacts.len(), 1);
                let dry = &frame.dry_contacts[0];
                assert_eq!(dry.port, 1);
                assert!(dry.interface_temperature_k.is_some());
                assert!((dry.cell_heat_j.iter().map(|(_, q)| q).sum::<f64>() - dry.friction_heat_j).abs()
                    <= profile.controls.total_energy_tolerance_j);
                assert!(frame.total_energy_residual_j.abs() <= profile.controls.total_energy_tolerance_j);
                deposited += frame.deposited_heat_j;
                dry_heat += frame.routed_dry_heat_j;
                film_heat += frame.routed_plastic_heat_j + frame.routed_relaxation_heat_j;
            }
            assert!(dry_heat > 0.0 && film_heat > 0.0);
            assert!((renderer.thermal().energy_j() - initial_heat - deposited).abs() < 2e-9);
            assert!(renderer.thermal().temperature_k(CellId(30)).unwrap() > 293.15);
            assert!(renderer.thermal().temperature_k(CellId(40)).unwrap() > 293.15);
            assert!(renderer.thermal().temperature_k(FILM).unwrap() > 293.15);
            assert!(pressure.iter().any(|p| p.abs() > 1e-12));
            let temperature = renderer.thermal().temperature_k(FILM).unwrap();
            let final_frame = renderer.last_frame().unwrap().mechanics.clone();
            drop(renderer);

            let mut streamed = mixed_ready(cx, gate, &profile, cp, 128);
            let stream = profile.prepare(&streamed, gate).unwrap();
            let (wav, summary) = stream.write(&mut streamed, Cursor::new(Vec::new()), None, gate).unwrap();
            assert_eq!(summary.samples, 128);
            assert_eq!(wav.into_inner(), fs_couple::pcm_wav::encode_pcm16_wav(&pressure, 48_000, 0.01).unwrap().0);
            assert_eq!(streamed.performance.bank().system().last_frame(), Some(&final_frame));
            responses.push((temperature, pressure));
        }
        assert!(responses[0].0 > responses[1].0);
        assert!(responses[0].1.iter().zip(&responses[1].1).any(|(a, b)| a.to_bits() != b.to_bits()));
    });
}

#[test]
fn native_mixed_sources_cancel_inside_trial_without_publishing_either_history() {
    use std::sync::atomic::{AtomicBool, Ordering};
    scoped(|cx, gate| {
        let profile = Profile::parse(mixed_profile(false).as_bytes()).unwrap();
        let mut ready = mixed_ready(cx, gate, &profile, 20.0, 37);
        let prepared = profile.prepare(&ready, gate).unwrap();
        let network = prepared.network(gate).unwrap();
        let routes = prepared.film_routes(&ready, &network).unwrap();
        let dry = prepared.dry_routes(&ready).unwrap();
        let cancel_trial = AtomicBool::new(false);
        let mut renderer = ready.performance.with_combined_contact_thermal_network(network,
            film_bindings(&prepared), vec![], dry, vec![], profile.controls, gate,
            |endpoint| prepared.resolve(endpoint)).unwrap()
            .with_trial_heat_source(|endpoint, ports| {
                // Both mechanical source candidates and dry fluxes have already
                // been evaluated when the same trial's film conduction is read.
                if cancel_trial.load(Ordering::SeqCst) { gate.request(); }
                prepared.resolve_with_films(endpoint, ports, &routes)
            }).unwrap();
        renderer.block(&mut [0.0; 32]).unwrap();
        let count = renderer.samples_rendered();
        let cells = renderer.thermal().cells().clone();
        let frame = renderer.last_frame().unwrap().mechanics.clone();
        let system = renderer.performance().bank().system();
        let history = system.material_shear_state(0).cloned();
        let states = system.components().iter().map(|body| body.states().to_vec()).collect::<Vec<_>>();
        let dry_temperature = renderer.performance().bank().bow_interface(1).unwrap().temperature_kelvin();
        cancel_trial.store(true, Ordering::SeqCst);
        let mut output = [12345.0];
        assert!(renderer.block(&mut output).is_err());
        assert!(gate.is_requested());
        assert_eq!(output, [12345.0]);
        assert_eq!(renderer.samples_rendered(), count);
        assert_eq!(renderer.thermal().cells(), &cells);
        assert_eq!(renderer.last_frame().unwrap().mechanics, frame);
        let system = renderer.performance().bank().system();
        assert_eq!(system.material_shear_state(0), history.as_ref());
        assert_eq!(system.components().iter().map(|body| body.states().to_vec()).collect::<Vec<_>>(), states);
        assert_eq!(renderer.performance().bank().bow_interface(1).unwrap().temperature_kelvin(), dry_temperature);
    });
}
