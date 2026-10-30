//! Original-owner replay and material/upstream composition of the native reservoir.
//! Baseline geometry and conservation cases live in the sibling reservoir tests.
use super::*;

fn source(temperature_k: f64) -> String {
    let gas = reference(100_000.0);
    let density = gas.pressure / (gas.specific_gas_constant * temperature_k);
    fixture(gas, "dry-air 300 100000").replace(
        "initial rest\n",
        &format!(
            "reservoir_aperture 1 low {density:.17e} 0 0 0 {temperature_k} 3\n\
reservoir_cell 5 0 1\nreservoir_cell 6 0 1\nreservoir_cell 7 0 1\ninitial rest\n"
        ),
    )
}

#[test]
fn controlled_side_opening_replays_original_fluid_and_duct_and_survives_refused_trials() {
    let text = source(303.0);
    let mut coupled = read(&text, 8).into_renderer().system;
    let original = read(&text, 8).into_renderer().system;
    let mut fluid = original.zone;
    let mut duct = original.load;
    let gate = CancelGate::new_clock_free();
    for index in 0..12 {
        let pressure = if index < 6 { 100020.0 } else { 100000.0 };
        let temperature = if index < 6 { 303.0 } else { 300.0 };
        let frame = coupled.step_with_reservoir(pressure, temperature, &gate).unwrap();
        assert_eq!(frame.drive, NsfDgDrive::Reservoir {
            pressure_pa: pressure, temperature_k: temperature,
        });
        let state = fluid.reservoir_aperture_state_at(pressure, temperature).unwrap();
        let trial = fluid
            .preview_boundary_interval_with_reservoir(
                [0.0, frame.incoming_characteristic_m_s]
                    .map(NsfBoundaryDrive::IncomingCharacteristic),
                None,
                state,
                1,
                108,
                || false,
            )
            .unwrap();
        assert_eq!(trial.boundary_exchange(), &frame.faces);
        assert_eq!(
            trial.energy_ledger().open_energy_j,
            frame.flow.open_energy_j
        );
        assert_eq!(trial.reservoir_mass_out_kg(), frame.reservoir_mass_out_kg);
        let acoustic = duct.preview_step(frame.inlet_flow_m3_s).unwrap();
        assert_eq!(acoustic.output(), &frame.duct);
        trial.accept();
        acoustic.accept();
        assert_eq!(fluid.state_hash(), coupled.zone().state_hash());
        assert_eq!(duct.state(), coupled.load().state());
    }
    let before = (
        coupled.zone().state_hash(),
        coupled.load().state().clone(),
        coupled.accepted_steps(),
    );
    let refused = CancelGate::new_clock_free();
    refused.request();
    assert!(coupled.step_with_reservoir(100010.0, 301.0, &refused).is_err());
    assert_eq!(
        before,
        (
            coupled.zone().state_hash(),
            coupled.load().state().clone(),
            coupled.accepted_steps()
        )
    );
    assert!(coupled.zone().reservoir_aperture().is_some());
    coupled.step_with_reservoir(100010.0, 301.0, &gate).unwrap();
    assert_eq!(coupled.accepted_steps(), before.2 + 1);
}

#[test]
fn side_opening_retains_thermochemical_sources_and_responds_to_their_calorics() {
    let build = |cp_r| {
        let packs = support::packs("synthetic-alpha", 0.028, cp_r, 200.0);
        let bindings = [packs.binding()];
        let sources = ThermochemicalSources {
            species: &bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        };
        let declaration = "thermochemical 300 100000 mole 1 synthetic-alpha 1";
        let record = read_gas_record(
            &mut declaration.split_ascii_whitespace(),
            None,
            Some(&sources),
        )
        .unwrap();
        // Fix geometry and reservoir pressure/temperature independently of
        // calorics. Only the retained original NASA9 heat-capacity law changes.
        let density = record.state().pressure / (record.state().specific_gas_constant * 303.0);
        let text = fixture(reference(100_000.0), declaration).replace(
            "initial rest\n",
            &format!(
                "reservoir_aperture 1 low {density:.17e} 0 0 0 303 3\n\
reservoir_cell 5 0 1\nreservoir_cell 6 0 1\nreservoir_cell 7 0 1\ninitial rest\n"
            ),
        ).replace("boundary 1\nhold 32 -0.1", "boundary reservoir 2\nhold 16 100020 303\nhold 32 100000 300");
        let file = ResolvedDgFilePerformance::from_bytes_with_sources(
            text.as_bytes(),
            8,
            Arc::new(CancelGate::new_clock_free()),
            Some(&sources),
            |_| Err("no material wall attachment in this fixture".into()),
        )
        .unwrap();
        assert_eq!(file.gas_records(), &[record]);
        (
            file,
            [
                packs.species.content_hash(),
                packs.nasa9.content_hash(),
                packs.transport.content_hash(),
            ],
        )
        // Borrowed source packs, registry, and description all drop here.
    };
    let (baseline, baseline_packs) = build(3.5);
    let (higher_cp, higher_packs) = build(4.5);
    assert_eq!(
        baseline.renderer().system().zone().law().dx_m,
        higher_cp.renderer().system().zone().law().dx_m
    );
    assert_eq!(
        baseline.renderer().system().zone().law().dt_s,
        higher_cp.renderer().system().zone().law().dt_s
    );
    assert_eq!(
        baseline.renderer().system().zone().reservoir_aperture(),
        higher_cp.renderer().system().zone().reservoir_aperture()
    );
    assert_eq!(baseline_packs[0], higher_packs[0]);
    assert_eq!(baseline_packs[2], higher_packs[2]);
    assert_ne!(baseline_packs[1], higher_packs[1]);
    let (baseline_pressure, baseline) = render(baseline, 8);
    let (higher_pressure, higher_cp) = render(higher_cp, 8);
    assert_ne!(
        baseline_pressure, higher_pressure,
        "original calorics must reach the solved pressure with the aperture attached"
    );
    for (renderer, identities) in [(&baseline, baseline_packs), (&higher_cp, higher_packs)] {
        let material = renderer
            .system()
            .zone()
            .thermochemical_state([6, 0, 1])
            .unwrap()
            .unwrap();
        let retained = &material.source_packs()[0];
        assert_eq!(retained.species_pack, identities[0]);
        assert_eq!(retained.nasa9_pack, identities[1]);
        assert_eq!(retained.transport_pack, identities[2]);
        let opening = renderer.system().zone().reservoir_aperture().unwrap();
        assert_eq!(opening.state.temperature_k, 300.0);
        assert_eq!(opening.state.velocity_m_s, [0.0; 3]);
        let frame = renderer.last_frame().unwrap();
        assert_eq!(frame.drive, NsfDgDrive::Reservoir { pressure_pa: 100000.0, temperature_k: 300.0 });
        assert!(frame.flow.open_energy_j.abs() > 1e-15);
        assert_eq!(frame.flow.viscous_heating_j, 0.0);
        assert!(frame.pressure_residual_pa.abs() < 1e-6);
        let mass_balance = frame.flow_mass_change_kg
            + frame.reservoir_mass_out_kg
            + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
        assert!(mass_balance.abs() < 1e-16);
        let supplied = frame.flow.wall_heat_j
            + frame.flow.open_boundary_energy_j
            + frame.flow.open_energy_j
            + frame.flow.source_energy_j;
        assert!((frame.flow_energy_change_j - supplied).abs() < 2e-11);
    }
}

#[test]
fn side_opening_keeps_finite_upstream_feedback_and_all_owners_on_refusal() {
    let input = |termination| {
        let gas = reference(100_000.0);
        let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
        let radius = (9.0 * dx * dx / core::f64::consts::PI).sqrt();
        source(300.0)
            .replace("coupling 64 6912", "coupling 256 27648")
            .replace("observation inlet-pressure", &format!(
                "upstream_dg dg_tube 0.03 {radius:.17e} {termination} 1 2 32 50 2000 0.01 1e-11 512\nobservation inlet-pressure"))
            .replace("boundary 1", "boundary reservoir 1")
            .replace("hold 32 -0.1", "hold 32 100020 300")
    };
    let mut closed = read(&input("closed"), 8).into_renderer().system;
    let mut matched = read(&input("matched"), 8).into_renderer().system;
    let gate = CancelGate::new_clock_free();
    let mut closed_pressure = Vec::new();
    let mut matched_pressure = Vec::new();
    let mut aperture_exchanged_mass = false;
    for index in 0..28 {
        for (system, pressure) in [
            (&mut closed, &mut closed_pressure),
            (&mut matched, &mut matched_pressure),
        ] {
            let temperature = if index < 14 { 300.0 } else { 300.15 };
            let frame = system.step_with_reservoir(100020.0, temperature, &gate).unwrap();
            let upstream = frame.upstream.as_ref().unwrap();
            pressure.push(frame.inlet_pressure_pa);
            assert!(frame.pressure_residual_pa.abs() < 1e-6);
            assert!(upstream.pressure_residual_pa.abs() < 1e-6);
            assert_eq!(upstream.source_work_j, 0.0, "side blowing adds no axial series source");
            assert!(
                (frame.faces[upstream.face].mean_acoustic_pressure_pa
                    - upstream.inlet_pressure_pa)
                    .abs()
                    < 1e-6
            );
            aperture_exchanged_mass |= frame.reservoir_mass_out_kg.abs() > 1e-20;
            let mass_balance = frame.flow_mass_change_kg
                + frame.reservoir_mass_out_kg
                + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
            assert!(mass_balance.abs() < 1e-16);
            let supplied = frame.flow.wall_heat_j
                + frame.flow.open_boundary_energy_j
                + frame.flow.open_energy_j
                + frame.flow.source_energy_j;
            assert!((frame.flow_energy_change_j - supplied).abs() < 2e-11);
        }
    }
    assert!(aperture_exchanged_mass);
    assert_ne!(
        closed_pressure, matched_pressure,
        "upstream acoustic feedback must still reach the main bore with a side opening"
    );
    let before = (
        closed.zone().state_hash(),
        closed.zone().cumulative_ledger(),
        closed.zone().last_reservoir_mass_out_kg(),
        closed.load().state().clone(),
        closed.upstream_load().unwrap().state().clone(),
        closed.accepted_steps(),
        closed.zone().reservoir_aperture().cloned(),
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    // The finite but unreachable reservoir pressure enters original fluid
    // trials before refusal; accepted reservoir conditions must also roll back.
    for (pressure, temperature, control) in [(100020.0, 300.15, &cancelled),
        (1e6, 300.15, &gate), (100020.0, 1000.0, &gate)] {
        assert!(
            closed
                .step_with_reservoir(pressure, temperature, control)
                .is_err()
        );
        assert_eq!(
            before,
            (
                closed.zone().state_hash(),
                closed.zone().cumulative_ledger(),
                closed.zone().last_reservoir_mass_out_kg(),
                closed.load().state().clone(),
                closed.upstream_load().unwrap().state().clone(),
                closed.accepted_steps(),
                closed.zone().reservoir_aperture().cloned(),
            )
        );
    }
    let frame = closed.step_with_reservoir(100020.0, 300.15, &gate).unwrap();
    assert_eq!(closed.accepted_steps(), before.5 + 1);
    assert!(frame.upstream.is_some());
    assert!(closed.zone().reservoir_aperture().is_some());
}

#[test]
fn reservoir_schedule_is_partition_invariant_and_owns_pressure_and_temperature() {
    let text = source(300.0).replace("boundary 1\nhold 32 -0.1",
        "boundary reservoir 3\nhold 8 100020 300\nhold 16 100020 300.15\nhold 32 100000 300");
    let file = read(&text, 8);
    let history = file.renderer().drive_history();
    for (index, pressure, temperature) in [(0, 100020.0, 300.0),
        (8, 100020.0, 300.15), (16, 100000.0, 300.0)]
    {
        assert_eq!(history[index], NsfDgDrive::Reservoir {
            pressure_pa: pressure, temperature_k: temperature,
        });
    }
    let (pressure, whole) = render(file, 8);
    let (partitioned, split) = render(read(&text, 3), 3);
    assert_eq!(pressure, partitioned);
    assert_eq!(whole.system().zone().state_hash(), split.system().zone().state_hash());
    assert!(pressure.iter().any(|p| p.abs() > 1e-6), "side pressure must reach the bore");
    let opening = whole.system().zone().reservoir_aperture().unwrap();
    assert_eq!(opening.state.velocity_m_s, [0.0; 3], "no scheduled jet velocity");
    assert_eq!(opening.state.temperature_k, 300.0);
    assert_eq!(opening.state.density_kg_m3, reference(100000.0).density);
    let unheated = text.replace("hold 16 100020 300.15", "hold 16 100020 300");
    assert_ne!(pressure, render(read(&unheated, 8), 8).0,
        "held temperature must change original pressure, not only metadata");
    let no_opening = fixture(reference(100000.0), "dry-air 300 100000")
        .replace("boundary 1\nhold 32 -0.1", "boundary reservoir 1\nhold 32 100020 300");
    for invalid in [no_opening,
        text.replace("hold 8 100020 300", "hold 8 0 300"),
        text.replace("hold 8 100020 300", "hold 8 100020 -1"),
        text.replace("hold 8 100020 300", "hold 8 100020")]
    {
        assert!(ResolvedDgFilePerformance::from_bytes(invalid.as_bytes(), 8,
            Arc::new(CancelGate::new_clock_free())).is_err());
    }
}

#[test]
fn published_transverse_source_keeps_both_acoustic_loads_in_native_playback() {
    let text = include_str!("../../../../../examples/resolved-dg-transverse-drive.performance")
        .replace("output 2048 20", "output 32 20")
        .replace("hold 128", "hold 8")
        .replace("hold 256", "hold 16")
        .replace("hold 2048", "hold 32");
    let (_, renderer) = render(read(&text, 4), 4);
    let frame = renderer.last_frame().unwrap();
    assert!(renderer.system().upstream_load().is_some());
    assert_eq!(frame.upstream.as_ref().unwrap().source_work_j, 0.0);
    assert_eq!(frame.drive, NsfDgDrive::Reservoir { pressure_pa: 100000.0, temperature_k: 300.0 });
    assert_eq!(renderer.system().zone().steps(), 32);
    assert!(renderer.reservoir_receiver().is_some());
    assert!(renderer.reservoir_observer_gas().is_some());
    assert_eq!(frame.reservoir_flows.len(), 1);
    assert_eq!(frame.reservoir_flows[0].cell, [5, 2, 1]);
}
