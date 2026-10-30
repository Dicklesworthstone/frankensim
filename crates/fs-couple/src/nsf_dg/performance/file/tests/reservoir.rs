//! Native side-reservoir records retain the original ghost and flux owners.
//! The small pressure differences below are numerical integration fixtures,
//! not a qualified flute aperture, jet or radiation boundary.
use super::*;

#[test]
fn reference_reservoir_state_uses_the_admitted_material_state() {
    let states = [
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap(),
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 315.0, 75000.0).unwrap(),
        GasState::try_new(&GasSpec::water_vapor_nist(), 350.0, 75000.0).unwrap(),
    ];
    for gas in states {
        let mut row = Row { fields: "reference 3".split_ascii_whitespace(), line: 1 };
        let state = row.reservoir_state(&gas).ok().expect("source reference");
        assert_eq!(state.density_kg_m3, gas.density);
        assert_eq!(state.temperature_k, gas.temperature);
        assert_eq!(state.velocity_m_s, [0.0; 3]);
        assert_eq!(row.number::<usize>().ok(), Some(3));
        assert!(row.end().is_ok());
        let mut physical = Row {
            fields: "reference 0.005 0.0095 0.0005 0.002".split_ascii_whitespace(), line: 1,
        };
        assert_eq!(physical.reservoir_state(&gas).ok(), Some(state));
        let bounds = std::array::from_fn::<_, 4, _>(|_| physical.scalar().ok().expect("physical bounds"));
        assert_eq!(bounds, [0.005, 0.0095, 0.0005, 0.002]);
        assert!(physical.end().is_ok());
        // An independent authored source keeps its physical state after the
        // surrounding reference material changes.
        let mut explicit = Row { fields: "1.5 0.1 -0.2 0.3 320 2".split_ascii_whitespace(), line: 1 };
        assert_eq!(explicit.reservoir_state(&gas).ok(), Some(NsfOpenState {
            density_kg_m3: 1.5, velocity_m_s: [0.1, -0.2, 0.3], temperature_k: 320.0,
        }));
        assert_eq!(explicit.number::<usize>().ok(), Some(2));
        assert!(explicit.end().is_ok());
    }
}

#[test]
fn reference_reservoir_refuses_ambiguous_or_nonphysical_records() {
    let gas = reference(100000.0);
    for text in ["reference", "reference 3 extra", "reference NaN", "reference -1",
        "ref 3", "NaN 0 0 0 300 3", "-1 0 0 0 300 3", "1 0 0 0 0 3"] {
        let mut row = Row { fields: text.split_ascii_whitespace(), line: 1 };
        let result = row.reservoir_state(&gas).and_then(|_| row.number::<usize>())
            .and_then(|_| row.end());
        assert!(result.is_err(), "admitted {text}");
    }
}

#[test]
fn reference_reservoir_pressure_edits_keep_geometry_and_baseline_replay() {
    let text = input(100000.0).replace(&record(100000.0), "reservoir_aperture 1 low reference 3");
    let (explicit, old) = render(read(&input(100000.0), 8), 8);
    let (bound, reference) = render(read(&text, 3), 3);
    assert_eq!(explicit, bound);
    assert_eq!(old.system().zone().state_hash(), reference.system().zone().state_hash());
    for pressure in [75000.0, 125000.0] {
        let edited = text.replace("gas dry-air 300 100000", &format!("gas dry-air 300 {pressure}"));
        let file = read(&edited, 8);
        let zone = file.renderer().system().zone();
        let gas = file.gas_records()[0].state();
        let opening = zone.reservoir_aperture().unwrap();
        assert_eq!(opening.state.density_kg_m3, gas.density);
        assert_eq!(opening.state.temperature_k, gas.temperature);
        assert_eq!(opening.cells, vec![[7, 0, 1], [5, 0, 1], [6, 0, 1]]);
        assert_eq!(zone.law().dx_m, old.system().zone().law().dx_m);
        let physical = edited.replace(
            "reservoir_aperture 1 low reference 3\nreservoir_cell 7 0 1\nreservoir_cell 5 0 1\nreservoir_cell 6 0 1",
            &format!("reservoir_aperture_m 1 low reference {:.17e} {:.17e} {:.17e} {:.17e}",
                5.0*zone.law().dx_m, 8.0*zone.law().dx_m, zone.law().dx_m, 2.0*zone.law().dx_m),
        );
        let metric = read(&physical, 8);
        let metric_zone = metric.renderer().system().zone();
        assert_eq!(metric_zone.reservoir_aperture().unwrap().state, opening.state);
        assert_eq!(metric_zone.last_reservoir_flows(), zone.last_reservoir_flows());
    }
}

#[test]
fn reference_reservoir_follows_complete_thermochemical_sources() {
    let mut densities = Vec::new();
    for molar_mass in [0.028, 0.044] {
        let packs = support::packs("synthetic-alpha", molar_mass, 3.5, 200.0);
        let bindings = [packs.binding()];
        let sources = ThermochemicalSources { species: &bindings, water_species: None,
            phase: ThermochemicalPhase::FrozenGas };
        let declaration = "thermochemical 300 100000 mole 1 synthetic-alpha 1";
        let text = fixture(reference(100000.0), declaration).replace("initial rest",
            "reservoir_aperture 1 low reference 1\nreservoir_cell 6 0 1\ninitial rest");
        let file = ResolvedDgFilePerformance::from_bytes_with_sources(text.as_bytes(), 8,
            Arc::new(CancelGate::new_clock_free()), Some(&sources),
            |_| Err("no material wall in this fixture".into())).unwrap();
        let original = file.gas_records()[0].state();
        let aperture = file.renderer().system().zone().reservoir_aperture().unwrap();
        assert_eq!(aperture.state.density_kg_m3, original.density);
        assert_eq!(aperture.state.temperature_k, original.temperature);
        densities.push(aperture.state.density_kg_m3);
    }
    assert_ne!(densities[0], densities[1]);
}

fn record(pressure_pa: f64) -> String {
    format!(
        "reservoir_aperture 1 low {:.17e} 0 0 0 300 3",
        reference(pressure_pa).density
    )
}

fn input(pressure_pa: f64) -> String {
    fixture(reference(100000.0), "dry-air 300 100000")
        .replace(
            "initial rest",
            &format!(
                "{}\nreservoir_cell 7 0 1\nreservoir_cell 5 0 1\nreservoir_cell 6 0 1\ninitial rest",
                record(pressure_pa)
            ),
        )
        .replace("boundary 1", "boundary pressure 1")
        .replace("hold 32 -0.1", "hold 32 0")
}

#[test]
fn original_reservoir_state_reaches_flow_and_duct_with_callback_identity() {
    let file = {
        let text = input(100002.0);
        read(&text, 8)
    }; // The source text can go away; the grid owns its original declaration.
    let expected = NsfReservoirAperture {
        axis: 1,
        high: false,
        cells: vec![[7, 0, 1], [5, 0, 1], [6, 0, 1]],
        state: NsfOpenState {
            density_kg_m3: reference(100002.0).density,
            velocity_m_s: [0.0; 3],
            temperature_k: 300.0,
        },
    };
    let zone = file.renderer().system().zone();
    assert_eq!(zone.reservoir_aperture(), Some(&expected));
    assert_eq!(zone.last_reservoir_mass_out_kg(), Some(0.0));
    assert_eq!(zone.steps(), 0);
    let (density, velocity, temperature) = zone.cell([6, 0, 1]);
    assert!((density - 1.0).abs() < 2e-14);
    assert!(velocity.iter().all(|v| v.abs() < 2e-14));
    assert_eq!(temperature, 300.0);

    let (whole, complete) = render(file, 8);
    let (parts, partitioned) = render(read(&input(100002.0), 3), 3);
    assert_eq!(whole, parts);
    assert_eq!(
        complete.system().zone().state_hash(),
        partitioned.system().zone().state_hash()
    );
    assert_eq!(
        complete.system().load().state(),
        partitioned.system().load().state()
    );
    assert_eq!(
        complete.system().zone().reservoir_aperture(),
        Some(&expected)
    );
    assert!(whole.iter().any(|p| p.abs() > 1e-8));

    // Only exterior density changes: temperature, mesh, mechanical clock,
    // material law and zero characteristic-face pressure remain identical.
    let (higher, changed) = render(read(&input(100004.0), 8), 8);
    assert_ne!(whole, higher);
    assert_ne!(
        complete.system().zone().cell([6, 0, 1]),
        changed.system().zone().cell([6, 0, 1])
    );
    let sealed = fixture(reference(100000.0), "dry-air 300 100000")
        .replace("boundary 1", "boundary pressure 1")
        .replace("hold 32 -0.1", "hold 32 0");
    let (closed_pressure, closed) = render(read(&sealed, 8), 8);
    assert_ne!(whole, closed_pressure);
    assert!(closed.system().zone().reservoir_aperture().is_none());
    assert_eq!(closed.last_frame().unwrap().reservoir_mass_out_kg, 0.0);
}

#[test]
fn accepted_native_frames_close_mass_and_energy_with_the_separate_reservoir() {
    let mut renderer = read(&input(100002.0), 1).into_renderer();
    let mass_before = renderer.system().zone().total_mass_kg();
    let energy_before = renderer.system().zone().total_energy_j();
    let mut reservoir_mass = 0.0;
    let mut characteristic_mass = 0.0;
    let mut reservoir_energy = 0.0;
    for step in 1..=SAMPLES {
        let mut pressure = [0.0];
        renderer.block(&mut pressure).unwrap();
        assert!(pressure[0].is_finite());
        let frame = renderer.last_frame().unwrap();
        assert_eq!(frame.step, step as u64);
        assert!(frame.reservoir_mass_out_kg.is_finite());
        assert!(frame.flow.open_energy_j.is_finite());
        assert_eq!(
            Some(frame.reservoir_mass_out_kg),
            renderer.system().zone().last_reservoir_mass_out_kg()
        );
        let faces = frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
        let mass_error = frame.flow_mass_change_kg + faces + frame.reservoir_mass_out_kg;
        assert!(mass_error.abs() / mass_before < 2e-11, "{frame:?}");
        let energy_error = frame.flow_energy_change_j
            - frame.flow.open_energy_j
            - frame.flow.open_boundary_energy_j
            - frame.flow.wall_heat_j
            - frame.flow.source_energy_j;
        assert!(energy_error.abs() / energy_before < 2e-11, "{frame:?}");
        assert_eq!(frame.flow.wall_heat_j, 0.0);
        assert_eq!(frame.flow.source_energy_j, 0.0);
        reservoir_mass += frame.reservoir_mass_out_kg;
        characteristic_mass += faces;
        reservoir_energy += frame.flow.open_energy_j;
    }
    let mass_change = renderer.system().zone().total_mass_kg() - mass_before;
    assert!((mass_change + characteristic_mass + reservoir_mass).abs() / mass_before < 2e-11);
    assert!(reservoir_mass.abs() > 1e-12 * mass_before);
    assert!(reservoir_energy.abs() > 1e-12 * energy_before);
    // Omitting the side aperture is an observable accounting error; it cannot
    // be hidden in either of the two characteristic-face mass accounts.
    assert!((mass_change + characteristic_mass).abs() > 1e-12 * mass_before);
}

#[test]
fn native_reservoir_refuses_malformed_states_and_ambiguous_wall_topology() {
    let text = input(100002.0);
    let source = record(100002.0);
    let rho = reference(100002.0).density;
    let mut invalid = vec![
        text.replace(&source, "reservoir_aperture 1 low 1 0 0 0 300"),
        text.replace(&source, &format!("{source} extra")),
        text.replace(&source, "reservoir_aperture 1 low NaN 0 0 0 300 3"),
        text.replace(&source, "reservoir_aperture 1 low -1 0 0 0 300 3"),
        text.replace(
            &source,
            &format!("reservoir_aperture 1 low {rho:.17e} 0 0 0 700 3"),
        ),
        text.replace(
            &source,
            &format!("reservoir_aperture 1 low {rho:.17e} 10000 0 0 300 3"),
        ),
        text.replace("reservoir_aperture 1 low", "reservoir_aperture 1 middle"),
        text.replace("reservoir_aperture 1 low", "reservoir_aperture 3 low"),
        text.replace("reservoir_aperture 1 low", "reservoir_aperture 0 low"),
        text.replace("300 3\nreservoir_cell", "300 0\nreservoir_cell"),
        text.replace("300 3\nreservoir_cell", "300 109\nreservoir_cell"),
        text.replace("reservoir_cell 7 0 1", "reservoir_cell 5 0 1"),
        text.replace("solids 0", "solids 1\nsolid 6 1 1"),
        text.replace("grid 12 3 3", "grid 12 1 3"),
    ];
    // Out-of-grid, wrong face, characteristic collar, and transverse corner.
    for cell in ["12 0 1", "6 1 1", "1 0 1", "6 0 0"] {
        invalid.push(text.replace("reservoir_cell 7 0 1", &format!("reservoir_cell {cell}")));
    }
    for (case, source) in invalid.iter().enumerate() {
        assert!(
            ResolvedDgFilePerformance::from_bytes(
                source.as_bytes(),
                8,
                Arc::new(CancelGate::new_clock_free())
            )
            .is_err(),
            "invalid reservoir case {case} was admitted"
        );
    }
}
