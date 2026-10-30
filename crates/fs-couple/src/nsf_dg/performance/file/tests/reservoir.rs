//! Native side-reservoir records retain the original ghost and flux owners.
//! The small pressure differences below are numerical integration fixtures,
//! not a qualified flute aperture, jet or radiation boundary.
use super::*;

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
