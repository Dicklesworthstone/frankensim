//! Manufactured slit, free jet gap and solid wedge on the original resolved
//! render path. These are integration checks, not an edge-tone or specimen
//! qualification: the two-cell slit is intentionally a coarse mesh rung.
use super::*;

const SOURCE: &str = include_str!("../../../../../examples/resolved-dg-slit-edge.performance");
const ESCAPE_OBSERVER: &str = "reservoir_observation baffled -0.05125 -0.00725 0\n";
const ESCAPE_FEEDBACK: &str = "reservoir_feedback -5 5 16 32\n";

fn diagnostic_prefix() -> String {
    // The geometry diagnostics retain their original inlet-pressure observation.
    // The canonical file test below exercises the shared exterior and its reaction.
    SOURCE
        .replace(ESCAPE_OBSERVER, "")
        .replace(ESCAPE_FEEDBACK, "")
        .replace("output 2048 20", "output 32 20")
        .replace(
            "observation baffled 0 0 0.2 48 96 2000",
            "observation inlet-pressure",
        )
        .replace("observer_gas stokes-kirchhoff 4096 1000000 50000000\n", "")
        .replace(
            "boundary pressure 4\nhold 64 0.02\nhold 512 0.2\nhold 1792 2\nhold 2048 0",
            "boundary pressure 1\nhold 32 2",
        )
}

#[test]
fn slit_edge_source_binds_the_actual_solid_geometry_and_both_openings() {
    let file = read(SOURCE, 8);
    assert_eq!(file.info().dimensions, [24, 12, 5]);
    assert_eq!(file.info().solid_cells, 375);
    assert_eq!(file.info().initial_state_cells, 0);
    assert!(file.info().physical_intervals > file.info().samples);
    let system = file.renderer().system();
    let zone = system.zone();
    assert_eq!(zone.cell_count(), 1440);
    assert_eq!(zone.law().dx_m, 0.0005);
    assert_eq!(system.limits().substeps, 4);
    let expected_dt = 4.0 * zone.law().dt_s;
    assert!((system.load().propagator().dt() - expected_dt).abs() < 1e-18);
    assert!((system.area_m2() - 30.0 * 0.0005_f64.powi(2)).abs() < 1e-20);
    assert_eq!(
        system.coupled_patch(),
        Some(NsfCharacteristicPatch {
            start: [0, 0],
            end: [6, 5]
        })
    );
    assert_eq!(zone.reservoir_aperture().unwrap().cells.len(), 27);
    let reaction = system.reservoir_radiation().unwrap();
    assert_eq!(reaction.origin_m, [-0.05125, -0.00725]);
    assert_eq!(reaction.cells, [16, 32]);
    assert_eq!(reaction.terminal_radial_rings, 48);
    assert_eq!(reaction.terminal_angular_points, 96);
    assert_eq!(reaction.pressure_bracket_pa, [-5.0, 5.0]);
    assert!(system.load().independent_mouth_response().is_ok());
    assert!(file.renderer().reservoir_receiver().is_some());
    assert_eq!(
        file.renderer().reservoir_observer_gas().unwrap().filters.len(),
        27
    );
    // Source-owned material queries distinguish the actual interior solids
    // from fluid cells. The windway is not a prescribed velocity or jet card.
    for cell in [[2, 4, 2], [7, 7, 2], [14, 6, 2], [21, 8, 2]] {
        assert!(zone.thermochemical_state(cell).is_err(), "solid {cell:?}");
    }
    for cell in [[2, 5, 2], [7, 6, 2], [10, 6, 2], [20, 4, 2], [23, 4, 2]] {
        assert!(zone.thermochemical_state(cell).is_ok(), "fluid {cell:?}");
    }
    assert!(file.renderer().baffled_receiver().is_some());
    assert!(file.renderer().observer_gas().is_some());
    // Construct and advance the same exterior-observed source used by the
    // native player. Its 0.2 m propagation delay makes initial PCM silence
    // physically expected, so this does not assert an invented onset time.
    let mut renderer = file.into_renderer();
    let mut pressure = [0.0; 1];
    renderer.block(&mut pressure).unwrap();
    assert!(pressure[0].is_finite());
    assert!(renderer.system().accepted_steps() > 1);
    let mouth = renderer.last_frame().unwrap().reservoir_radiation.as_ref().unwrap();
    assert_eq!(
        mouth.outward_flow_m3_s,
        mouth.resolved_exchange.normal_volume_flow_m3_s
    );
    assert!(mouth.pressure_residual_pa.abs() < 1e-5);
    assert!((mouth.resolved_exchange.area_m2 - 27.0 * 0.0005_f64.powi(2)).abs() < 1e-20);
    assert_eq!(
        renderer.last_frame().unwrap().inlet_flow_m3_s,
        renderer
            .last_frame()
            .unwrap()
            .coupled_patch
            .unwrap()
            .normal_volume_flow_m3_s
    );
}

#[test]
fn pressure_driven_slit_and_edge_affect_the_coupled_field_with_callback_identity() {
    let source = diagnostic_prefix();
    let (pressure, complete) = render(read(&source, 8), 8);
    let (partitioned_pressure, partitioned) = render(read(&source, 3), 3);
    assert_eq!(pressure, partitioned_pressure);
    assert_eq!(
        complete.system().zone().state_hash(),
        partitioned.system().zone().state_hash()
    );
    assert_eq!(
        complete.system().load().state(),
        partitioned.system().load().state()
    );
    assert!(pressure.iter().any(|p| p.abs() > 1e-10));
    let slit_velocity = complete.system().zone().cell([7, 5, 2]).1;
    assert!(
        slit_velocity[0] > 0.0,
        "the pressure boundary must drive flow through the slit"
    );

    // Remove only the downstream obstacle. The nozzle, pressure, gas,
    // acoustic bore and both exterior controls are identical.
    let no_edge = source
        .lines()
        .filter(|line| {
            let mut fields = line.split_ascii_whitespace();
            fields.next() != Some("solid") || fields.next().unwrap().parse::<usize>().unwrap() < 14
        })
        .collect::<Vec<_>>()
        .join("\n")
        .replace("solids 375", "solids 300")
        + "\n";
    let (without_edge, unobstructed) = render(read(&no_edge, 8), 8);
    assert_ne!(
        pressure, without_edge,
        "the resolved solid edge must reach the bore response"
    );
    assert_ne!(
        complete.system().zone().cell([13, 6, 2]),
        unobstructed.system().zone().cell([13, 6, 2])
    );
    assert!(
        unobstructed
            .system()
            .zone()
            .thermochemical_state([14, 6, 2])
            .is_ok()
    );
}

#[test]
fn slit_edge_intervals_keep_bore_and_escape_mass_and_energy_accounts_separate() {
    let mut renderer = read(&diagnostic_prefix(), 1).into_renderer();
    let mass = renderer.system().zone().total_mass_kg();
    let energy = renderer.system().zone().total_energy_j();
    for _ in 0..8 {
        renderer.block(&mut [0.0]).unwrap();
        let frame = renderer.last_frame().unwrap();
        let mass_error = frame.flow_mass_change_kg
            + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>()
            + frame.reservoir_mass_out_kg;
        assert!(mass_error.abs() / mass < 2e-11, "{frame:?}");
        let energy_error = frame.flow_energy_change_j
            - frame.flow.open_boundary_energy_j
            - frame.flow.open_energy_j
            - frame.flow.wall_heat_j
            - frame.flow.source_energy_j;
        assert!(energy_error.abs() / energy < 2e-11, "{frame:?}");
        assert_eq!(frame.flow.source_energy_j, 0.0);
        assert_eq!(frame.flow.wall_heat_j, 0.0);
        assert_eq!(
            frame.inlet_flow_m3_s,
            frame.coupled_patch.unwrap().normal_volume_flow_m3_s
        );
    }
}
