//! G0/G3: independent authored openings replay through the same native owner.
//! These are manufactured composition fixtures, not measured flute anatomy.
use super::*;

fn input() -> String {
    fixture(reference(100000.0), "dry-air 300 100000")
        .replace("initial rest", "reservoir_apertures 2\nreservoir_aperture 1 low reference 2\nreservoir_cell 5 0 1\nreservoir_cell 6 0 1\nreservoir_aperture 1 high reference 3\nreservoir_cell 7 2 1\nreservoir_cell 8 2 1\nreservoir_cell 9 2 1\ninitial rest")
        .replace("boundary 1\nhold 32 -0.1", "boundary aperture 0 3\nhold 8 100002 300\nhold 16 100002 300.01\nhold 32 100000 300")
}

#[test]
fn indexed_breath_and_distinct_escape_keep_geometry_material_and_callback_replay() {
    let text = input();
    let file = read(&text, 8);
    let zone = file.renderer().system().zone();
    assert_eq!(zone.reservoir_aperture_count(), 2);
    assert!(zone.reservoir_aperture().is_none());
    for (index, high, cells) in [(0, false, 2), (1, true, 3)] {
        let source = zone.reservoir_aperture_at(index).unwrap();
        assert_eq!(source.high, high);
        assert_eq!(source.cells.len(), cells);
        assert_eq!(source.state.density_kg_m3, reference(100000.0).density);
        assert!(zone.last_reservoir_flows_at(index).iter().all(|flow| flow.aperture == index));
    }
    for (index, pressure_pa, temperature_k) in [(0, 100002.0, 300.0),
        (8, 100002.0, 300.01), (16, 100000.0, 300.0)]
    {
        assert_eq!(file.renderer().drive_history()[index],
            NsfDgDrive::ReservoirAperture { aperture: 0, pressure_pa, temperature_k });
    }
    let initial_escape = zone.reservoir_aperture_at(1).unwrap().clone();
    let dx = zone.law().dx_m;
    let metric = text
        .replace("reservoir_aperture 1 low reference 2\nreservoir_cell 5 0 1\nreservoir_cell 6 0 1",
            &format!("reservoir_aperture_m 1 low reference {:.17e} {:.17e} {dx:.17e} {:.17e}", 5.0*dx, 7.0*dx, 2.0*dx))
        .replace("reservoir_aperture 1 high reference 3\nreservoir_cell 7 2 1\nreservoir_cell 8 2 1\nreservoir_cell 9 2 1",
            &format!("reservoir_aperture_m 1 high reference {:.17e} {:.17e} {dx:.17e} {:.17e}", 7.0*dx, 10.0*dx, 2.0*dx));
    let (whole, a) = render(file, 8);
    let (parts, b) = render(read(&metric, 3), 3);
    assert_eq!(whole, parts);
    assert_eq!(a.system().zone().state_hash(), b.system().zone().state_hash());
    assert_eq!(a.system().load().state(), b.system().load().state());
    assert_eq!(a.system().zone().reservoir_aperture_at(1), Some(&initial_escape));
    let frame = a.last_frame().unwrap();
    assert_eq!(frame.reservoir_exchanges.len(), 2);
    assert_eq!(frame.reservoir_flows.iter().filter(|flow| flow.aperture == 0).count(), 2);
    assert_eq!(frame.reservoir_flows.iter().filter(|flow| flow.aperture == 1).count(), 3);
    let balance = frame.flow_mass_change_kg + frame.reservoir_mass_out_kg
        + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
    assert!(balance.abs() < 2e-11 * a.system().zone().total_mass_kg());
    let unheated = text.replace("hold 16 100002 300.01", "hold 16 100002 300");
    assert_ne!(whole, render(read(&unheated, 8), 8).0);
}

#[test]
fn single_reservoir_index_zero_keeps_legacy_physical_replay() {
    let text = input().replace("reservoir_apertures 2\n", "")
        .replace("reservoir_aperture 1 high reference 3\nreservoir_cell 7 2 1\nreservoir_cell 8 2 1\nreservoir_cell 9 2 1\n", "");
    let (indexed, a) = render(read(&text, 8), 8);
    let (legacy, b) = render(read(&text.replace("boundary aperture 0", "boundary reservoir"), 3), 3);
    assert_eq!(indexed, legacy);
    assert_eq!(a.system().zone().state_hash(), b.system().zone().state_hash());
    assert_eq!(a.system().load().state(), b.system().load().state());
}

