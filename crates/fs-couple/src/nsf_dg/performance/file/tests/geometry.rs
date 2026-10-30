//! SI mouth geometry reaches the existing fluid/acoustic/render owners.
use super::*;

#[test]
fn physical_geometry_and_explicit_cells_have_identical_original_trajectories() {
    let gas = reference(100_000.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let baseline = fixture(gas, "dry-air 300 100000");
    let physical = baseline.replace(
        "solids 0",
        &format!(
            "solid_geometry 2 216\n\
solid_box add {:.17e} 0 0 {:.17e} {:.17e} {:.17e}\n\
solid_box subtract {:.17e} {:.17e} 0 {:.17e} {:.17e} {:.17e}",
            5.0 * dx,
            6.0 * dx,
            3.0 * dx,
            3.0 * dx,
            5.0 * dx,
            dx,
            6.0 * dx,
            2.0 * dx,
            3.0 * dx,
        ),
    );
    let explicit = baseline.replace(
        "solids 0",
        "solids 6\n\
solid 5 0 0\nsolid 5 2 0\nsolid 5 0 1\nsolid 5 2 1\nsolid 5 0 2\nsolid 5 2 2",
    );
    let physical = read(&physical, 8);
    let explicit = read(&explicit, 8);
    assert_eq!(physical.info().solid_cells, 6);
    assert_eq!(
        physical.renderer().system().zone().state_hash(),
        explicit.renderer().system().zone().state_hash()
    );
    let (a, physical) = render(physical, 8);
    let (b, explicit) = render(explicit, 8);
    assert_eq!(
        a, b,
        "same original NSF mask, source and DG load must reproduce the same pressure"
    );
    assert_eq!(
        physical.system().zone().state_hash(),
        explicit.system().zone().state_hash()
    );
    assert_eq!(
        physical.system().load().state(),
        explicit.system().load().state()
    );
}

#[test]
fn si_slit_edge_example_admits_original_resolved_geometry_without_a_jet_card() {
    let file = read(
        include_str!("../../../../../examples/resolved-dg-si-slit-edge.performance"),
        8,
    );
    assert_eq!(file.info().dimensions, [24, 12, 5]);
    assert_eq!(file.info().solid_cells, 375);
    assert_eq!(file.renderer().system().zone().law().dx_m, 0.0005);
    assert!(
        file.renderer()
            .system()
            .zone()
            .reservoir_aperture()
            .is_some()
    );
}
