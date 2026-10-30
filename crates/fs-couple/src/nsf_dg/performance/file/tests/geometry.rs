//! SI mouth geometry reaches the existing fluid/acoustic/render owners.
use super::*;

#[test]
fn physical_source_bore_and_reservoir_apertures_replay_the_original_cell_owners() {
    let gas = reference(100_000.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    for high in [false, true] {
        for key in ["source_patch", "coupled_patch", "reservoir_aperture"] {
            let mut base = fixture(gas, "dry-air 300 100000").replace(
                "characteristic 0 high",
                if high {
                    "characteristic 0 high"
                } else {
                    "characteristic 0 low"
                },
            );
            let (indexed, physical) = if key == "reservoir_aperture" {
                let side = if high { "high" } else { "low" };
                let normal = if high { 2 } else { 0 };
                (
                    format!(
                        "reservoir_aperture 1 {side} {:.17e} 0 0 0 300 3\nreservoir_cell 5 {normal} 1\nreservoir_cell 6 {normal} 1\nreservoir_cell 7 {normal} 1\n",
                        gas.density
                    ),
                    format!(
                        "reservoir_aperture_m 1 {side} {:.17e} 0 0 0 300 {:.17e} {:.17e} {:.17e} {:.17e}\n",
                        gas.density,
                        5.0 * dx,
                        8.0 * dx,
                        dx,
                        2.0 * dx
                    ),
                )
            } else {
                if key == "coupled_patch" {
                    let whole_radius = (9.0 * dx * dx / core::f64::consts::PI).sqrt();
                    let patch_radius = (dx * dx / core::f64::consts::PI).sqrt();
                    base = base.replace(
                        &format!("{whole_radius:.17e}"),
                        &format!("{patch_radius:.17e}"),
                    );
                }
                (
                    format!("{key} 1 2 1 2 pressure 0\n"),
                    format!(
                        "{key}_m {:.17e} {:.17e} {:.17e} {:.17e} pressure 0\n",
                        dx,
                        2.0 * dx,
                        dx,
                        2.0 * dx
                    ),
                )
            };
            let a = read(
                &base.replace("initial rest", &(indexed + "initial rest")),
                8,
            );
            let b = read(
                &base.replace("initial rest", &(physical + "initial rest")),
                3,
            );
            assert_eq!(
                a.renderer().system().zone().state_hash(),
                b.renderer().system().zone().state_hash()
            );
            assert_eq!(
                a.renderer().system().zone().characteristic_patch(),
                b.renderer().system().zone().characteristic_patch()
            );
            assert_eq!(
                a.renderer().system().zone().reservoir_aperture(),
                b.renderer().system().zone().reservoir_aperture()
            );
            let (a, owner_a) = render(a, 8);
            let (b, owner_b) = render(b, 3);
            assert_eq!(a, b, "{key}, high={high}");
            assert_eq!(
                owner_a.system().zone().state_hash(),
                owner_b.system().zone().state_hash()
            );
            assert_eq!(
                owner_a.system().load().state(),
                owner_b.system().load().state()
            );
        }
    }
}

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

#[test]
fn curved_si_geometry_replays_the_original_explicit_cell_owners() {
    let gas = reference(100_000.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let baseline = fixture(gas, "dry-air 300 100000");
    let declarations = format!(
        "solid_geometry 3 324\n\
solid_cylinder add {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}\n\
solid_cylinder subtract {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}\n\
solid_sphere add {:.17e} {:.17e} {:.17e} {:.17e}",
        5.0 * dx,
        1.5 * dx,
        1.5 * dx,
        7.0 * dx,
        1.5 * dx,
        1.5 * dx,
        1.15 * dx,
        5.0 * dx,
        1.5 * dx,
        1.5 * dx,
        7.0 * dx,
        1.5 * dx,
        1.5 * dx,
        0.55 * dx,
        5.5 * dx,
        1.5 * dx,
        1.5 * dx,
        0.55 * dx,
    );
    let physical = baseline.replace("solids 0", &declarations);
    let explicit = baseline.replace(
        "solids 0",
        "solids 9\n\
solid 5 1 0\nsolid 5 0 1\nsolid 5 1 1\nsolid 5 2 1\nsolid 5 1 2\n\
solid 6 1 0\nsolid 6 0 1\nsolid 6 2 1\nsolid 6 1 2",
    );
    let physical = read(&physical, 8);
    let explicit = read(&explicit, 3);
    assert_eq!(physical.info().solid_cells, 9);
    assert_eq!(
        physical.renderer().system().zone().state_hash(),
        explicit.renderer().system().zone().state_hash()
    );
    let (a, physical) = render(physical, 8);
    let (b, explicit) = render(explicit, 3);
    assert_eq!(a, b);
    assert_eq!(
        physical.system().zone().state_hash(),
        explicit.system().zone().state_hash()
    );
    assert_eq!(
        physical.system().load().state(),
        explicit.system().load().state()
    );
}
