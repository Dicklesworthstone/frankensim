//! G1: nonlinear compact links and ownership of collapsed component sides.

use super::*;

fn potential(t: f64) -> f64 {
    let x = t - 300.0;
    x + 0.005 * x * x
}

fn from_potential(p: f64) -> f64 {
    300.0 + 2.0 * p / (1.0 + (1.0 + 0.02 * p).sqrt())
}

fn branch_reference(junction: f64, reservoir: f64, width: f64, area: f64, r: f64) -> (f64, f64) {
    // Q=(Tj-Tsurface)/R and integral(k dT)=width*Q/area. Eliminating Q
    // gives 0.005*dT²+(k_reservoir+beta)*dT-beta*(Tj-Treservoir)=0.
    // This whole-slab quadratic does not use the production half-cell law.
    let beta = width / (area * r);
    let linear = 1.0 + 0.01 * (reservoir - 300.0) + beta;
    let drive = beta * (junction - reservoir);
    let drop = 2.0 * drive / (linear + (linear * linear + 0.02 * drive).sqrt());
    let surface = reservoir + drop;
    let heat = (junction - surface) / r;
    let cell = from_potential(0.5 * (potential(surface) + potential(reservoir)));
    (heat, cell)
}

#[test]
fn nonlinear_compact_links_match_whole_slab_oracles_and_area_shares() {
    // Two unequal face cells make each K/W junction resistor split 1:2 in
    // parallel. The package thickness and its tiny k do not add a bulk path.
    let domain =
        VoxelDomain::graded_from_fn([vec![0.2, 0.4], vec![0.3], vec![0.5, 0.2, 0.8]], |p| {
            Voxel::Solid(u16::from(p[2] > 0.5 && p[2] < 0.7))
        })
        .unwrap();
    let solids = [
        SolidMaterial::new("nonlinear outside slabs", 1.0).with_conductivity_table(&[
            (250.0, 0.5),
            (290.0, 0.9),
            (330.0, 1.3),
            (370.0, 1.7),
            (450.0, 2.5),
        ]),
        SolidMaterial::new("collapsed package", 1e-6),
    ];
    let area = (0.2 + 0.4) * 0.3;
    let (board_r, case_r) = (2.0, 5.0);
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[4] = ThermalFace::Temperature(300.0);
    faces[5] = ThermalFace::Temperature(310.0);
    for target in [280.0, 400.0] {
        let (board_w, board_cell) = branch_reference(target, 300.0, 0.5, area, board_r);
        let (case_w, case_cell) = branch_reference(target, 310.0, 0.8, area, case_r);
        let power = board_w + case_w;
        let mut setup = ThermalSetup::new(faces);
        setup.compact_components.push(CompactComponent {
            lo: [0, 0, 1],
            hi: [2, 1, 2],
            board_face: Face3::ZMin,
            power_w: power,
            junction_to_case_k_w: case_r,
            junction_to_board_k_w: board_r,
            heat_capacity_j_k: None,
        });
        let solution = solve_energy(
            &domain,
            &unit_fluid(),
            &solids,
            &FlowField::quiescent(&domain),
            &setup,
            &EnergyConfig::default(),
            &CancelGate::new(),
        )
        .unwrap();
        let junction = &solution.junctions[0];
        assert!(
            (junction.temperature_k - target).abs() < 1e-7,
            "{junction:?}"
        );
        assert!((junction.board_w - board_w).abs() < 1e-7);
        assert!((junction.case_w - case_w).abs() < 1e-7);
        assert!((junction.board_w + junction.case_w - power).abs() < 1e-7);
        assert_eq!(junction.links.len(), 4);
        for &(cell, conductance) in &junction.links {
            let (heat, temperature) = if domain.coords(cell)[2] == 0 {
                (board_w, board_cell)
            } else {
                (case_w, case_cell)
            };
            assert!((solution.temperature[cell] - temperature).abs() < 1e-7);
            let actual = conductance * (junction.temperature_k - solution.temperature[cell]);
            let share = heat * domain.face_area(cell, 2) / area;
            assert!(
                (actual - share).abs() < 1e-7,
                "target={target}, cell={cell}"
            );
        }
        for x in 0..2 {
            assert!((solution.temperature[domain.index(x, 0, 1)] - target).abs() < 1e-7);
        }
        let balance = solution.report.balance;
        assert!((balance.source_w - power).abs() < 1e-12);
        assert!((balance.boundary_outflow_w - power).abs() < 1e-7);
        assert!(balance.relative_residual < 1e-9, "{balance:?}");
    }
}

#[test]
fn collapsed_component_sides_do_not_invent_ordinary_boundary_heat() {
    let domain = VoxelDomain::from_fn(1, 1, 3, 1.0, |_| Voxel::Solid(0)).unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::HeatFlux(5.0);
    let mut setup = ThermalSetup::new(faces);
    setup.fixed_temperature = vec![(0, 300.0), (2, 300.0)];
    setup.compact_components.push(CompactComponent {
        lo: [0, 0, 1],
        hi: [1, 1, 2],
        board_face: Face3::ZMin,
        power_w: 4.0,
        junction_to_case_k_w: 1.0,
        junction_to_board_k_w: 1.0,
        heat_capacity_j_k: None,
    });
    let solution = solve_energy(
        &domain,
        &unit_fluid(),
        &[SolidMaterial::new("solid", 1.0)],
        &FlowField::quiescent(&domain),
        &setup,
        &EnergyConfig::default(),
        &CancelGate::new(),
    )
    .unwrap();
    // Each link is 1 K/W plus a 0.5 K/W outside half cell. Two links carry
    // 2 W each at Tj=303 K. Compact sides are owned as adiabatic; the global
    // 5 W/m² ordinary face condition must not appear on the collapsed cell.
    let junction = &solution.junctions[0];
    assert!((junction.temperature_k - 303.0).abs() < 1e-9);
    assert!((junction.board_w - 2.0).abs() < 1e-9);
    assert!((junction.case_w - 2.0).abs() < 1e-9);
    let balance = solution.report.balance;
    assert!((balance.fixed_cell_injection_w + 4.0).abs() < 1e-9);
    assert_eq!(balance.source_w, 4.0);
    assert_eq!(balance.boundary_outflow_w, 0.0);
    assert!(balance.residual_w.abs() < 1e-9, "{balance:?}");
}
