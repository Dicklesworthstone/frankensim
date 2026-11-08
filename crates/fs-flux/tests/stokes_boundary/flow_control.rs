//! Integral flow constraints leave the boundary profile to the original solver.

use super::*;

#[test]
fn controlled_flow_recovers_graded_shear_pressure_and_physical_reactions() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let samples = Cell::new(0);
    let viscosity = |_: usize, point: [f64; 2]| {
        samples.set(samples.get() + 1);
        0.7 + 0.9 * point[1]
    };
    let velocity = |point: [f64; 2]| [0.2 + 0.6 * point[1], 0.0];
    let datum = 3.0;
    let pressure_drop = 1.8;
    let traction = |point: [f64; 2], normal: [f64; 2]| {
        let shear = (0.7 + 0.9 * point[1]) * 0.6;
        [
            shear * normal[1] - datum * normal[0],
            shear * normal[0] - datum * normal[1],
        ]
    };
    let edges: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .collect();
    let inlet: Vec<_> = edges
        .iter()
        .copied()
        .filter(|&edge| mesh.edges[edge].normal[0] < -0.5)
        .collect();
    let boundaries: Vec<_> = edges
        .iter()
        .map(|&edge| {
            if mesh.edges[edge].normal[1].abs() > 0.5 {
                StokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    // -div(2 nu D(u)) + grad(p) = -0.9*0.6 - pressure_drop.
    let prepared = system
        .assemble_stokes_mixed_variable_viscosity(
            FluxParams::default(),
            &viscosity,
            &|_| [-0.9 * 0.6 - pressure_drop, 0.0],
            &boundaries,
        )
        .unwrap()
        .prepare_flow_control(&[&inlet])
        .unwrap();
    let queried = samples.get();
    let solved = prepared.solve(&[-0.5], 1e-12).unwrap();
    let (error, divergence) = system.velocity_error(&solved.flow.field.x, &velocity);
    assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
    assert!((solved.pressure_offsets_m2_s2[0] - pressure_drop).abs() < 3e-9);
    assert!((solved.outward_volume_rates_m2_s[0] + 0.5).abs() < 1e-12);
    assert_eq!(
        solved.flow.pressure_reference,
        StokesPressureReference::TractionDatum
    );
    for (pressure, point) in solved.flow.pressure.iter().zip(&mesh.centroids) {
        assert!((pressure - (datum + pressure_drop * (1.0 - point[0]))).abs() < 3e-9);
    }
    let inlet_load = solved
        .flow
        .reactions
        .resultant(&inlet, [0.0; 2], &[None])
        .unwrap();
    assert!((inlet_load.stress.force[0] + datum + pressure_drop).abs() < 3e-9);
    assert!((inlet_load.stress.force[1] - 0.6 * (0.7 + 0.45)).abs() < 3e-9);
    assert!((inlet_load.stress.moment - 0.5 * (datum + pressure_drop)).abs() < 3e-9);
    let total = solved
        .flow
        .reactions
        .resultant(&edges, [0.0; 2], &[None])
        .unwrap();
    assert!((total.stress.force[0] + 0.9 * 0.6 + pressure_drop).abs() < 3e-9);
    assert!((total.stress.moment - 0.5 * (0.9 * 0.6 + pressure_drop)).abs() < 3e-9);
    let reduced = prepared.solve(&[-0.25], 1e-12).unwrap();
    assert!((reduced.outward_volume_rates_m2_s[0] + 0.25).abs() < 1e-12);
    assert!(reduced.pressure_offsets_m2_s2[0] < solved.pressure_offsets_m2_s2[0]);
    let stopped = prepared.solve(&[0.0], 1e-12).unwrap();
    assert!(stopped.outward_volume_rates_m2_s[0].abs() < 1e-12);
    assert!(
        system.kinetic_energy(&stopped.flow.field.x) > 0.0,
        "zero net port flow must not erase wall-driven resolved motion"
    );
    assert_eq!(
        samples.get(),
        queried,
        "warm flow control never re-queries viscosity"
    );
}

#[test]
fn controlled_transient_flow_keeps_weighted_inertia_and_separate_pressure_datums() {
    let mesh = separate_domains();
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |point: [f64; 2], normal: [f64; 2]| {
        let datum = if point[0] < 2.0 { 3.0 } else { 7.0 };
        normal.map(|n| -datum * n)
    };
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            if mesh.edges[edge].normal[1].abs() > 0.5 {
                StokesBoundary::NavierSlip {
                    edge,
                    value: &zero,
                    tangential_resistance_m_s: 0.0,
                }
            } else {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let mut ports = [Vec::new(), Vec::new()];
    let mut previous = vec![0.0; system.n_u];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        let component = usize::from(geometry.mid[0] > 2.0);
        previous[2 * edge] = [0.2, -0.1][component] * geometry.normal[0];
        if mesh.is_boundary(edge) && geometry.normal[0] < -0.5 {
            ports[component].push(edge);
        }
    }
    let prepared = system
        .assemble_stokes_mixed_bdf1_variable_density(
            FluxParams::default(),
            &|_, point| if point[0] < 2.0 { 2.0 } else { 5.0 },
            &|_, _| 0.7,
            &zero,
            &boundaries,
            &previous,
            0.1,
        )
        .unwrap()
        .prepare_flow_control(&[&ports[0], &ports[1]])
        .unwrap();
    let solved = prepared.solve(&[-0.7, 0.4], 1e-12).unwrap();
    let exact = |point: [f64; 2]| {
        if point[0] < 2.0 {
            [0.7, 0.0]
        } else {
            [-0.4, 0.0]
        }
    };
    let (error, divergence) = system.velocity_error(&solved.flow.field.x, &exact);
    assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
    for (actual, expected) in solved.pressure_offsets_m2_s2.iter().zip([10.0, -15.0]) {
        assert!((actual - expected).abs() < 3e-9);
    }
    for (pressure, point) in solved.flow.pressure.iter().zip(&mesh.centroids) {
        let expected = if point[0] < 2.0 {
            3.0 + 10.0 * (1.0 - point[0])
        } else {
            7.0 - 15.0 * (4.0 - point[0])
        };
        assert!((pressure - expected).abs() < 3e-9);
    }
    let energy = prepared.kinetic_energy(&solved.flow.field.x).unwrap();
    assert!((energy - 0.5 * (2.0 * 0.7_f64.powi(2) + 5.0 * 0.4_f64.powi(2))).abs() < 3e-10);
}

#[test]
fn redundant_pressure_control_or_a_strong_velocity_port_is_refused() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let mut inlet = Vec::new();
    let mut outlet = Vec::new();
    let mut wall = Vec::new();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            if mesh.edges[edge].normal[1].abs() > 0.5 {
                wall.push(edge);
                StokesBoundary::Velocity { edge, value: &zero }
            } else {
                if mesh.edges[edge].normal[0] < 0.0 {
                    inlet.push(edge);
                } else {
                    outlet.push(edge);
                }
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let assembly = || {
        system
            .assemble_stokes_mixed(FluxParams::default(), &zero, &boundaries)
            .unwrap()
    };
    assert!(matches!(
        assembly().prepare_flow_control(&[&inlet, &outlet]),
        Err(StokesBoundaryError::InvalidInput(_))
    ));
    assert!(matches!(
        assembly().prepare_flow_control(&[&wall]),
        Err(StokesBoundaryError::InvalidInput(_))
    ));
    assert!(matches!(
        assembly().prepare_flow_control(&[&inlet, &inlet]),
        Err(StokesBoundaryError::DuplicateEdge { .. })
    ));
}
