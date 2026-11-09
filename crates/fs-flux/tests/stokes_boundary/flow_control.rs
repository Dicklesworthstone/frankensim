//! Integral flow constraints leave the boundary profile to the original solver.

use super::*;
use fs_flux::StokesResistancePort;

#[path = "flow_control/windkessel.rs"]
mod windkessel;

/// G1: independently integrate plug-flow momentum and its pressure work. The
/// channel has perfect-slip walls; dissipation is solely the external resistor.
#[test]
fn resistance_ports_close_transient_momentum_and_external_power_under_reversal() {
    let (length, height, dt) = (2.0, 1.5, 0.1);
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(length, height, 2, 2));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let mut inlet = Vec::new();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| {
            let n = mesh.edges[edge].normal;
            if n[0] < -0.5 {
                inlet.push(edge);
            }
            if n[1].abs() > 0.5 {
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
    let prepared = system
        .assemble_stokes_mixed_bdf1(
            FluxParams::default(),
            &zero,
            &boundaries,
            &vec![0.0; system.n_u],
            dt,
        )
        .unwrap()
        .prepare_flow_control(&[&inlet])
        .unwrap();
    for resistance in [0.0, 5.0, 200.0] {
        for reservoir in [-2.0, 0.0, 2.0] {
            let load = StokesResistancePort {
                reservoir_pressure_m2_s2: reservoir,
                resistance_s_inv: resistance,
            };
            let solved = prepared.solve_resistance_ports(&[load], 1e-12).unwrap();
            let speed = reservoir / (length / dt + resistance * height);
            let (error, divergence) =
                system.velocity_error(&solved.flow.field.x, &|_| [speed, 0.0]);
            assert!(error < 2e-10 && divergence < 2e-10, "{error}, {divergence}");
            let q = -height * speed;
            assert!((solved.outward_volume_rates_m2_s[0] - q).abs() < 1e-11);
            assert!((solved.pressure_offsets_m2_s2[0] - length * speed / dt).abs() < 1e-10);
            let loss = resistance * q * q;
            assert!((solved.resistance_dissipation_m4_s3[0] - loss).abs() < 1e-11);
            let kinetic_rate = 0.5 * length * height * speed * speed / dt;
            // Backward Euler's numerical damping equals the kinetic gain for
            // this start from rest. Neither term belongs in physical R*Q².
            assert!((-reservoir * q - 2.0 * kinetic_rate - loss).abs() < 2e-11);
            let reaction = solved
                .flow
                .reactions
                .resultant(&inlet, [0.0; 2], &[None])
                .unwrap();
            assert!((reaction.imposed.force[0] + height * length * speed / dt).abs() < 2e-10);
        }
    }
}

/// G0/G3: an independent full-field rank-one Robin assembly checks all
/// cross-port terms, material variation, boundary work and cache reuse.
#[test]
fn resistance_ports_match_full_robin_operator_and_interact_across_ports() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system = FluxSystem::new(&mesh);
    let queries = Cell::new(0);
    let viscosity = |_: usize, point: [f64; 2]| {
        queries.set(queries.get() + 1);
        0.5 + 0.4 * point[1]
    };
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let mut ports = [Vec::new(), Vec::new()];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| {
            let g = &mesh.edges[edge];
            if g.normal[0] < -0.5 {
                ports[usize::from(g.mid[1] > 0.5)].push(edge);
            }
            if g.normal[1].abs() > 0.5 {
                StokesBoundary::Velocity { edge, value: &zero }
            } else {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let params = FluxParams::default();
    let assembly = system
        .assemble_stokes_mixed_variable_viscosity(params, &viscosity, &zero, &boundaries)
        .unwrap();
    let original = assembly.matrix().clone();
    let original_rhs = assembly.rhs().to_vec();
    let prepared = assembly
        .prepare_flow_control(&[&ports[0], &ports[1]])
        .unwrap();
    let calls = queries.get();
    let loads = [
        StokesResistancePort {
            reservoir_pressure_m2_s2: 2.0,
            resistance_s_inv: 3.0,
        },
        StokesResistancePort {
            reservoir_pressure_m2_s2: -1.0,
            resistance_s_inv: 7.0,
        },
    ];
    let solved = prepared.solve_resistance_ports(&loads, 1e-12).unwrap();
    let mut matrix = fs_sparse::Coo::new(system.n, system.n);
    let mut rhs = original_rhs;
    for row in 0..system.n {
        let (columns, values) = original.row(row);
        for (&col, &value) in columns.iter().zip(values) {
            matrix.push(row, col, value);
        }
    }
    let mut external_power = 0.0;
    for (port, load) in ports.iter().zip(loads) {
        for &row_edge in port {
            let row_weight = mesh.edges[row_edge].len;
            rhs[2 * row_edge] -= load.reservoir_pressure_m2_s2 * row_weight;
            for &col_edge in port {
                let resistance = load.resistance_s_inv * row_weight * mesh.edges[col_edge].len;
                matrix.push(2 * row_edge, 2 * col_edge, resistance);
                external_power += resistance
                    * solved.flow.field.x[2 * row_edge]
                    * solved.flow.field.x[2 * col_edge];
            }
        }
    }
    let direct = system.solve_linear(&matrix.assemble(), &rhs, params);
    assert!(direct.converged(params.tol));
    for (&native, &oracle) in solved.flow.field.x.iter().zip(&direct.x) {
        assert!((native - oracle).abs() < 2e-9, "{native} versus {oracle}");
    }
    assert!(
        (external_power - solved.resistance_dissipation_m4_s3.iter().sum::<f64>()).abs() < 1e-10
    );
    let mut changed = loads;
    changed[1].reservoir_pressure_m2_s2 = 1.0;
    let second = prepared.solve_resistance_ports(&changed, 1e-12).unwrap();
    assert!(
        (solved.outward_volume_rates_m2_s[0] - second.outward_volume_rates_m2_s[0]).abs() > 1e-4
    );
    assert_eq!(
        queries.get(),
        calls,
        "cached loads do not re-query viscosity"
    );
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        changed[0].resistance_s_inv = invalid;
        assert!(matches!(
            prepared.solve_resistance_ports(&changed, 1e-12),
            Err(StokesBoundaryError::InvalidInput(_))
        ));
    }
    assert!(prepared.solve_resistance_ports(&loads[..1], 1e-12).is_err());
    assert!(prepared.solve_resistance_ports(&loads, f64::NAN).is_err());
    let again = prepared.solve_resistance_ports(&loads, 1e-12).unwrap();
    assert_eq!(
        again.flow.field.x, solved.flow.field.x,
        "refusal never mutates the cache"
    );
}

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
