//! Reservoir laws select pressure when every natural port is controlled.

use super::*;
use fs_flux::{StokesWindkesselHistory, StokesWindkesselPort};

fn plug_moments(mesh: &TriMesh, speed: f64) -> Vec<f64> {
    let mut moments = vec![0.0; 2 * mesh.edges.len()];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        moments[2 * edge] = speed * geometry.normal[0];
    }
    moments
}

#[test]
fn every_channel_port_obeys_reservoir_pressure_conservation_and_power() {
    let (length, height, dt) = (2.0, 1.5, 0.1);
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(length, height, 2, 2));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let datum = 3.0;
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -datum * n);
    let mut ports = [Vec::new(), Vec::new()];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            let n = mesh.edges[edge].normal;
            if n[1].abs() > 0.5 {
                StokesBoundary::NavierSlip {
                    edge,
                    value: &zero,
                    tangential_resistance_m_s: 0.0,
                }
            } else {
                ports[usize::from(n[0] > 0.0)].push(edge);
                StokesBoundary::TotalTraction { edge, value: &traction }
            }
        })
        .collect();
    let previous_q = 0.2;
    let previous = plug_moments(&mesh, previous_q / height);
    let assembly = || {
        system.assemble_stokes_mixed_bdf1(
            FluxParams::default(), &zero, &boundaries, &previous, dt,
        ).unwrap()
    };
    // Prescribed flows still require their independent pressure datum.
    assert!(matches!(
        assembly().prepare_flow_control(&[&ports[0], &ports[1]]),
        Err(StokesBoundaryError::InvalidInput(_)),
    ));
    let prepared = assembly()
        .prepare_pressure_response(&[&ports[0], &ports[1]])
        .unwrap();
    let fluid_inertance = length / height;
    for (left, right) in [(0.0, 0.0), (0.0, 7.0), (3.0, 7.0)] {
        for drive in [-4.0, 0.0, 4.0] {
            let q = (fluid_inertance * previous_q / dt + drive)
                / (fluid_inertance / dt + left + right);
            for shift in [0.0, 17.0] {
                let loads = [
                    StokesResistancePort {
                        reservoir_pressure_m2_s2: shift + drive,
                        resistance_s_inv: left,
                    },
                    StokesResistancePort {
                        reservoir_pressure_m2_s2: shift,
                        resistance_s_inv: right,
                    },
                ];
                let solved = prepared.solve_resistance_ports(&loads, 1e-10).unwrap();
                let (error, divergence) =
                    system.velocity_error(&solved.flow.field.x, &|_| [q / height, 0.0]);
                assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
                assert!((solved.outward_volume_rates_m2_s[0] + q).abs() < 3e-10);
                assert!((solved.outward_volume_rates_m2_s[1] - q).abs() < 3e-10);
                assert!(solved.outward_volume_rates_m2_s.iter().sum::<f64>().abs() < 1e-11);
                let p = [shift + drive - left * q, shift + right * q];
                for (actual, expected) in solved.pressure_offsets_m2_s2.iter().zip(p) {
                    assert!((actual - expected).abs() < 3e-9);
                }
                for (actual, point) in solved.flow.pressure.iter().zip(&mesh.centroids) {
                    let expected = datum + p[0] + (p[1] - p[0]) * point[0] / length;
                    assert!((actual - expected).abs() < 3e-9);
                }
                assert_eq!(solved.flow.pressure_reference, StokesPressureReference::TractionDatum);
                let fluid_change = 0.5 * fluid_inertance * (q * q - previous_q * previous_q) / dt;
                let bdf_loss = 0.5 * fluid_inertance * (q - previous_q).powi(2) / dt;
                let external_loss = solved.resistance_dissipation_m4_s3.iter().sum::<f64>();
                assert!((external_loss - (left + right) * q * q).abs() < 3e-10);
                assert!((drive * q - fluid_change - bdf_loss - external_loss).abs() < 3e-9);
            }
        }
    }
}

#[test]
fn single_zero_flow_port_retains_its_reservoir_pressure_without_mobility_rank() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let mut port = Vec::new();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            if mesh.edges[edge].normal[0] > 0.5 {
                port.push(edge);
                StokesBoundary::TotalTraction { edge, value: &traction }
            } else {
                StokesBoundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let prepared = system
        .assemble_stokes_mixed(FluxParams::default(), &zero, &boundaries)
        .unwrap()
        .prepare_pressure_response(&[&port])
        .unwrap();
    for reservoir in [-2.0, 3.0] {
        for resistance in [0.0, 4.0] {
            let solved = prepared.solve_resistance_ports(
                &[StokesResistancePort {
                    reservoir_pressure_m2_s2: reservoir,
                    resistance_s_inv: resistance,
                }],
                1e-10,
            ).unwrap();
            let (error, divergence) = system.velocity_error(&solved.flow.field.x, &zero);
            assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
            assert!(solved.outward_volume_rates_m2_s[0].abs() < 1e-11);
            for pressure in &solved.flow.pressure {
                assert!((pressure - reservoir).abs() < 3e-9);
            }
        }
    }
}

#[test]
fn both_windkessel_ports_advance_storage_and_close_total_channel_power() {
    let (length, height, dt) = (2.0, 1.5, 0.1);
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(length, height, 2, 2));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let mut ports = [Vec::new(), Vec::new()];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            let n = mesh.edges[edge].normal;
            if n[1].abs() > 0.5 {
                StokesBoundary::NavierSlip {
                    edge,
                    value: &zero,
                    tangential_resistance_m_s: 0.0,
                }
            } else {
                ports[usize::from(n[0] > 0.0)].push(edge);
                StokesBoundary::TotalTraction { edge, value: &traction }
            }
        })
        .collect();
    let (qp, qo) = (0.2, -0.1);
    let previous = plug_moments(&mesh, qp / height);
    let older = plug_moments(&mesh, qo / height);
    let previous_x = [0.3, -0.4];
    let older_x = [0.2, -0.1];
    let loads = [
        StokesWindkesselPort {
            reservoir_pressure_m2_s2: 4.0,
            proximal_resistance_s_inv: 0.3,
            distal_resistance_s_inv: 2.1,
            compliance_s2: 0.6,
            inertance: 0.2,
        },
        StokesWindkesselPort {
            reservoir_pressure_m2_s2: -1.0,
            proximal_resistance_s_inv: 0.5,
            distal_resistance_s_inv: 1.7,
            compliance_s2: 0.8,
            inertance: 0.4,
        },
    ];
    let recent = StokesWindkesselHistory {
        velocity: &previous,
        compliance_pressures_m2_s2: &previous_x,
    };
    for old_dt in [None, Some(0.07)] {
        let (a0, a1, a2) = old_dt.map_or((1.0, -1.0, 0.0), |h| {
            let r = dt / h;
            ((1.0 + 2.0 * r) / (1.0 + r), -(1.0 + r), r * r / (1.0 + r))
        });
        let old = old_dt.map(|h| (
            StokesWindkesselHistory {
                velocity: &older,
                compliance_pressures_m2_s2: &older_x,
            },
            h,
        ));
        let solved = system.bdf2_stokes_mixed_variable_step_windkessel_ports(
            FluxParams::default(), &zero, &boundaries, &[&ports[0], &ports[1]],
            &loads, recent, old, dt, 1e-10, 1e-11,
        ).unwrap();
        let denominator: Vec<_> = loads.iter()
            .map(|load| load.compliance_s2 * a0 / dt + 1.0 / load.distal_resistance_s_inv)
            .collect();
        let offset: Vec<_> = loads.iter().enumerate().map(|(i, load)| {
            -load.compliance_s2 * (a1 * previous_x[i] + a2 * older_x[i])
                / (dt * denominator[i])
        }).collect();
        let fluid_inertance = length / height;
        let total_inertance = fluid_inertance + loads[0].inertance + loads[1].inertance;
        let q = (5.0 - total_inertance * (a1 * qp + a2 * qo) / dt + offset[0] - offset[1])
            / (total_inertance * a0 / dt + 0.8 + 1.0 / denominator[0] + 1.0 / denominator[1]);
        let x = [offset[0] - q / denominator[0], offset[1] + q / denominator[1]];
        let derivative_q = (a0 * q + a1 * qp + a2 * qo) / dt;
        let (error, divergence) = system.velocity_error(&solved.flow.field.x, &|_| [q / height, 0.0]);
        assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
        for (actual, expected) in solved.outward_volume_rates_m2_s.iter().zip([-q, q]) {
            assert!((actual - expected).abs() < 3e-10);
        }
        assert!(solved.outward_volume_rates_m2_s.iter().sum::<f64>().abs() < 1e-11);
        let fluid_change = 0.5 * fluid_inertance * (q * q - qp * qp) / dt;
        let fluid_defect = fluid_inertance * q * derivative_q - fluid_change;
        let mut total_power = fluid_change + fluid_defect;
        for (i, balance) in solved.balances.iter().enumerate() {
            assert!((balance.compliance_pressure_m2_s2 - x[i]).abs() < 3e-10);
            assert!(balance.power_balance_residual_m4_s3.abs() < 3e-9);
            total_power += loads[i].reservoir_pressure_m2_s2 * [-q, q][i]
                + balance.proximal_dissipation_m4_s3
                + balance.distal_dissipation_m4_s3
                + balance.stored_energy_change_rate_m4_s3
                + balance.bdf_energy_defect_m4_s3;
        }
        assert!(total_power.abs() < 3e-9, "total power defect {total_power}");
    }
}
