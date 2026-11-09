//! G1: resolved continuum inertia and dynamic hydraulic storage at one endpoint.

use super::*;
use fs_flux::{StokesWindkesselHistory, StokesWindkesselPort};

fn plug_moments(mesh: &TriMesh, speed: f64) -> Vec<f64> {
    let mut moments = vec![0.0; 2 * mesh.edges.len()];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        moments[2 * edge] = speed * geometry.normal[0];
    }
    moments
}

/// Independently integrate the perfect-slip channel: its fluid inertance is
/// length/height and Q=height*u. The original coupled two-state ODE includes
/// that inertia PLUS the external series inertance, with no imposed profile.
#[test]
fn resolved_rcr_inertance_matches_both_bdf_histories_and_total_energy() {
    let (length, height, dt) = (2.0, 1.5, 0.1);
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(length, height, 2, 2));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let datum = 3.0;
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -datum * n);
    let mut outlet = Vec::new();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| {
            let n = mesh.edges[edge].normal;
            if n[0] > 0.5 {
                outlet.push(edge);
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
    let (qp, xp, qo, xo) = (0.2, -0.3, -0.1, 0.4);
    let previous = plug_moments(&mesh, qp / height);
    let older = plug_moments(&mesh, qo / height);
    let snapshots = (previous.clone(), older.clone());
    let previous_x = [xp];
    let older_x = [xo];
    let recent = StokesWindkesselHistory {
        velocity: &previous,
        compliance_pressures_m2_s2: &previous_x,
    };
    for old_dt in [None, Some(0.1), Some(0.07), Some(0.2)] {
        let (a0, a1, a2) = old_dt.map_or((1.0, -1.0, 0.0), |h| {
            let r = dt / h;
            ((1.0 + 2.0 * r) / (1.0 + r), -(1.0 + r), r * r / (1.0 + r))
        });
        let old = old_dt.map(|h| {
            (
                StokesWindkesselHistory {
                    velocity: &older,
                    compliance_pressures_m2_s2: &older_x,
                },
                h,
            )
        });
        for reservoir in [-4.0, 0.0, 4.0] {
            for inertance in [0.0, 0.4] {
                let load = StokesWindkesselPort {
                    reservoir_pressure_m2_s2: reservoir,
                    proximal_resistance_s_inv: 0.3,
                    distal_resistance_s_inv: 2.1,
                    compliance_s2: 0.6,
                    inertance,
                };
                let solved = system
                    .bdf2_stokes_mixed_variable_step_windkessel_ports(
                        FluxParams::default(),
                        &zero,
                        &boundaries,
                        &[&outlet],
                        &[load],
                        recent,
                        old,
                        dt,
                        1e-10,
                        1e-11,
                    )
                    .unwrap();
                let fluid_inertance = length / height;
                let total_inertance = fluid_inertance + inertance;
                let a = total_inertance * a0 / dt + load.proximal_resistance_s_inv;
                let b = load.compliance_s2 * a0 / dt + 1.0 / load.distal_resistance_s_inv;
                let rhs1 = -reservoir - total_inertance * (a1 * qp + a2 * qo) / dt;
                let rhs2 = -load.compliance_s2 * (a1 * xp + a2 * xo) / dt;
                let q = (rhs1 * b - rhs2) / (a * b + 1.0);
                let x = (rhs1 + a * rhs2) / (a * b + 1.0);
                let derivative_q = (a0 * q + a1 * qp + a2 * qo) / dt;
                let p = -fluid_inertance * derivative_q;
                let (error, divergence) =
                    system.velocity_error(&solved.flow.field.x, &|_| [q / height, 0.0]);
                assert!(error < 2e-10 && divergence < 2e-10, "{error}, {divergence}");
                assert!((solved.outward_volume_rates_m2_s[0] - q).abs() < 2e-10);
                assert!((solved.pressure_offsets_m2_s2[0] - p).abs() < 2e-9);
                for (actual, point) in solved.flow.pressure.iter().zip(&mesh.centroids) {
                    assert!((actual - datum - p * point[0] / length).abs() < 2e-9);
                }
                let balance = solved.balances[0];
                assert!((balance.compliance_pressure_m2_s2 - x).abs() < 2e-10);
                assert!(
                    (solved.kinetic_energy_m4_s2 - 0.5 * fluid_inertance * q * q).abs() < 2e-10
                );
                assert!(balance.power_balance_residual_m4_s3.abs() < 2e-9);
                let fluid_change = 0.5 * fluid_inertance * (q * q - qp * qp) / dt;
                let fluid_bdf_defect = fluid_inertance * q * derivative_q - fluid_change;
                let full_balance = reservoir * q
                    + balance.proximal_dissipation_m4_s3
                    + balance.distal_dissipation_m4_s3
                    + fluid_change
                    + balance.stored_energy_change_rate_m4_s3
                    + fluid_bdf_defect
                    + balance.bdf_energy_defect_m4_s3;
                assert!(
                    full_balance.abs() < 2e-9,
                    "total power defect {full_balance}"
                );
                assert_eq!(previous, snapshots.0);
                assert_eq!(older, snapshots.1);
            }
        }
    }
    let load = StokesWindkesselPort {
        reservoir_pressure_m2_s2: 0.0,
        proximal_resistance_s_inv: 0.0,
        distal_resistance_s_inv: 1.0,
        compliance_s2: 1.0,
        inertance: 0.0,
    };
    let incomplete = StokesWindkesselHistory {
        velocity: &previous,
        compliance_pressures_m2_s2: &[],
    };
    assert!(
        system
            .bdf2_stokes_mixed_variable_step_windkessel_ports(
                FluxParams::default(),
                &zero,
                &boundaries,
                &[&outlet],
                &[load],
                incomplete,
                None,
                dt,
                1e-10,
                1e-11,
            )
            .is_err()
    );
    assert!(
        system
            .bdf2_stokes_mixed_variable_step_windkessel_ports(
                FluxParams::default(),
                &zero,
                &boundaries,
                &[&outlet],
                &[load],
                recent,
                Some((recent, dt / 3.0)),
                dt,
                1e-10,
                1e-11,
            )
            .is_err()
    );
}
