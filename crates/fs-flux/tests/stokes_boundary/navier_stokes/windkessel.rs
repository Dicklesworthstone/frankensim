//! G1/G3: live laminar convection, constitutive updates and dynamic port laws.

use super::*;
use fs_flux::{StokesWindkesselHistory, StokesWindkesselPort};

#[test]
fn nonlinear_windkessel_resolves_convection_rheology_and_relaxed_pressure_together() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let (dt, old_dt) = (0.03, 0.02);
    let base = |p: [f64; 2]| [0.4 + 0.3 * p[1], 0.1];
    let viscosity = |rate: f64| 0.2 + 1.2 / (1.0 + (2.0 * rate).powi(2)).powf(0.25);
    for scale in [-1.0, 1.0] {
        let (recent_scale, older_scale) = (0.7 * scale, 0.4 * scale);
        let target = |p| base(p).map(|u| scale * u);
        let previous = moments(&mesh, &|p| base(p).map(|u| recent_scale * u));
        let older = moments(&mesh, &|p| base(p).map(|u| older_scale * u));
        let snapshots = (previous.clone(), older.clone());
        let (xp, xo) = (0.3, 0.4);
        let x_previous = [xp];
        let x_older = [xo];
        for bdf2 in [false, true] {
            let (a0, a1, a2) = if bdf2 {
                let r = dt / old_dt;
                ((1.0 + 2.0 * r) / (1.0 + r), -(1.0 + r), r * r / (1.0 + r))
            } else {
                (1.0, -1.0, 0.0)
            };
            let derivative_scale = (a0 * scale + a1 * recent_scale + a2 * older_scale) / dt;
            let q = 0.55 * scale;
            let (rp, rd, c, inertia, added_pressure) = (0.3, 2.1, 0.6, 0.4, 1.4);
            let x = (q - c * (a1 * xp + a2 * xo) / dt) / (c * a0 / dt + 1.0 / rd);
            let reservoir = added_pressure - rp * q - inertia * 0.55 * derivative_scale - x;
            let load = StokesWindkesselPort {
                reservoir_pressure_m2_s2: reservoir,
                proximal_resistance_s_inv: rp,
                distal_resistance_s_inv: rd,
                compliance_s2: c,
                inertance: inertia,
            };
            let mu = viscosity(0.3 * scale.abs());
            let actual_traction = |p: [f64; 2], n: [f64; 2]| {
                let pressure = 7.0 + 2.0 * p[0] - 3.0 * p[1];
                [
                    -pressure * n[0] + 0.3 * scale * mu * n[1],
                    0.3 * scale * mu * n[0] - pressure * n[1],
                ]
            };
            let original_traction = |p: [f64; 2], n: [f64; 2]| {
                let mut t = actual_traction(p, n);
                if n[0] > 0.5 {
                    for i in 0..2 {
                        t[i] += added_pressure * n[i];
                    }
                }
                t
            };
            let mut outlet = Vec::new();
            let boundaries: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, geometry)| {
                    if geometry.normal[0] > 0.5 {
                        outlet.push(edge);
                    }
                    if geometry.normal[0].abs() > 0.5 {
                        NavierStokesBoundary::Open {
                            edge,
                            traction: &original_traction,
                            backflow: OpenBoundaryPolicy::NaturalTraction,
                        }
                    } else {
                        NavierStokesBoundary::Velocity {
                            edge,
                            value: &target,
                        }
                    }
                })
                .collect();
            // Nonzero u.grad(u)=(.03*scale²,0); pressure can no longer absorb
            // its removal while preserving both the native datum and RCR law.
            let force = |p| {
                let v = base(p);
                [
                    derivative_scale * v[0] + 0.03 * scale * scale + 2.0,
                    derivative_scale * v[1] - 3.0,
                ]
            };
            let old = bdf2.then_some((
                StokesWindkesselHistory {
                    velocity: &older,
                    compliance_pressures_m2_s2: &x_older,
                },
                old_dt,
            ));
            for relaxation in [1.0, 0.65] {
                let smallest_rate = Cell::new(f64::INFINITY);
                let largest_rate = Cell::new(0.0_f64);
                let law = |_: usize, _: [f64; 2], rate: f64| {
                    smallest_rate.set(smallest_rate.get().min(rate));
                    largest_rate.set(largest_rate.get().max(rate));
                    viscosity(rate)
                };
                let solved = system
                    .bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports(
                        params,
                        &law,
                        &force,
                        &boundaries,
                        &[&outlet],
                        &[load],
                        StokesWindkesselHistory {
                            velocity: &previous,
                            compliance_pressures_m2_s2: &x_previous,
                        },
                        old,
                        dt,
                        NavierStokesControls {
                            maximum_iterations: 160,
                            relative_tolerance: 1e-10,
                            relaxation,
                        },
                        1e-10,
                        1e-11,
                    )
                    .unwrap();
                let endpoint = &solved.flow.solution;
                let (error, divergence) = system.velocity_error(&endpoint.field.x, &target);
                assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
                assert!(solved.flow.iterations > 1 && solved.flow.nonlinear_residual < 1e-10);
                assert!(largest_rate.get() - smallest_rate.get() > 0.08);
                assert!((solved.pressure_offsets_m2_s2[0] - added_pressure).abs() < 4e-8);
                assert!((solved.outward_volume_rates_m2_s[0] - q).abs() < 3e-9);
                assert!((solved.balances[0].compliance_pressure_m2_s2 - x).abs() < 3e-9);
                assert!(solved.balances[0].power_balance_residual_m4_s3.abs() < 4e-8);
                for (pressure, point) in endpoint.pressure.iter().zip(&mesh.centroids) {
                    assert!((pressure - (7.0 + 2.0 * point[0] - 3.0 * point[1])).abs() < 4e-8);
                }
                let true_boundaries: Vec<_> = boundaries
                    .iter()
                    .map(|boundary| match *boundary {
                        NavierStokesBoundary::Open { edge, backflow, .. } => {
                            NavierStokesBoundary::Open {
                                edge,
                                traction: &actual_traction,
                                backflow,
                            }
                        }
                        other => other,
                    })
                    .collect();
                let history: Vec<_> = previous
                    .iter()
                    .zip(&older)
                    .map(|(recent, old)| -(a1 * recent + a2 * old) / a0)
                    .collect();
                let assembled = system
                    .assemble_navier_stokes_mixed_bdf1_variable_viscosity(
                        params,
                        &|_, _| mu,
                        &force,
                        &true_boundaries,
                        &endpoint.field.x[..system.n_u],
                        &history,
                        dt / a0,
                    )
                    .unwrap();
                let mut residual = vec![0.0; system.n];
                assembled.matrix().spmv(&endpoint.field.x, &mut residual);
                let norm = residual
                    .iter()
                    .zip(assembled.rhs())
                    .fold(0.0_f64, |norm, (ax, b)| norm.hypot(ax - b));
                assert!(norm < 5e-8, "independent nonlinear residual {norm}");
            }
            assert_eq!(previous, snapshots.0);
            assert_eq!(older, snapshots.1);
            assert!(matches!(
                system.bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports(
                    params,
                    &|_, _, rate| viscosity(rate),
                    &force,
                    &boundaries,
                    &[&outlet],
                    &[load],
                    StokesWindkesselHistory {
                        velocity: &previous,
                        compliance_pressures_m2_s2: &x_previous
                    },
                    old,
                    dt,
                    NavierStokesControls {
                        maximum_iterations: 1,
                        ..NavierStokesControls::default()
                    },
                    1e-10,
                    1e-11,
                ),
                Err(StokesBoundaryError::Rheology(
                    fs_flux::FluxRheologyError::NonlinearConvergence { .. }
                ))
            ));
        }
    }
}
