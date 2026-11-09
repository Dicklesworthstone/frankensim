//! G1/G3: actual Carreau source owners survive coupled convective/storage solves.

use super::*;
use fs_flux::{MaterialStokesWindkesselHistory, MaterialStokesWindkesselPort};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

#[test]
fn carreau_windkessel_preserves_accepted_sources_convection_and_si_power() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let (dt, old_dt) = (0.03, 0.02);
    let target = |p: [f64; 2]| [0.4 + 0.3 * p[1], 0.1];
    let base_moments: Vec<_> = mesh
        .edges
        .iter()
        .flat_map(|edge| {
            let u = target(edge.mid);
            let dy = mesh.verts[edge.verts.1][1] - mesh.verts[edge.verts.0][1];
            [
                u[0] * edge.normal[0] + u[1] * edge.normal[1],
                0.3 * dy * edge.normal[0] / 12.0,
            ]
        })
        .collect();
    let previous: Vec<_> = base_moments.iter().map(|v| 0.7 * v).collect();
    let older: Vec<_> = base_moments.iter().map(|v| 0.4 * v).collect();
    let snapshots = (previous.clone(), older.clone());
    let recent_x = [0.3];
    let old_x = [0.4];
    let history = MaterialStokesWindkesselHistory {
        velocity: &previous,
        compliance_pressures_pa: &recent_x,
    };
    for density in [2.0, 7.0] {
        let reference = transport(density, 1.6);
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-11,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        for power_index in [1.0, 0.4] {
            let law = ThermalCarreauYasuda::new(
                CarreauYasuda::new(1.6, 0.2, 3.0, power_index, 2.0, [0.0, 1e4]).unwrap(),
                TemperatureShift::new(
                    TemperatureShiftModel::Arrhenius {
                        activation_energy_j_mol: 12_000.0,
                    },
                    300.0,
                    [280.0, 340.0],
                )
                .unwrap(),
                ViscosityShiftMode::VerticalAndRate,
            );
            // At the original reference temperature a_T=1. This is the
            // independent declared Carreau equation, not the solver response.
            let mu =
                0.2 + 1.4 * (1.0_f64 + (3.0_f64 * 0.3).powi(2)).powf((power_index - 1.0) / 2.0);
            for bdf2 in [false, true] {
                let (a0, a1, a2) = if bdf2 {
                    let r = dt / old_dt;
                    ((1.0 + 2.0 * r) / (1.0 + r), -(1.0 + r), r * r / (1.0 + r))
                } else {
                    (1.0, -1.0, 0.0)
                };
                let derivative = (a0 + a1 * 0.7 + a2 * 0.4) / dt;
                let old = bdf2.then_some((
                    MaterialStokesWindkesselHistory {
                        velocity: &older,
                        compliance_pressures_pa: &old_x,
                    },
                    old_dt,
                ));
                let (q, rp, rd, compliance, inertance, added_pressure) =
                    (0.55, 0.3, 2.1, 0.6, 0.4, 1.4);
                let x = (q - compliance * (a1 * recent_x[0] + a2 * old_x[0]) / dt)
                    / (compliance * a0 / dt + 1.0 / rd);
                let load = MaterialStokesWindkesselPort {
                    reservoir_pressure_pa: added_pressure
                        - rp * q
                        - inertance * 0.55 * derivative
                        - x,
                    proximal_resistance_pa_s_per_m2: rp,
                    distal_resistance_pa_s_per_m2: rd,
                    compliance_m2_per_pa: compliance,
                    inertance_pa_s2_per_m2: inertance,
                };
                let force = |p| {
                    let u = target(p);
                    [
                        density * (derivative * u[0] + 0.03) + 2.0,
                        density * derivative * u[1] - 3.0,
                    ]
                };
                let traction = |p: [f64; 2], n: [f64; 2]| {
                    let pressure = 7.0 + 2.0 * p[0]
                        - 3.0 * p[1]
                        - if n[0] > 0.5 { added_pressure } else { 0.0 };
                    [
                        -pressure * n[0] + 0.3 * mu * n[1],
                        0.3 * mu * n[0] - pressure * n[1],
                    ]
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
                            MaterialNavierStokesBoundary::Open {
                                edge,
                                traction: &traction,
                                backflow: OpenBoundaryPolicy::NaturalTraction,
                            }
                        } else {
                            MaterialNavierStokesBoundary::Velocity {
                                edge,
                                value: &target,
                            }
                        }
                    })
                    .collect();
                let controls = NavierStokesControls {
                    maximum_iterations: 160,
                    relative_tolerance: 1e-10,
                    relaxation: 0.65,
                };
                let queries = Cell::new(0);
                let smallest_rate = Cell::new(f64::INFINITY);
                let result = owner
                    .bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports(
                        &force,
                        &boundaries,
                        &[&outlet],
                        &[load],
                        history,
                        old,
                        dt,
                        controls,
                        1e-9,
                        1e-10,
                        |point, retained| {
                            assert_eq!(retained, &reference);
                            queries.set(queries.get() + 1);
                            smallest_rate.set(smallest_rate.get().min(point.shear_rate_per_s));
                            law.evaluate(300.0, point.shear_rate_per_s)?.bind(*retained)
                        },
                    )
                    .unwrap();
                let step = &result.flow.step;
                let (error, divergence) =
                    system.velocity_error(&step.solution.discrete().x, &target);
                assert!(error < 4e-9 && divergence < 4e-9, "{error}, {divergence}");
                assert!(step.iterations > 1 && step.nonlinear_residual < 1e-10);
                assert!(queries.get() > mesh.tris.len() && smallest_rate.get() < 0.22);
                let rates = system.shear_rates(&step.solution.discrete().x);
                for (cell, sample) in result.flow.samples().iter().enumerate() {
                    assert_eq!(sample.material().source(), &reference);
                    assert_eq!(sample.point().shear_rate_per_s, rates[cell]);
                    assert_eq!(sample.material().rheology().shear_rate_per_s(), rates[cell]);
                    assert_eq!(
                        sample.material().rheology().temperature().temperature_k(),
                        300.0
                    );
                    near(sample.transport().dynamic_viscosity_pa_s(), mu, 3e-9);
                }
                near(
                    step.solution.viscous_dissipation().total_power_w_per_m(),
                    0.09 * mu,
                    3e-9,
                );
                near(step.solution.kinetic_energy_j_per_m(), 0.16 * density, 3e-9);
                near(result.pressure_offsets_pa[0], added_pressure, 5e-8);
                near(result.outward_volume_rates_m2_s[0], q, 4e-9);
                near(result.balances[0].compliance_pressure_pa, x, 4e-9);
                near(
                    result.balances[0].compliance_energy_j_per_m,
                    0.5 * compliance * x * x,
                    4e-9,
                );
                near(
                    result.balances[0].proximal_dissipation_w_per_m,
                    rp * q * q,
                    4e-9,
                );
                near(
                    result.balances[0].distal_dissipation_w_per_m,
                    x * x / rd,
                    4e-9,
                );
                near(result.balances[0].power_balance_residual_w_per_m, 0.0, 5e-8);
                if power_index == 1.0 {
                    let fixed = owner
                        .bdf2_navier_stokes_mixed_variable_step_windkessel_ports(
                            &force,
                            &boundaries,
                            &[&outlet],
                            &[load],
                            history,
                            old,
                            dt,
                            controls,
                            1e-9,
                            1e-10,
                        )
                        .unwrap();
                    near(fixed.outward_volume_rates_m2_s[0], q, 4e-9);
                    near(fixed.pressure_offsets_pa[0], added_pressure, 5e-8);
                }
                let refused = owner
                    .bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports(
                        &force,
                        &boundaries,
                        &[&outlet],
                        &[load],
                        history,
                        old,
                        dt,
                        controls,
                        1e-9,
                        1e-10,
                        |point, retained| {
                            law.evaluate(400.0, point.shear_rate_per_s)?.bind(*retained)
                        },
                    );
                assert!(matches!(refused, Err(MaterialFluxError::Material(_))));
                assert_eq!(previous, snapshots.0);
                assert_eq!(older, snapshots.1);
                assert_eq!(recent_x, [0.3]);
                assert_eq!(old_x, [0.4]);
            }
        }
    }
}
