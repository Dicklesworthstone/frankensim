//! G1/G3: physical dynamic storage uses the original material reference exactly once.

use super::*;
use fs_flux::{MaterialStokesWindkesselHistory, MaterialStokesWindkesselPort};

#[test]
fn windkessel_keeps_si_inertia_compliance_energy_and_reference_invariance() {
    let (length, height, dt, old_dt) = (2.0, 1.5, 0.1, 0.07);
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(length, height, 2, 2));
    let system = FluxSystem::new(&mesh);
    let mut outlet = Vec::new();
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -3.0 * n);
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| {
            let n = mesh.edges[edge].normal;
            if n[0] > 0.5 {
                outlet.push(edge);
            }
            if n[1].abs() > 0.5 {
                MaterialStokesBoundary::NavierSlip {
                    edge,
                    value: &no_slip,
                    tangential_resistance_pa_s_per_m: 0.0,
                }
            } else {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let (qp, xp, qo, xo) = (0.2, -0.3, -0.1, 0.4);
    let previous = flow_control_history(&mesh, &|_| [qp / height, 0.0]);
    let older = flow_control_history(&mesh, &|_| [qo / height, 0.0]);
    let snapshot = (previous.clone(), older.clone());
    let previous_x = [xp];
    let older_x = [xo];
    let recent = MaterialStokesWindkesselHistory {
        velocity: &previous,
        compliance_pressures_pa: &previous_x,
    };
    let old = MaterialStokesWindkesselHistory {
        velocity: &older,
        compliance_pressures_pa: &older_x,
    };
    for density in [2.0, 7.0] {
        let source = Counted {
            material: transport(density, 0.4),
            queries: Cell::new(0),
        };
        let uniform =
            MaterialFluxSystem::new(&system, &source, MaterialFluxNumerics::default()).unwrap();
        for reference_density in [0.25, 11.0] {
            let reference = transport(reference_density, 0.9);
            let binding =
                MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default())
                    .unwrap();
            let field = binding
                .bind_stokes_materials(vec![&source; mesh.tris.len()])
                .unwrap();
            let calls = source.queries.get();
            for history in [None, Some((old, old_dt))] {
                let (a0, a1, a2) = if history.is_some() {
                    let r = dt / old_dt;
                    ((1.0 + 2.0 * r) / (1.0 + r), -(1.0 + r), r * r / (1.0 + r))
                } else {
                    (1.0, -1.0, 0.0)
                };
                for reservoir in [-4.0, 4.0] {
                    let load = MaterialStokesWindkesselPort {
                        reservoir_pressure_pa: reservoir,
                        proximal_resistance_pa_s_per_m2: 0.3,
                        distal_resistance_pa_s_per_m2: 2.1,
                        compliance_m2_per_pa: 0.6,
                        inertance_pa_s2_per_m2: 0.4,
                    };
                    let homogeneous = uniform
                        .bdf2_stokes_mixed_variable_step_windkessel_ports(
                            &zero_force,
                            &boundaries,
                            &[&outlet],
                            &[load],
                            recent,
                            history,
                            dt,
                            1e-9,
                            1e-10,
                        )
                        .unwrap();
                    let layered = field
                        .bdf2_stokes_mixed_variable_step_windkessel_ports(
                            &zero_force,
                            &boundaries,
                            &[&outlet],
                            &[load],
                            recent,
                            history,
                            dt,
                            1e-9,
                            1e-10,
                        )
                        .unwrap();
                    let fluid_inertance = density * length / height;
                    let total_inertance = fluid_inertance + load.inertance_pa_s2_per_m2;
                    let a = total_inertance * a0 / dt + load.proximal_resistance_pa_s_per_m2;
                    let b = load.compliance_m2_per_pa * a0 / dt
                        + 1.0 / load.distal_resistance_pa_s_per_m2;
                    let rhs1 = -reservoir - total_inertance * (a1 * qp + a2 * qo) / dt;
                    let rhs2 = -load.compliance_m2_per_pa * (a1 * xp + a2 * xo) / dt;
                    let q = (rhs1 * b - rhs2) / (a * b + 1.0);
                    let x = (rhs1 + a * rhs2) / (a * b + 1.0);
                    let p = -fluid_inertance * (a0 * q + a1 * qp + a2 * qo) / dt;
                    for (actual_q, actual_p, energy, heat, balance) in [
                        (
                            homogeneous.outward_volume_rates_m2_s[0],
                            homogeneous.pressure_offsets_pa[0],
                            homogeneous.flow.kinetic_energy_j_per_m(),
                            homogeneous.flow.viscous_dissipation().total_power_w_per_m(),
                            homogeneous.balances[0],
                        ),
                        (
                            layered.outward_volume_rates_m2_s[0],
                            layered.pressure_offsets_pa[0],
                            layered.flow.kinetic_energy_j_per_m(),
                            layered.flow.viscous_dissipation().total_power_w_per_m(),
                            layered.balances[0],
                        ),
                    ] {
                        near(actual_q, q, 3e-9);
                        near(actual_p, p, 3e-8);
                        near(energy, 0.5 * fluid_inertance * q * q, 3e-9);
                        near(heat, 0.0, 1e-15);
                        near(balance.compliance_pressure_pa, x, 3e-9);
                        near(
                            balance.compliance_energy_j_per_m,
                            0.5 * load.compliance_m2_per_pa * x * x,
                            3e-9,
                        );
                        near(
                            balance.inertance_energy_j_per_m,
                            0.5 * load.inertance_pa_s2_per_m2 * q * q,
                            3e-9,
                        );
                        near(
                            balance.proximal_dissipation_w_per_m,
                            load.proximal_resistance_pa_s_per_m2 * q * q,
                            3e-9,
                        );
                        near(
                            balance.distal_dissipation_w_per_m,
                            x * x / load.distal_resistance_pa_s_per_m2,
                            3e-9,
                        );
                        near(balance.port_power_w_per_m, p * q, 3e-8);
                        near(balance.power_balance_residual_w_per_m, 0.0, 3e-8);
                    }
                    for cell in 0..mesh.tris.len() {
                        assert!(std::ptr::eq(
                            &field.samples()[cell],
                            &layered.flow.samples()[cell]
                        ));
                    }
                }
            }
            let bad = MaterialStokesWindkesselPort {
                reservoir_pressure_pa: 0.0,
                proximal_resistance_pa_s_per_m2: 0.0,
                distal_resistance_pa_s_per_m2: 1.0,
                compliance_m2_per_pa: if reference_density < 1.0 {
                    f64::from_bits(1)
                } else {
                    1.0
                },
                inertance_pa_s2_per_m2: if reference_density > 1.0 {
                    f64::from_bits(1)
                } else {
                    0.0
                },
            };
            assert!(matches!(
                field.bdf2_stokes_mixed_variable_step_windkessel_ports(
                    &zero_force,
                    &boundaries,
                    &[&outlet],
                    &[bad],
                    recent,
                    None,
                    dt,
                    1e-9,
                    1e-10,
                ),
                Err(MaterialFluxError::Invalid(_))
            ));
            assert_eq!(source.queries.get(), calls);
        }
    }
    assert_eq!((previous, older), snapshot);
}

/// A stationary two-layer Couette field carries flow while its external
/// compliance charges. Both material viscosities keep their own shear heat;
/// neither storage work nor either resistor is deposited into that bulk heat.
#[test]
fn layered_windkessel_preserves_viscosity_sources_and_external_storage() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let lower = Counted {
        material: transport(2.0, 0.2),
        queries: Cell::new(0),
    };
    let upper = Counted {
        material: transport(7.0, 1.4),
        queries: Cell::new(0),
    };
    let reference = transport(11.0, 0.4);
    let binding =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let sources: Vec<&dyn FluidTransport> = mesh
        .centroids
        .iter()
        .map(|point| {
            if point[1] < 0.5 {
                &lower as &dyn FluidTransport
            } else {
                &upper as &dyn FluidTransport
            }
        })
        .collect();
    let field = binding.bind_stokes_materials(sources).unwrap();
    let calls = (lower.queries.get(), upper.queries.get());
    let (tau, offset, drop, dt, previous_x) = (0.1, 0.3, 1.8, 0.07, 1.4);
    let slope = [tau / 0.2, tau / 1.4];
    let velocity = |point: [f64; 2]| {
        [
            if point[1] < 0.5 {
                offset + slope[0] * point[1]
            } else {
                offset + 0.5 * slope[0] + slope[1] * (point[1] - 0.5)
            },
            0.0,
        ]
    };
    let traction = |_: [f64; 2], n: [f64; 2]| [tau * n[1] - 3.0 * n[0], tau * n[0] - 3.0 * n[1]];
    let inlet: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].normal[0] < -0.5)
        .collect();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| {
            if mesh.edges[edge].normal[1].abs() > 0.5 {
                MaterialStokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let q = -(offset + 0.375 * slope[0] + 0.125 * slope[1]);
    let (rp, rd, compliance) = (0.3, 2.1, 0.6);
    let x = (compliance * previous_x / dt + q) / (compliance / dt + 1.0 / rd);
    let load = MaterialStokesWindkesselPort {
        reservoir_pressure_pa: drop - rp * q - x,
        proximal_resistance_pa_s_per_m2: rp,
        distal_resistance_pa_s_per_m2: rd,
        compliance_m2_per_pa: compliance,
        inertance_pa_s2_per_m2: 0.4,
    };
    let moments = flow_control_history(&mesh, &velocity);
    let solved = field
        .bdf2_stokes_mixed_variable_step_windkessel_ports(
            &|_| [-drop, 0.0],
            &boundaries,
            &[&inlet],
            &[load],
            MaterialStokesWindkesselHistory {
                velocity: &moments,
                compliance_pressures_pa: &[previous_x],
            },
            None,
            dt,
            1e-9,
            1e-10,
        )
        .unwrap();
    let (error, divergence) = system.velocity_error(&solved.flow.discrete().x, &velocity);
    assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
    near(solved.pressure_offsets_pa[0], drop, 3e-8);
    near(solved.outward_volume_rates_m2_s[0], q, 3e-9);
    near(solved.balances[0].compliance_pressure_pa, x, 3e-9);
    near(
        solved.balances[0].stored_energy_change_rate_w_per_m,
        0.5 * compliance * (x * x - previous_x * previous_x) / dt,
        3e-8,
    );
    near(
        solved.flow.viscous_dissipation().total_power_w_per_m(),
        0.5 * tau * tau * (1.0 / 0.2 + 1.0 / 1.4),
        3e-9,
    );
    near(solved.balances[0].power_balance_residual_w_per_m, 0.0, 3e-8);
    assert_eq!((lower.queries.get(), upper.queries.get()), calls);
}


/// Both reservoir connections retain SI normalization even when the source
/// densities differ from the common numerical reference density.
#[test]
fn all_port_resistances_preserve_si_pressure_inertia_and_source_ownership() {
    use fs_flux::MaterialStokesResistancePort;

    let (length, height, dt, previous_q) = (2.0, 1.5, 0.1, 0.2);
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(length, height, 2, 2));
    let system = FluxSystem::new(&mesh);
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -3.0 * n);
    let mut ports = [Vec::new(), Vec::new()];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            let normal = mesh.edges[edge].normal;
            if normal[1].abs() > 0.5 {
                MaterialStokesBoundary::NavierSlip {
                    edge,
                    value: &no_slip,
                    tangential_resistance_pa_s_per_m: 0.0,
                }
            } else {
                ports[usize::from(normal[0] > 0.0)].push(edge);
                MaterialStokesBoundary::TotalTraction { edge, value: &traction }
            }
        })
        .collect();
    let previous = flow_control_history(&mesh, &|_| [previous_q / height, 0.0]);
    for density in [2.0, 7.0] {
        let source = Counted {
            material: transport(density, 0.4),
            queries: Cell::new(0),
        };
        let uniform =
            MaterialFluxSystem::new(&system, &source, MaterialFluxNumerics::default()).unwrap();
        for reference_density in [0.25, 11.0] {
            let reference = transport(reference_density, 0.9);
            let binding =
                MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default())
                    .unwrap();
            let field = binding.bind_stokes_materials(vec![&source; mesh.tris.len()]).unwrap();
            let calls = source.queries.get();
            for drive in [-30.0, 30.0] {
                let loads = [
                    MaterialStokesResistancePort {
                        reservoir_pressure_pa: drive,
                        resistance_pa_s_per_m2: 3.0,
                    },
                    MaterialStokesResistancePort {
                        reservoir_pressure_pa: 0.0,
                        resistance_pa_s_per_m2: 7.0,
                    },
                ];
                let homogeneous = uniform.bdf2_stokes_mixed_variable_step_resistance_ports(
                    &zero_force, &boundaries, &[&ports[0], &ports[1]], &loads,
                    &previous, None, dt, 1e-9,
                ).unwrap();
                let layered = field.bdf2_stokes_mixed_variable_step_resistance_ports(
                    &zero_force, &boundaries, &[&ports[0], &ports[1]], &loads,
                    &previous, None, dt, 1e-9,
                ).unwrap();
                let fluid_inertance = density * length / height;
                let q = (fluid_inertance * previous_q / dt + drive)
                    / (fluid_inertance / dt + 10.0);
                for (flows, pressures, energy, heat, losses) in [
                    (
                        &homogeneous.outward_volume_rates_m2_s,
                        &homogeneous.pressure_offsets_pa,
                        homogeneous.flow.kinetic_energy_j_per_m(),
                        homogeneous.flow.viscous_dissipation().total_power_w_per_m(),
                        &homogeneous.resistance_dissipation_w_per_m,
                    ),
                    (
                        &layered.outward_volume_rates_m2_s,
                        &layered.pressure_offsets_pa,
                        layered.flow.kinetic_energy_j_per_m(),
                        layered.flow.viscous_dissipation().total_power_w_per_m(),
                        &layered.resistance_dissipation_w_per_m,
                    ),
                ] {
                    near(flows[0], -q, 3e-9);
                    near(flows[1], q, 3e-9);
                    near(flows.iter().sum(), 0.0, 3e-9);
                    near(pressures[0], drive - 3.0 * q, 3e-8);
                    near(pressures[1], 7.0 * q, 3e-8);
                    near(energy, 0.5 * fluid_inertance * q * q, 3e-9);
                    near(heat, 0.0, 1e-15);
                    near(losses.iter().sum(), 10.0 * q * q, 3e-9);
                    let fluid_change =
                        0.5 * fluid_inertance * (q * q - previous_q * previous_q) / dt;
                    let bdf_loss = 0.5 * fluid_inertance * (q - previous_q).powi(2) / dt;
                    near(drive * q, fluid_change + bdf_loss + losses.iter().sum::<f64>(), 3e-8);
                }
                for cell in 0..mesh.tris.len() {
                    assert!(std::ptr::eq(&field.samples()[cell], &layered.flow.samples()[cell]));
                }
            }
            assert_eq!(source.queries.get(), calls);
        }
    }
}
