//! G1/G3: original cell material states drive inertia, stress and SI reporting.

use super::*;

struct Counted<G> {
    material: G,
    queries: Cell<usize>,
}

impl<G: FluidTransport> FluidTransport for Counted<G> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        self.material.transport_state()
    }
}

#[test]
fn layered_source_density_and_viscosity_share_inertia_stress_heat_and_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let spec = GasSpec::dry_air_ussa1976();
    let lower = Counted {
        material: GasState::try_new(&spec, 293.15, 101_325.0).unwrap(),
        queries: Cell::new(0),
    };
    let upper = Counted {
        material: GasState::try_new(&spec, 353.15, 303_975.0).unwrap(),
        queries: Cell::new(0),
    };
    let states = [
        lower.material.transport_state().unwrap(),
        upper.material.transport_state().unwrap(),
    ];
    let rho = states.map(|state| state.density_kg_m3());
    let mu = states.map(|state| state.dynamic_viscosity_pa_s());
    assert!(rho[1] > 2.0 * rho[0] && mu[1] != mu[0]);
    let shear_stress_rate = 0.8 * mu[0];
    let slope = mu.map(|viscosity| shear_stress_rate / viscosity);
    let offset = [0.7, 0.7 + 0.5 * slope[0]];
    let rate = |y: f64| {
        if y < 0.5 {
            offset[0] + slope[0] * y
        } else {
            offset[1] + slope[1] * (y - 0.5)
        }
    };
    let gradient = 1.3;
    let slip_resistance = 1.7 * mu[0];
    let all_edges: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .collect();
    // u=t*rate(y), p=5-gradient*x. Both rho and mu jump at y=1/2,
    // but mu*du/dy=t*shear_stress_rate is continuous across that interface.
    // The flow is parallel to the density layers, so fixed density is physical.
    let momentum_rate: f64 = (0..2)
        .map(|layer| rho[layer] * (0.5 * offset[layer] + 0.125 * slope[layer]))
        .sum();
    let energy_coefficient: f64 = (0..2)
        .map(|layer| {
            0.5 * rho[layer]
                * (0.5 * offset[layer].powi(2)
                    + 0.25 * offset[layer] * slope[layer]
                    + slope[layer].powi(2) / 24.0)
        })
        .sum();
    for reference_density in [0.6, 6.2] {
        // The normalization source is deliberately different from both layers.
        // Neither its density nor its viscosity may replace a cell's material.
        let reference = transport(reference_density, 0.7);
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
        let before = [lower.queries.get(), upper.queries.get()];
        let field = binding.bind_stokes_materials(sources).unwrap();
        assert_eq!(field.reference_density_kg_m3(), reference_density);
        assert_eq!(lower.queries.get() - before[0], 4);
        assert_eq!(upper.queries.get() - before[1], 4);
        let admitted_queries = [lower.queries.get(), upper.queries.get()];
        let mut previous = vec![0.0; system.n_u];
        let mut time = 0.0;
        for dt in [0.17, 0.11] {
            time += dt;
            let velocity = |point: [f64; 2]| [time * rate(point[1]), 0.0];
            let wall = |point: [f64; 2]| {
                let sign = if point[1] > 0.5 { 1.0 } else { -1.0 };
                [
                    time * (rate(point[1]) + sign * shear_stress_rate / slip_resistance),
                    0.0,
                ]
            };
            let force = |point: [f64; 2]| {
                let layer = usize::from(point[1] >= 0.5);
                [rho[layer] * rate(point[1]) - gradient, 0.0]
            };
            let traction = |point: [f64; 2], normal: [f64; 2]| {
                let pressure = 5.0 - gradient * point[0];
                let shear = time * shear_stress_rate;
                [
                    shear * normal[1] - pressure * normal[0],
                    shear * normal[0] - pressure * normal[1],
                ]
            };
            let boundaries: Vec<_> = all_edges
                .iter()
                .map(|&edge| {
                    if mesh.edges[edge].normal[1].abs() > 0.5 {
                        MaterialStokesBoundary::NavierSlip {
                            edge,
                            value: &wall,
                            tangential_resistance_pa_s_per_m: slip_resistance,
                        }
                    } else {
                        MaterialStokesBoundary::TotalTraction {
                            edge,
                            value: &traction,
                        }
                    }
                })
                .collect();
            let history = previous.clone();
            let solved = field
                .bdf1_stokes_mixed_step(&force, &boundaries, &previous, dt)
                .unwrap();
            assert_eq!(previous, history);
            let (error, divergence) = system.velocity_error(&solved.discrete().x, &velocity);
            assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
            assert_eq!(
                solved.pressure_reference(),
                StokesPressureReference::TractionDatum
            );
            let mut momentum = 0.0;
            let mut expected_heat = 0.0;
            for (cell, point) in mesh.centroids.iter().enumerate() {
                let layer = usize::from(point[1] >= 0.5);
                near(solved.pressure_pa()[cell], 5.0 - gradient * point[0], 3e-8);
                assert_eq!(
                    solved.pressure_reference_at(cell),
                    Some(StokesPressureReference::TractionDatum)
                );
                assert_eq!(*solved.samples()[cell].transport(), states[layer]);
                assert!(std::ptr::eq(
                    &field.samples()[cell],
                    &solved.samples()[cell]
                ));
                let source: &dyn FluidTransport = if layer == 0 { &lower } else { &upper };
                assert!(std::ptr::addr_eq(
                    *solved.samples()[cell].material(),
                    source
                ));
                momentum += rho[layer]
                    * mesh.areas[cell]
                    * system.velocity_at(&solved.discrete().x, *point)[0];
                let heating = mu[layer] * (time * slope[layer]).powi(2);
                near(
                    solved.viscous_dissipation().cell_heating_w_m3()[cell],
                    heating,
                    1e-11,
                );
                expected_heat += mesh.areas[cell] * heating;
            }
            near(momentum, time * momentum_rate, 3e-9);
            near(
                solved.kinetic_energy_j_per_m(),
                time * time * energy_coefficient,
                3e-9,
            );
            near(
                solved.viscous_dissipation().total_power_w_per_m(),
                expected_heat,
                1e-11,
            );
            let load = solved
                .boundary_loads(&all_edges, [0.0; 2], &[None])
                .unwrap();
            near(load.stress.force_n_per_m[0], -gradient, 3e-8);
            near(load.stress.force_n_per_m[1], 0.0, 3e-8);
            near(load.stress.moment_n, 0.5 * gradient, 3e-8);
            near(
                load.slip_dissipation_w_per_m,
                2.0 * (time * shear_stress_rate).powi(2) / slip_resistance,
                1e-11,
            );
            assert_eq!([lower.queries.get(), upper.queries.get()], admitted_queries);
            previous.copy_from_slice(solved.velocity_dofs());
        }
    }
}

struct SwitchSource {
    state: FluidTransportState,
    queries: Cell<usize>,
    refuse: Cell<bool>,
}

fn source_refusal() -> MaterialError {
    MaterialError::State {
        what: "original source outside admitted conditions".into(),
    }
}

impl FluidTransport for SwitchSource {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        if self.refuse.get() {
            Err(source_refusal())
        } else {
            Ok(self.state)
        }
    }
}

#[test]
fn material_field_retains_source_errors_and_frozen_snapshots() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 0.4);
    let binding =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let source = SwitchSource {
        state: transport(7.0, 0.9),
        queries: Cell::new(0),
        refuse: Cell::new(false),
    };
    assert!(matches!(
        binding.bind_stokes_materials(vec![&source]),
        Err(MaterialFluxError::Invalid(_))
    ));
    assert_eq!(
        source.queries.get(),
        0,
        "shape admission precedes source evaluation"
    );
    source.refuse.set(true);
    let refusal = binding.bind_stokes_materials(vec![&source, &source]);
    assert!(
        matches!(refusal, Err(MaterialFluxError::Material(error)) if error == source_refusal())
    );
    assert_eq!(
        source.queries.get(),
        1,
        "no query after the first original refusal"
    );
    source.refuse.set(false);
    let field = binding
        .bind_stokes_materials(vec![&source, &source])
        .unwrap();
    assert_eq!(source.queries.get(), 3);
    source.refuse.set(true);
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| MaterialStokesBoundary::Velocity {
            edge,
            value: &no_slip,
        })
        .collect();
    let previous = vec![0.0; system.n_u];
    let solved = field
        .bdf1_stokes_mixed_step(&zero_force, &boundaries, &previous, 0.1)
        .unwrap();
    assert_eq!(
        source.queries.get(),
        3,
        "prepared timesteps must not requery a material"
    );
    drop(field);
    assert!(
        solved
            .velocity_dofs()
            .iter()
            .all(|value| value.abs() < 1e-12)
    );
    for sample in solved.samples() {
        assert!(std::ptr::eq(*sample.material(), &source));
        assert_eq!(*sample.transport(), source.state);
    }
    near(solved.kinetic_energy_j_per_m(), 0.0, 1e-15);
}

#[test]
fn spatial_binding_refuses_unrepresentable_common_normalization() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let reference = FluidTransportState::try_new(1e-100, 1.0, 1e-100, 1e-100).unwrap();
    let source = FluidTransportState::try_new(1e300, 1.0, 1e300, 1e300).unwrap();
    let binding =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    assert!(matches!(
        binding.bind_stokes_materials(vec![source, source]),
        Err(MaterialFluxError::Invalid(_))
    ));
}

#[test]
fn sourced_bdf2_keeps_two_histories_and_physical_normalization() {
    for previous_dt in [0.25, 0.5] {
        sourced_bdf2_physical_normalization(previous_dt);
    }
}

fn sourced_bdf2_physical_normalization(previous_dt: f64) {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let spec = GasSpec::dry_air_ussa1976();
    let lower = Counted {
        material: GasState::try_new(&spec, 293.15, 101_325.0).unwrap(),
        queries: Cell::new(0),
    };
    let upper = Counted {
        material: GasState::try_new(&spec, 353.15, 303_975.0).unwrap(),
        queries: Cell::new(0),
    };
    let densities = [
        lower.material.transport_state().unwrap().density_kg_m3(),
        upper.material.transport_state().unwrap().density_kg_m3(),
    ];
    // u(t)=t² on both layers. Histories use their actual times; both constant
    // and doubled intervals must recover u(1)=1 from endpoint u'=2.
    let mut previous = vec![0.0; system.n_u];
    let mut older = vec![0.0; system.n_u];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        previous[2 * edge] = 0.25 * geometry.normal[0];
        older[2 * edge] = (0.5 - previous_dt).powi(2) * geometry.normal[0];
    }
    let snapshots = [previous.clone(), older.clone()];
    let traction = |point: [f64; 2], normal: [f64; 2]| {
        let pressure = 5.0 - 0.8 * point[0];
        normal.map(|n| -pressure * n)
    };
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| MaterialStokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let force = |point: [f64; 2]| [2.0 * densities[usize::from(point[1] >= 0.5)] - 0.8, 0.0];
    for reference_density in [2.0, 9.0] {
        let reference = transport(reference_density, 0.7);
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
        let queries = [lower.queries.get(), upper.queries.get()];
        let solved = field
            .bdf2_stokes_mixed_variable_step(
                &force,
                &boundaries,
                &previous,
                Some((&older, previous_dt)),
                0.5,
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&solved.discrete().x, &|_| [1.0, 0.0]);
        assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
        for (cell, point) in mesh.centroids.iter().enumerate() {
            near(solved.pressure_pa()[cell], 5.0 - 0.8 * point[0], 3e-9);
            assert!(std::ptr::eq(
                &field.samples()[cell],
                &solved.samples()[cell]
            ));
        }
        near(
            solved.kinetic_energy_j_per_m(),
            0.25 * (densities[0] + densities[1]),
            3e-10,
        );
        near(
            solved.viscous_dissipation().total_power_w_per_m(),
            0.0,
            1e-15,
        );
        assert_eq!([lower.queries.get(), upper.queries.get()], queries);

        // The homogeneous source adapter reaches the same native time formula.
        let homogeneous = binding
            .bdf2_stokes_mixed_variable_step(
                &|_| [2.0 * reference_density - 0.8, 0.0],
                &boundaries,
                &previous,
                Some((&older, previous_dt)),
                0.5,
            )
            .unwrap();
        assert!(
            system
                .velocity_error(&homogeneous.discrete().x, &|_| [1.0, 0.0])
                .0
                < 3e-10
        );
        near(
            homogeneous.kinetic_energy_j_per_m(),
            0.5 * reference_density,
            3e-10,
        );
    }
    assert_eq!(previous, snapshots[0]);
    assert_eq!(older, snapshots[1]);
}

#[test]
fn layered_flow_control_retains_source_viscosity_energy_pressure_and_loads() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 4));
    let system = FluxSystem::new(&mesh);
    let lower = Counted {
        material: transport(2.0, 1.0),
        queries: Cell::new(0),
    };
    let upper = Counted {
        material: transport(5.0, 4.0),
        queries: Cell::new(0),
    };
    // Continuous shear stress 0.8: u=0.8y below and u=0.3+0.2y above.
    // The exact integral is 0.325 m²/s; only the walls prescribe velocity.
    let velocity = |point: [f64; 2]| {
        [
            if point[1] < 0.5 {
                0.8 * point[1]
            } else {
                0.3 + 0.2 * point[1]
            },
            0.0,
        ]
    };
    let traction = |_: [f64; 2], n: [f64; 2]| [0.8 * n[1] - 3.2 * n[0], 0.8 * n[0] - 3.2 * n[1]];
    let inlet: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].normal[0] < -0.5)
        .collect();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
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
    for reference_density in [1.3, 8.0] {
        let reference = transport(reference_density, 0.7);
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
        let queries = [lower.queries.get(), upper.queries.get()];
        let controlled = field
            .solve_stokes_mixed_flow_control(
                &|_| [-1.7, 0.0],
                &boundaries,
                &[&inlet],
                &[-0.325],
                1e-12,
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&controlled.flow.discrete().x, &velocity);
        assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
        near(controlled.pressure_offsets_pa[0], 1.7, 3e-9);
        near(controlled.outward_volume_rates_m2_s[0], -0.325, 1e-12);
        for (cell, point) in mesh.centroids.iter().enumerate() {
            near(
                controlled.flow.pressure_pa()[cell],
                3.2 + 1.7 * (1.0 - point[0]),
                3e-9,
            );
            assert!(std::ptr::eq(
                &field.samples()[cell],
                &controlled.flow.samples()[cell]
            ));
        }
        // Independent integrals of rho*u²/2 and mu*(du/dy)² on the two layers.
        near(
            controlled.flow.kinetic_energy_j_per_m(),
            337.0 / 1200.0,
            3e-10,
        );
        near(
            controlled.flow.viscous_dissipation().total_power_w_per_m(),
            0.4,
            3e-10,
        );
        let load = controlled
            .flow
            .boundary_loads(&inlet, [0.0; 2], &[None])
            .unwrap();
        near(load.stress.force_n_per_m[0], -4.9, 3e-9);
        near(load.stress.force_n_per_m[1], 0.8, 3e-9);
        near(load.stress.moment_n, 2.45, 3e-9);
        assert_eq!([lower.queries.get(), upper.queries.get()], queries);
    }
}

#[test]
fn homogeneous_flow_control_uses_integrated_area_and_physical_pressure_units() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 2.0, 2, 4));
    let system = FluxSystem::new(&mesh);
    let velocity = |point: [f64; 2]| [0.25 * point[1], 0.0];
    let traction = |_: [f64; 2], n: [f64; 2]| [0.5 * n[1] - 3.2 * n[0], 0.5 * n[0] - 3.2 * n[1]];
    let inlet: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].normal[0] < -0.5)
        .collect();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
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
    for density in [0.6, 6.2] {
        let source = transport(density, 2.0);
        let binding =
            MaterialFluxSystem::new(&system, &source, MaterialFluxNumerics::default()).unwrap();
        let controlled = binding
            .solve_stokes_mixed_flow_control(
                &|_| [-1.7, 0.0],
                &boundaries,
                &[&inlet],
                &[-0.5],
                1e-12,
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&controlled.flow.discrete().x, &velocity);
        assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
        near(controlled.pressure_offsets_pa[0], 1.7, 3e-9);
        near(controlled.outward_volume_rates_m2_s[0], -0.5, 1e-12);
        near(
            controlled.flow.kinetic_energy_j_per_m(),
            density / 12.0,
            3e-10,
        );
        near(
            controlled.flow.viscous_dissipation().total_power_w_per_m(),
            0.25,
            3e-10,
        );
        let load = controlled
            .flow
            .boundary_loads(&inlet, [0.0; 2], &[None])
            .unwrap();
        near(load.stress.force_n_per_m[0], -9.8, 3e-9);
    }
}

// Exact BDM1 moments of a continuous piecewise-affine horizontal field.
// The interface is aligned to the mesh, so each edge trace is affine.
fn flow_control_history(mesh: &TriMesh, velocity: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mut moments = vec![0.0; 2 * mesh.edges.len()];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        let a = velocity(mesh.verts[geometry.verts.0]);
        let b = velocity(mesh.verts[geometry.verts.1]);
        let normal_a = a[0] * geometry.normal[0] + a[1] * geometry.normal[1];
        let normal_b = b[0] * geometry.normal[0] + b[1] * geometry.normal[1];
        moments[2 * edge] = 0.5 * (normal_a + normal_b);
        moments[2 * edge + 1] = (normal_b - normal_a) / 12.0;
    }
    moments
}

#[test]
fn flow_driven_startup_pressure_retains_previous_momentum() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let inlet: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].normal[0] < -0.5)
        .collect();
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -3.2 * n);
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            if mesh.edges[edge].normal[1].abs() > 0.5 {
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
    let resting = vec![0.0; system.n_u];
    let target_velocity = 0.7;
    for (density, dt) in [(0.6, 0.12), (6.2, 0.3)] {
        let source = transport(density, 0.8);
        let binding =
            MaterialFluxSystem::new(&system, &source, MaterialFluxNumerics::default()).unwrap();
        let started = binding
            .bdf1_stokes_mixed_flow_control(
                &zero_force,
                &boundaries,
                &[&inlet],
                &[-target_velocity],
                &resting,
                dt,
                1e-12,
            )
            .unwrap();
        let exact_drop = density * target_velocity / dt;
        let (error, divergence) =
            system.velocity_error(&started.flow.discrete().x, &|_| [target_velocity, 0.0]);
        assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
        near(started.pressure_offsets_pa[0], exact_drop, 3e-9);
        near(
            started.outward_volume_rates_m2_s[0],
            -target_velocity,
            1e-12,
        );
        near(
            started.flow.kinetic_energy_j_per_m(),
            0.5 * density * target_velocity.powi(2),
            3e-10,
        );
        for (cell, point) in mesh.centroids.iter().enumerate() {
            near(
                started.flow.pressure_pa()[cell],
                3.2 + exact_drop * (1.0 - point[0]),
                3e-9,
            );
        }
        let load = started
            .flow
            .boundary_loads(&inlet, [0.0; 2], &[None])
            .unwrap();
        near(load.stress.force_n_per_m[0], -(3.2 + exact_drop), 3e-9);

        // Holding the flow needs no second acceleration impulse: the accepted
        // velocity supplies its actual momentum to the next pressure solve.
        let held = binding
            .bdf1_stokes_mixed_flow_control(
                &zero_force,
                &boundaries,
                &[&inlet],
                &[-target_velocity],
                started.flow.velocity_dofs(),
                0.5 * dt,
                1e-12,
            )
            .unwrap();
        near(held.pressure_offsets_pa[0], 0.0, 3e-9);
        assert!(
            system
                .velocity_error(&held.flow.discrete().x, &|_| [target_velocity, 0.0])
                .0
                < 3e-10
        );
    }
    assert!(resting.iter().all(|&value| value == 0.0));
}

#[test]
fn pulsating_port_pressure_has_second_order_variable_step_inertia() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 2.3;
    let viscosity = 0.7;
    let source = transport(density, viscosity);
    let binding =
        MaterialFluxSystem::new(&system, &source, MaterialFluxNumerics::default()).unwrap();
    let time: f64 = 0.7;
    let omega: f64 = 3.0;
    let offset = (omega * time).sin();
    let velocity = |point: [f64; 2]| [offset + 0.2 + 0.4 * point[1], 0.0];
    let shear = viscosity * 0.4;
    let traction =
        |_: [f64; 2], n: [f64; 2]| [shear * n[1] - 3.2 * n[0], shear * n[0] - 3.2 * n[1]];
    let inlet: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].normal[0] < -0.5)
        .collect();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
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
    let mut pressure_errors = Vec::new();
    for dt in [0.02, 0.01, 0.005] {
        let previous_dt = 2.0 * dt;
        let previous_offset = (omega * (time - dt)).sin();
        let older_offset = (omega * (time - dt - previous_dt)).sin();
        let previous = flow_control_history(&mesh, &|point| {
            [previous_offset + 0.2 + 0.4 * point[1], 0.0]
        });
        let older =
            flow_control_history(&mesh, &|point| [older_offset + 0.2 + 0.4 * point[1], 0.0]);
        let snapshots = [previous.clone(), older.clone()];
        let controlled = binding
            .bdf2_stokes_mixed_variable_step_flow_control(
                &zero_force,
                &boundaries,
                &[&inlet],
                &[-(offset + 0.4)],
                &previous,
                Some((&older, previous_dt)),
                dt,
                1e-12,
            )
            .unwrap();
        // Independent derivative of the quadratic interpolant at the three
        // actual times; no equal-step approximation or startup is admissible.
        let derivative = (2.0 * dt + previous_dt) / (dt * (dt + previous_dt)) * offset
            - (dt + previous_dt) / (dt * previous_dt) * previous_offset
            + dt / (previous_dt * (dt + previous_dt)) * older_offset;
        near(
            controlled.pressure_offsets_pa[0],
            density * derivative,
            3e-8,
        );
        let (error, divergence) = system.velocity_error(&controlled.flow.discrete().x, &velocity);
        assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
        near(
            controlled.flow.viscous_dissipation().total_power_w_per_m(),
            viscosity * 0.4_f64.powi(2),
            3e-10,
        );
        let pressure_error =
            (controlled.pressure_offsets_pa[0] - density * omega * (omega * time).cos()).abs();
        pressure_errors.push(pressure_error);
        assert_eq!(previous, snapshots[0]);
        assert_eq!(older, snapshots[1]);
    }
    for pair in pressure_errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.6, "{pressure_errors:?}");
    }
}

#[test]
fn layered_transient_flow_control_keeps_weighted_history_energy_and_si_reactions() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 4));
    let system = FluxSystem::new(&mesh);
    let lower = Counted {
        material: transport(2.0, 1.0),
        queries: Cell::new(0),
    };
    let upper = Counted {
        material: transport(5.0, 4.0),
        queries: Cell::new(0),
    };
    let shape = |y: f64| if y < 0.5 { 0.8 * y } else { 0.3 + 0.2 * y };
    let traction = |_: [f64; 2], n: [f64; 2]| [0.8 * n[1] - 3.2 * n[0], 0.8 * n[0] - 3.2 * n[1]];
    let inlet: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].normal[0] < -0.5)
        .collect();
    for reference_density in [1.3, 8.0] {
        let reference = transport(reference_density, 0.7);
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
        let queries = [lower.queries.get(), upper.queries.get()];
        // BE exactly differentiates U=t; unequal-step BDF2 exactly
        // differentiates U=t² for step ratios 2, 1 and 1/2.
        for previous_dt in [None, Some(0.04), Some(0.08), Some(0.16)] {
            let (time, dt, offset, derivative) = if previous_dt.is_some() {
                (0.4_f64, 0.08_f64, 0.16_f64, 0.8_f64)
            } else {
                (0.2, 0.2, 0.2, 1.0)
            };
            let previous_offset = if previous_dt.is_some() {
                (time - dt).powi(2)
            } else {
                time - dt
            };
            let previous =
                flow_control_history(&mesh, &|point| [previous_offset + shape(point[1]), 0.0]);
            let older = previous_dt.map(|interval| {
                let older_offset = (time - dt - interval).powi(2);
                flow_control_history(&mesh, &|point| [older_offset + shape(point[1]), 0.0])
            });
            let snapshots = (previous.clone(), older.clone());
            let older_history = older.as_deref().zip(previous_dt);
            let velocity = |point: [f64; 2]| [offset + shape(point[1]), 0.0];
            let force = |point: [f64; 2]| {
                let rho = if point[1] < 0.5 { 2.0 } else { 5.0 };
                [(rho - 3.0) * derivative, 0.0]
            };
            let boundaries: Vec<_> = (0..mesh.edges.len())
                .filter(|&edge| mesh.is_boundary(edge))
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
            let controlled = if previous_dt.is_some() {
                field.bdf2_stokes_mixed_variable_step_flow_control(
                    &force,
                    &boundaries,
                    &[&inlet],
                    &[-(offset + 0.325)],
                    &previous,
                    older_history,
                    dt,
                    1e-12,
                )
            } else {
                field.bdf1_stokes_mixed_flow_control(
                    &force,
                    &boundaries,
                    &[&inlet],
                    &[-(offset + 0.325)],
                    &previous,
                    dt,
                    1e-12,
                )
            }
            .unwrap();
            let pressure_drop = 3.0 * derivative;
            near(controlled.pressure_offsets_pa[0], pressure_drop, 3e-9);
            let (error, divergence) =
                system.velocity_error(&controlled.flow.discrete().x, &velocity);
            assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
            for (cell, point) in mesh.centroids.iter().enumerate() {
                near(
                    controlled.flow.pressure_pa()[cell],
                    3.2 + pressure_drop * (1.0 - point[0]),
                    3e-9,
                );
                assert!(std::ptr::eq(
                    &field.samples()[cell],
                    &controlled.flow.samples()[cell]
                ));
            }
            near(
                controlled.flow.kinetic_energy_j_per_m(),
                337.0 / 1200.0 + 53.0 / 40.0 * offset + 7.0 / 4.0 * offset.powi(2),
                3e-10,
            );
            near(
                controlled.flow.viscous_dissipation().total_power_w_per_m(),
                0.4,
                3e-10,
            );
            let load = controlled
                .flow
                .boundary_loads(&inlet, [0.0; 2], &[None])
                .unwrap();
            near(load.stress.force_n_per_m[0], -(3.2 + pressure_drop), 3e-9);
            near(load.stress.force_n_per_m[1], 0.8, 3e-9);
            near(load.stress.moment_n, 0.5 * (3.2 + pressure_drop), 3e-9);

            // A rejected endpoint cannot replace either accepted history or
            // re-query/mutate a source. Exercise both temporal admission and
            // a flow-target refusal after the private response preparation.
            if let Some(older_velocity) = older.as_deref() {
                assert!(
                    field
                        .bdf2_stokes_mixed_variable_step_flow_control(
                            &force,
                            &boundaries,
                            &[&inlet],
                            &[-(offset + 0.325)],
                            &previous,
                            Some((older_velocity, dt / 3.0)),
                            dt,
                            1e-12,
                        )
                        .is_err()
                );
            }
            assert!(
                field
                    .bdf2_stokes_mixed_variable_step_flow_control(
                        &force,
                        &boundaries,
                        &[&inlet],
                        &[f64::NAN],
                        &previous,
                        older_history,
                        dt,
                        1e-12,
                    )
                    .is_err()
            );
            assert_eq!(previous, snapshots.0);
            assert_eq!(older, snapshots.1);
            assert_eq!([lower.queries.get(), upper.queries.get()], queries);
        }
    }
}
