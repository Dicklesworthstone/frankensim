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
