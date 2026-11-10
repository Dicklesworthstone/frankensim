//! Exact SI constitutive owners at a pressure-driven nonlinear endpoint.

use super::*;

#[test]
fn bdf2_material_rheology_preserves_endpoint_sources_pressure_and_instantaneous_heat() {
    // The legacy equal-step entry and explicit equal/growing/shrinking
    // intervals must all retain the same physical endpoint and source owners.
    for previous_timestep_s in [None, Some(0.1), Some(0.05), Some(0.2)] {
        check_bdf2_material_sources(previous_timestep_s);
    }
}

fn check_bdf2_material_sources(previous_timestep_s: Option<f64>) {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let dt = 0.1;
    let time = 3.0 * dt;
    let amplitude = |t: f64| 1.0 + 0.2 * t + 0.3 * t * t;
    let scale = amplitude(time);
    let base = |p: [f64; 2]| [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1];
    // Exact normal moments of the affine base field. No initialization solve
    // or constitutive query is used to manufacture either velocity history.
    let base_moments: Vec<_> = mesh
        .edges
        .iter()
        .flat_map(|edge| {
            let u = base(edge.mid);
            let a = mesh.verts[edge.verts.0];
            let b = mesh.verts[edge.verts.1];
            let delta = [b[0] - a[0], b[1] - a[1]];
            [
                u[0] * edge.normal[0] + u[1] * edge.normal[1],
                ((0.2 * delta[0] + 0.3 * delta[1]) * edge.normal[0]
                    - 0.2 * delta[1] * edge.normal[1])
                    / 12.0,
            ]
        })
        .collect();
    let previous: Vec<_> = base_moments
        .iter()
        .map(|u| amplitude(time - dt) * u)
        .collect();
    let older: Vec<_> = base_moments
        .iter()
        .map(|u| amplitude(time - dt - previous_timestep_s.unwrap_or(dt)) * u)
        .collect();
    let snapshots = [previous.clone(), older.clone()];
    let target = |p| base(p).map(|u| scale * u);
    let viscosity = |rate: f64| 1.2 / (1.0 + 0.8 * rate);
    let endpoint_rate = 0.5 * scale;
    for density in [2.0, 7.0] {
        let reference = transport(density, viscosity(0.0));
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-11,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let force = |p: [f64; 2]| {
            let v = base(p);
            [
                density * ((0.2 + 0.6 * time) * v[0] + scale * scale * (0.04 * p[0] + 0.11)) + 2.0,
                density * ((0.2 + 0.6 * time) * v[1] + scale * scale * (0.04 * p[1] - 0.02)) - 3.0,
            ]
        };
        let traction = |p: [f64; 2], n: [f64; 2]| {
            let pressure = 7.0 + 2.0 * p[0] - 3.0 * p[1];
            let mu = viscosity(endpoint_rate);
            [
                (0.4 * scale * mu - pressure) * n[0] + 0.3 * scale * mu * n[1],
                0.3 * scale * mu * n[0] + (-0.4 * scale * mu - pressure) * n[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
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
        let reads = Rc::new(Cell::new(0));
        let mut first_rates = vec![None; mesh.tris.len()];
        let mut queries = 0;
        let controls = NavierStokesControls {
            relative_tolerance: 1e-10,
            ..NavierStokesControls::default()
        };
        let resolve = |point: MaterialFluxRheologyPoint, retained: &FluidTransportState| {
            assert_eq!(retained, &reference);
            first_rates[point.cell].get_or_insert(point.shear_rate_per_s);
            queries += 1;
            Ok(RheologyResponse {
                transport: transport(density, viscosity(point.shear_rate_per_s)),
                point,
                temperature_k: 310.0,
                reads: Rc::clone(&reads),
            })
        };
        let result = if let Some(previous_dt) = previous_timestep_s {
            owner.bdf2_navier_stokes_mixed_rheology_variable_step(
                &force,
                &boundaries,
                &previous,
                Some((&older, previous_dt)),
                dt,
                controls,
                resolve,
            )
        } else {
            owner.bdf2_navier_stokes_mixed_rheology_step(
                &force,
                &boundaries,
                &previous,
                Some(&older),
                dt,
                controls,
                resolve,
            )
        }
        .unwrap();
        let solution = &result.step.solution;
        let (error, divergence) = system.velocity_error(&solution.discrete().x, &target);
        assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
        for (pressure, p) in solution.pressure_pa().iter().zip(&mesh.centroids) {
            near(*pressure, 7.0 + 2.0 * p[0] - 3.0 * p[1], 4e-8);
        }
        for initial_rate in first_rates {
            near(initial_rate.unwrap(), 0.5 * amplitude(time - dt), 1e-12);
        }
        let rates = system.shear_rates(&solution.discrete().x);
        for (cell, sample) in result.samples().iter().enumerate() {
            assert_eq!(sample.point(), sample.material().point);
            assert_eq!(sample.point().shear_rate_per_s, rates[cell]);
            near(
                sample.transport().dynamic_viscosity_pa_s(),
                viscosity(endpoint_rate),
                3e-8,
            );
            near(sample.material().temperature_k, 310.0, 0.0);
        }
        // q=2 mu D:D=mu gamma² is a rate at the accepted endpoint. It is
        // neither a BDF2 history difference nor a step-integrated heat source.
        let heat = viscosity(endpoint_rate) * endpoint_rate * endpoint_rate;
        for actual in solution.viscous_dissipation().cell_heating_w_m3() {
            near(*actual, heat, 4e-8);
        }
        near(
            solution.viscous_dissipation().total_power_w_per_m(),
            heat,
            4e-8,
        );
        assert_eq!(
            reads.get(),
            queries,
            "reporting must retain original source responses"
        );
        assert_eq!(solution.transport(), &reference);
    }
    assert_eq!(previous, snapshots[0]);
    assert_eq!(older, snapshots[1]);
}

#[test]
fn mixed_material_rheology_retains_endpoint_sources_and_coherent_pressure_inertia() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let dt = 0.025;
    let target = |p: [f64; 2]| [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1];
    let previous = vec![0.0; system.n_u];
    // gamma = sqrt(4*.2² + .3²) = .5. This synthetic source is shear thinning,
    // while temperature is retained source context, not an evolved field.
    let viscosity = |rate: f64| 1.2 / (1.0 + 0.8 * rate);
    for (density, policy, steady) in [
        (2.0, OpenBoundaryPolicy::NaturalTraction, false),
        (7.0, OpenBoundaryPolicy::DissipativeBackflow, false),
        (2.0, OpenBoundaryPolicy::NaturalTraction, true),
        (7.0, OpenBoundaryPolicy::DissipativeBackflow, true),
    ] {
        let reference = transport(density, viscosity(0.0));
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-11,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let inverse_dt = if steady { 0.0 } else { 1.0 / dt };
        let force = |p: [f64; 2]| {
            let u = target(p);
            [
                density * (u[0] * inverse_dt + 0.04 * p[0] + 0.11) + 2.0,
                density * (u[1] * inverse_dt + 0.04 * p[1] - 0.02) - 3.0,
            ]
        };
        let traction = |p: [f64; 2], n: [f64; 2]| {
            let pressure = 7.0 + 2.0 * p[0] - 3.0 * p[1];
            let mu = viscosity(0.5);
            let u = target(p);
            let resistance = match policy {
                OpenBoundaryPolicy::NaturalTraction => 0.0,
                OpenBoundaryPolicy::DissipativeBackflow => {
                    density * (u[0] * n[0] + u[1] * n[1]).min(0.0)
                }
            };
            [
                (0.4 * mu - pressure) * n[0] + 0.3 * mu * n[1] - resistance * u[0],
                0.3 * mu * n[0] + (-0.4 * mu - pressure) * n[1] - resistance * u[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[0].abs() > 0.5 {
                    MaterialNavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: policy,
                    }
                } else {
                    MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &target,
                    }
                }
            })
            .collect();
        let reads = Rc::new(Cell::new(0));
        let mut queries = 0;
        let mut last_query = vec![None; mesh.tris.len()];
        let controls = NavierStokesControls {
            relative_tolerance: 1e-10,
            relaxation: 0.8,
            ..NavierStokesControls::default()
        };
        let resolve = |point: MaterialFluxRheologyPoint, retained: &FluidTransportState| {
            assert_eq!(retained, &reference);
            queries += 1;
            last_query[point.cell] = Some(point);
            Ok(RheologyResponse {
                transport: transport(density, viscosity(point.shear_rate_per_s)),
                point,
                temperature_k: 310.0,
                reads: Rc::clone(&reads),
            })
        };
        let result = if steady {
            owner.solve_navier_stokes_mixed_rheology(
                &force,
                &boundaries,
                &previous,
                controls,
                resolve,
            )
        } else {
            owner.bdf1_navier_stokes_mixed_rheology_step(
                &force,
                &boundaries,
                &previous,
                dt,
                controls,
                resolve,
            )
        }
        .unwrap();
        let step = &result.step;
        let (error, divergence) = system.velocity_error(&step.solution.discrete().x, &target);
        assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
        assert!(step.iterations > 1 && step.nonlinear_residual <= 1e-10);
        assert_eq!(
            step.solution.pressure_reference(),
            StokesPressureReference::TractionDatum
        );
        assert!(
            step.open_boundary_policies
                .iter()
                .all(|(_, actual)| *actual == policy)
        );
        for (pressure, point) in step.solution.pressure_pa().iter().zip(&mesh.centroids) {
            near(*pressure, 7.0 + 2.0 * point[0] - 3.0 * point[1], 4e-8);
        }
        assert_eq!(step.solution.transport(), &reference);
        // Include the normal strains as well as shear: gamma²=.25.
        // The zero-shear reference viscosity would give the wrong heat.
        let heating = step.solution.viscous_dissipation();
        for &cell_heat in heating.cell_heating_w_m3() {
            near(cell_heat, viscosity(0.5) * 0.25, 5e-8);
        }
        near(heating.total_power_w_per_m(), viscosity(0.5) * 0.25, 5e-8);
        assert_eq!(result.samples().len(), mesh.tris.len());
        assert_eq!(
            reads.get(),
            queries,
            "each retained source is read only once"
        );
        assert!(queries > mesh.tris.len());
        assert!(
            queries <= 2 * mesh.tris.len() * step.iterations,
            "quadrature must reuse cell sources"
        );
        let rates = system.shear_rates(&step.solution.discrete().x);
        for (cell, sample) in result.samples().iter().enumerate() {
            let point = sample.point();
            assert_eq!(point.cell, cell);
            assert_eq!(point.position_m, mesh.centroids[cell]);
            assert_eq!(Some(point), last_query[cell]);
            assert_eq!(point, sample.material().point);
            assert_eq!(point.shear_rate_per_s, rates[cell]);
            assert_eq!(sample.transport(), &sample.material().transport);
            near(point.shear_rate_per_s, 0.5, 3e-8);
            near(
                sample.transport().dynamic_viscosity_pa_s(),
                viscosity(0.5),
                3e-8,
            );
            near(sample.material().temperature_k, 310.0, 0.0);
        }
        // The accepted source viscosity also owns dimensional pressure/shear
        // reactions. This affine field has nonzero normal viscous stress.
        let origin = [0.2, 0.3];
        for (edge_id, edge) in mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge_id, edge)| mesh.is_boundary(*edge_id) && edge.normal[1].abs() > 0.5)
        {
            let load = step
                .solution
                .boundary_loads(&[edge_id], origin, &[None])
                .unwrap();
            let mu = viscosity(0.5);
            let pressure = 7.0 + 2.0 * edge.mid[0] - 3.0 * edge.mid[1];
            let expected = [
                edge.len * ((pressure - 0.4 * mu) * edge.normal[0] - 0.3 * mu * edge.normal[1]),
                edge.len * (-0.3 * mu * edge.normal[0] + (pressure + 0.4 * mu) * edge.normal[1]),
            ];
            let a = mesh.verts[edge.verts.0];
            let b = mesh.verts[edge.verts.1];
            let delta = [b[0] - a[0], b[1] - a[1]];
            let moment = (edge.mid[0] - origin[0]) * expected[1]
                - (edge.mid[1] - origin[1]) * expected[0]
                + edge.len
                    * (2.0 * delta[0] - 3.0 * delta[1])
                    * (delta[0] * edge.normal[1] - delta[1] * edge.normal[0])
                    / 12.0;
            near(load.stress.force_n_per_m[0], expected[0], 1e-7);
            near(load.stress.force_n_per_m[1], expected[1], 1e-7);
            near(load.stress.moment_n, moment, 1e-7);
            near(
                load.stress.rigid_power_w_per_m([0.3, -0.2], 0.4).unwrap(),
                0.3 * expected[0] - 0.2 * expected[1] + 0.4 * moment,
                1e-7,
            );
            assert_eq!(load.pressure_offsets_pa, vec![(0, None)]);
        }
        // Analytic integral of |u|²/2 over the unit square, independent of the
        // FEEC mass matrix. Every dimensional term uses this same density.
        let x_square = 0.2_f64.powi(2) / 3.0
            + 0.3_f64.powi(2) / 3.0
            + 0.4_f64.powi(2)
            + 0.2 * 0.3 / 2.0
            + 0.2 * 0.4
            + 0.3 * 0.4;
        let y_square = 0.2_f64.powi(2) / 3.0 - 0.2 * 0.1 + 0.1_f64.powi(2);
        near(
            step.solution.kinetic_energy_j_per_m(),
            0.5 * density * (x_square + y_square),
            3e-8,
        );
        assert_eq!(
            reads.get(),
            queries,
            "reporting must not requery the source"
        );
        // Identical physical loads with the initial source viscosity cannot
        // reproduce the pressure-driven endpoint of the constitutive solve.
        let frozen = if steady {
            owner.solve_navier_stokes_mixed(
                &force,
                &boundaries,
                &previous,
                NavierStokesControls::default(),
            )
        } else {
            owner.bdf1_navier_stokes_mixed_step(
                &force,
                &boundaries,
                &previous,
                dt,
                NavierStokesControls::default(),
            )
        }
        .unwrap();
        assert!(
            system
                .velocity_error(&frozen.solution.discrete().x, &target)
                .0
                > 1e-4
        );
    }
    assert!(previous.iter().all(|value| *value == 0.0));
}

#[test]
fn mixed_material_rheology_preserves_endpoint_source_errors_and_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 1.0);
    let owner =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let previous = vec![0.0; system.n_u];
    let older = previous.clone();
    let target = |point: [f64; 2]| [0.3 * point[1], 0.0];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| MaterialNavierStokesBoundary::Velocity {
            edge,
            value: &target,
        })
        .collect();
    let controls = NavierStokesControls {
        maximum_iterations: 1,
        ..NavierStokesControls::default()
    };
    let refusal = MaterialError::State {
        what: "endpoint strain outside original source support".into(),
    };
    let mut queries = 0;
    // A budget of one forces the independent endpoint assembly. Only that
    // nonzero-strain response refuses, after the initial frozen solve succeeds.
    let result = owner.bdf2_navier_stokes_mixed_rheology_step(
        &zero_force,
        &boundaries,
        &previous,
        Some(&older),
        0.1,
        controls,
        |point, _| {
            queries += 1;
            if point.shear_rate_per_s > 1e-6 {
                Err(refusal.clone())
            } else {
                Ok(reference)
            }
        },
    );
    assert!(matches!(result, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert!(queries > mesh.tris.len() && queries <= 2 * mesh.tris.len());
    queries = 0;
    let result = owner.bdf2_navier_stokes_mixed_rheology_variable_step(
        &zero_force,
        &boundaries,
        &previous,
        Some((&older, 0.2)),
        0.1,
        controls,
        |point, _| {
            queries += 1;
            if point.shear_rate_per_s > 1e-6 {
                Err(refusal.clone())
            } else {
                Ok(reference)
            }
        },
    );
    assert!(matches!(result, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert!(queries > mesh.tris.len() && queries <= 2 * mesh.tris.len());
    for previous_dt in [0.0, 0.04] {
        let result = owner.bdf2_navier_stokes_mixed_rheology_variable_step(
            &zero_force,
            &boundaries,
            &previous,
            Some((&older, previous_dt)),
            0.1,
            controls,
            |_, _| -> Result<FluidTransportState, MaterialError> {
                panic!("inadmissible variable timestep queried the material source")
            },
        );
        assert!(matches!(
            result,
            Err(MaterialFluxError::Boundary(
                StokesBoundaryError::InvalidInput(_)
            ))
        ));
    }
    let invalid_older = vec![f64::NAN; system.n_u];
    let result = owner.bdf2_navier_stokes_mixed_rheology_step(
        &zero_force,
        &boundaries,
        &previous,
        Some(&invalid_older),
        0.1,
        controls,
        |_, _| -> Result<FluidTransportState, MaterialError> {
            panic!("invalid BDF2 history queried material source")
        },
    );
    assert!(matches!(
        result,
        Err(MaterialFluxError::Boundary(
            StokesBoundaryError::InvalidInput(_)
        ))
    ));
    assert_eq!(
        previous, older,
        "source refusal must preserve both histories"
    );
    let result = owner.bdf1_navier_stokes_mixed_rheology_step(
        &zero_force,
        &boundaries,
        &previous,
        0.1,
        controls,
        |point, _| {
            Ok(if point.shear_rate_per_s > 1e-6 {
                transport(3.0, 1.0)
            } else {
                reference
            })
        },
    );
    assert!(matches!(result, Err(MaterialFluxError::Invalid(_))));
    let result = owner.bdf1_navier_stokes_mixed_rheology_step(
        &zero_force,
        &boundaries,
        &previous,
        0.1,
        controls,
        |_, _| Ok(reference),
    );
    assert!(matches!(
        result,
        Err(MaterialFluxError::Boundary(StokesBoundaryError::Rheology(
            FluxRheologyError::NonlinearConvergence { iterations: 1, .. }
        )))
    ));
    queries = 0;
    let result = owner.bdf1_navier_stokes_mixed_rheology_step(
        &zero_force,
        &boundaries,
        &previous,
        f64::NAN,
        controls,
        |_, _| {
            queries += 1;
            Ok(reference)
        },
    );
    assert!(matches!(
        result,
        Err(MaterialFluxError::Boundary(
            StokesBoundaryError::InvalidInput(_)
        ))
    ));
    assert_eq!(
        queries, 0,
        "invalid history/time must refuse before any source query"
    );
    let result = owner.solve_navier_stokes_mixed_rheology(
        &zero_force,
        &boundaries,
        &previous,
        controls,
        |point, _| {
            if point.shear_rate_per_s > 1e-6 {
                Err(refusal.clone())
            } else {
                Ok(reference)
            }
        },
    );
    assert!(matches!(result, Err(MaterialFluxError::Material(error)) if error == refusal));
    let result = owner.solve_navier_stokes_mixed_rheology(
        &zero_force,
        &boundaries,
        &previous,
        controls,
        |point, _| {
            Ok(if point.shear_rate_per_s > 1e-6 {
                transport(3.0, 1.0)
            } else {
                reference
            })
        },
    );
    assert!(matches!(result, Err(MaterialFluxError::Invalid(_))));
    let invalid_guess = vec![f64::NAN; system.n_u];
    let result = owner.solve_navier_stokes_mixed_rheology(
        &zero_force,
        &boundaries,
        &invalid_guess,
        controls,
        |_, _| -> Result<FluidTransportState, MaterialError> {
            panic!("invalid steady guess queried material source")
        },
    );
    assert!(matches!(
        result,
        Err(MaterialFluxError::Boundary(
            StokesBoundaryError::InvalidInput(_)
        ))
    ));
    assert!(previous.iter().all(|value| *value == 0.0));
}

fn rheology_control_moments(mesh: &TriMesh, velocity: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    // Exact original BDM normal moments for an affine velocity; history is
    // authored independently of the nonlinear solve and constitutive source.
    mesh.edges
        .iter()
        .flat_map(|edge| {
            let a = velocity(mesh.verts[edge.verts.0]);
            let b = velocity(mesh.verts[edge.verts.1]);
            let na = a[0] * edge.normal[0] + a[1] * edge.normal[1];
            let nb = b[0] * edge.normal[0] + b[1] * edge.normal[1];
            [0.5 * (na + nb), (nb - na) / 12.0]
        })
        .collect()
}

#[test]
fn prescribed_transient_flow_closes_history_shear_convection_pressure_and_reactions() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let base = |p: [f64; 2]| [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1];
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, geometry)| mesh.is_boundary(*edge) && geometry.normal[0] < -0.5)
        .map(|(edge, _)| edge)
        .collect();
    let time = 0.3;
    let dt = 0.1;
    // BE differentiates linear amplitude exactly; BDF2 differentiates quadratic
    // amplitude exactly at each actual old time, including growth/shrinkage.
    for previous_dt in [None, Some(0.05), Some(0.1), Some(0.2)] {
        let amplitude = |t: f64| {
            1.0 + 0.2 * t
                + if previous_dt.is_some() {
                    0.3 * t * t
                } else {
                    0.0
                }
        };
        let scale = amplitude(time);
        let derivative = 0.2
            + if previous_dt.is_some() {
                0.6 * time
            } else {
                0.0
            };
        let previous =
            rheology_control_moments(&mesh, &|p| base(p).map(|u| amplitude(time - dt) * u));
        let older = previous_dt.map(|k| {
            rheology_control_moments(&mesh, &|p| base(p).map(|u| amplitude(time - dt - k) * u))
        });
        let snapshots = (previous.clone(), older.clone());
        let target = |p| base(p).map(|u| scale * u);
        let endpoint_rate = 0.5 * scale;
        let viscosity = |rate: f64| 1.2 / (1.0 + 0.8 * rate);
        for (density, policy) in [
            (2.0, OpenBoundaryPolicy::NaturalTraction),
            (7.0, OpenBoundaryPolicy::DissipativeBackflow),
        ] {
            let reference = transport(density, viscosity(0.0));
            let owner = MaterialFluxSystem::new(
                &system,
                &reference,
                MaterialFluxNumerics {
                    linear_tolerance: 1e-11,
                    ..MaterialFluxNumerics::default()
                },
            )
            .unwrap();
            // For u=A(t)*base, convection is A²*[.04x+.11,.04y-.02].
            // The body source omits rho*A' in x; the required pressure drop
            // must supply that exact inertial contribution. Viscosity is
            // spatially uniform at this affine endpoint, so div(2mu D)=0.
            let force = |p: [f64; 2]| {
                let v = base(p);
                [
                    density * (derivative * (v[0] - 1.0) + scale * scale * (0.04 * p[0] + 0.11)),
                    density * (derivative * v[1] + scale * scale * (0.04 * p[1] - 0.02)) - 3.0,
                ]
            };
            let traction = |p: [f64; 2], n: [f64; 2]| {
                let pressure = 7.0 - 3.0 * p[1];
                let mu = viscosity(endpoint_rate);
                let u = target(p);
                let backflow = match policy {
                    OpenBoundaryPolicy::NaturalTraction => 0.0,
                    OpenBoundaryPolicy::DissipativeBackflow => {
                        density * (u[0] * n[0] + u[1] * n[1]).min(0.0)
                    }
                };
                [
                    (0.4 * scale * mu - pressure) * n[0] + 0.3 * scale * mu * n[1]
                        - backflow * u[0],
                    0.3 * scale * mu * n[0] + (-0.4 * scale * mu - pressure) * n[1]
                        - backflow * u[1],
                ]
            };
            let boundaries: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, geometry)| {
                    if geometry.normal[0].abs() > 0.5 {
                        MaterialNavierStokesBoundary::Open {
                            edge,
                            traction: &traction,
                            backflow: policy,
                        }
                    } else {
                        MaterialNavierStokesBoundary::Velocity {
                            edge,
                            value: &target,
                        }
                    }
                })
                .collect();
            let reads = Rc::new(Cell::new(0));
            let mut first_rates = vec![None; mesh.tris.len()];
            let mut queries = 0;
            let resolve = |point: MaterialFluxRheologyPoint, retained: &FluidTransportState| {
                assert_eq!(retained, &reference);
                first_rates[point.cell].get_or_insert(point.shear_rate_per_s);
                queries += 1;
                Ok(RheologyResponse {
                    transport: transport(density, viscosity(point.shear_rate_per_s)),
                    point,
                    temperature_k: 310.0,
                    reads: Rc::clone(&reads),
                })
            };
            let controls = NavierStokesControls {
                maximum_iterations: 100,
                relative_tolerance: 1e-10,
                relaxation: 0.8,
            };
            let result = if previous_dt.is_some() {
                owner.bdf2_navier_stokes_mixed_rheology_variable_step_flow_control(
                    &force,
                    &boundaries,
                    &[&inlet],
                    &[-0.55 * scale],
                    &previous,
                    older.as_deref().zip(previous_dt),
                    dt,
                    controls,
                    1e-12,
                    resolve,
                )
            } else {
                owner.bdf1_navier_stokes_mixed_rheology_flow_control(
                    &force,
                    &boundaries,
                    &[&inlet],
                    &[-0.55 * scale],
                    &previous,
                    dt,
                    controls,
                    1e-12,
                    resolve,
                )
            }
            .unwrap();
            let solution = &result.flow.step.solution;
            let (error, divergence) = system.velocity_error(&solution.discrete().x, &target);
            assert!(error < 4e-9 && divergence < 4e-9, "{error}, {divergence}");
            near(result.pressure_offsets_pa[0], density * derivative, 8e-8);
            near(result.outward_volume_rates_m2_s[0], -0.55 * scale, 1e-11);
            assert!(result.flow.step.nonlinear_residual <= controls.relative_tolerance);
            assert!(result.flow.step.relative_update <= controls.relative_tolerance);
            assert!(
                result
                    .flow
                    .step
                    .open_boundary_policies
                    .iter()
                    .all(|(_, actual)| *actual == policy)
            );
            for (pressure, p) in solution.pressure_pa().iter().zip(&mesh.centroids) {
                near(
                    *pressure,
                    7.0 - 3.0 * p[1] + density * derivative * (1.0 - p[0]),
                    8e-8,
                );
            }
            for first in first_rates {
                near(first.unwrap(), 0.5 * amplitude(time - dt), 1e-12);
            }
            let rates = system.shear_rates(&solution.discrete().x);
            for (cell, sample) in result.flow.samples().iter().enumerate() {
                assert_eq!(sample.point(), sample.material().point);
                assert_eq!(
                    sample.point().shear_rate_per_s.to_bits(),
                    rates[cell].to_bits()
                );
                near(
                    sample.transport().dynamic_viscosity_pa_s(),
                    viscosity(endpoint_rate),
                    3e-8,
                );
                assert_eq!(sample.transport().density_kg_m3(), density);
                assert_eq!(sample.material().temperature_k, 310.0);
            }
            // Integral |base|²=131/300; physical heat is instantaneous mu*gamma².
            near(
                solution.kinetic_energy_j_per_m(),
                0.5 * density * scale * scale * 131.0 / 300.0,
                5e-8,
            );
            near(
                solution.viscous_dissipation().total_power_w_per_m(),
                viscosity(endpoint_rate) * endpoint_rate.powi(2),
                4e-8,
            );
            let mu = viscosity(endpoint_rate);
            let load = solution.boundary_loads(&inlet, [0.0; 2], &[None]).unwrap();
            near(
                load.stress.force_n_per_m[0],
                0.4 * mu * scale - 5.5 - density * derivative,
                1e-7,
            );
            near(load.stress.force_n_per_m[1], 0.3 * mu * scale, 1e-7);
            near(
                load.stress.moment_n,
                2.5 + 0.5 * density * derivative - 0.2 * mu * scale,
                1e-7,
            );
            assert_eq!(
                reads.get(),
                queries,
                "reporting retains the accepted source samples"
            );
            assert_eq!(previous, snapshots.0);
            assert_eq!(older, snapshots.1);
        }
    }
}

#[test]
fn pulsating_material_flow_pressure_retains_be_and_variable_bdf2_temporal_orders() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 2.3;
    let viscosity = |rate: f64| 1.2 / (1.0 + 0.8 * rate);
    let reference = transport(density, viscosity(0.0));
    let owner = MaterialFluxSystem::new(
        &system,
        &reference,
        MaterialFluxNumerics {
            linear_tolerance: 1e-12,
            ..MaterialFluxNumerics::default()
        },
    )
    .unwrap();
    let time: f64 = 0.7;
    let omega: f64 = 3.0;
    let offset = (omega * time).sin();
    let target = |p: [f64; 2]| [offset + 0.2 + 0.4 * p[1], 0.0];
    let shear = viscosity(0.4) * 0.4;
    let traction =
        |_: [f64; 2], n: [f64; 2]| [shear * n[1] - 3.2 * n[0], shear * n[0] - 3.2 * n[1]];
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, geometry)| mesh.is_boundary(*edge) && geometry.normal[0] < -0.5)
        .map(|(edge, _)| edge)
        .collect();
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
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
        relative_tolerance: 1e-10,
        ..NavierStokesControls::default()
    };
    let mut errors = [Vec::new(), Vec::new()];
    for dt in [0.02, 0.01, 0.005] {
        let previous_dt = 2.0 * dt;
        let recent = (omega * (time - dt)).sin();
        let old = (omega * (time - dt - previous_dt)).sin();
        let previous = rheology_control_moments(&mesh, &|p| [recent + 0.2 + 0.4 * p[1], 0.0]);
        let older = rheology_control_moments(&mesh, &|p| [old + 0.2 + 0.4 * p[1], 0.0]);
        let snapshots = [previous.clone(), older.clone()];
        for (method, errors) in errors.iter_mut().enumerate() {
            let resolve = |point: MaterialFluxRheologyPoint, _: &FluidTransportState| {
                Ok(transport(density, viscosity(point.shear_rate_per_s)))
            };
            let result = if method == 0 {
                owner.bdf1_navier_stokes_mixed_rheology_flow_control(
                    &zero_force,
                    &boundaries,
                    &[&inlet],
                    &[-(offset + 0.4)],
                    &previous,
                    dt,
                    controls,
                    1e-12,
                    resolve,
                )
            } else {
                owner.bdf2_navier_stokes_mixed_rheology_variable_step_flow_control(
                    &zero_force,
                    &boundaries,
                    &[&inlet],
                    &[-(offset + 0.4)],
                    &previous,
                    Some((&older, previous_dt)),
                    dt,
                    controls,
                    1e-12,
                    resolve,
                )
            }
            .unwrap();
            let derivative = if method == 0 {
                (offset - recent) / dt
            } else {
                (2.0 * dt + previous_dt) / (dt * (dt + previous_dt)) * offset
                    - (dt + previous_dt) / (dt * previous_dt) * recent
                    + dt / (previous_dt * (dt + previous_dt)) * old
            };
            near(result.pressure_offsets_pa[0], density * derivative, 8e-8);
            let solution = &result.flow.step.solution;
            let (error, divergence) = system.velocity_error(&solution.discrete().x, &target);
            assert!(error < 4e-9 && divergence < 4e-9, "{error}, {divergence}");
            near(
                solution.viscous_dissipation().total_power_w_per_m(),
                viscosity(0.4) * 0.16,
                4e-8,
            );
            errors.push(
                (result.pressure_offsets_pa[0] - density * omega * (omega * time).cos()).abs(),
            );
        }
        assert_eq!(previous, snapshots[0]);
        assert_eq!(older, snapshots[1]);
    }
    for (method, samples) in errors.iter().enumerate() {
        for pair in samples.windows(2) {
            assert!(
                pair[0] / pair[1] > if method == 0 { 1.9 } else { 3.6 },
                "{method}: {samples:?}"
            );
        }
    }
}

#[test]
fn transient_flow_control_preserves_sources_and_both_histories_on_refusal() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 1.0);
    let owner =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let target = |p: [f64; 2]| [0.2 + 0.4 * p[1], 0.0];
    let previous = rheology_control_moments(&mesh, &|p| [0.1 + 0.2 * p[1], 0.0]);
    let older = rheology_control_moments(&mesh, &|p| [0.05 + 0.1 * p[1], 0.0]);
    let snapshots = [previous.clone(), older.clone()];
    let calls = Cell::new(0);
    let force = |_: [f64; 2]| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    let traction = |_: [f64; 2], n: [f64; 2]| {
        calls.set(calls.get() + 1);
        [0.4 * n[1] - 3.2 * n[0], 0.4 * n[0] - 3.2 * n[1]]
    };
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, geometry)| mesh.is_boundary(*edge) && geometry.normal[0] < -0.5)
        .map(|(edge, _)| edge)
        .collect();
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
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
    let controls = NavierStokesControls::default();
    let mut queries = 0;
    let failed = owner.bdf2_navier_stokes_mixed_rheology_variable_step_flow_control(
        &force,
        &boundaries,
        &[&inlet],
        &[-0.4],
        &previous,
        Some((&older, 0.01)),
        0.1,
        controls,
        1e-12,
        |_, _| {
            queries += 1;
            Ok(reference)
        },
    );
    assert!(matches!(
        failed,
        Err(MaterialFluxError::Boundary(
            StokesBoundaryError::InvalidInput(_)
        ))
    ));
    assert_eq!(
        (queries, calls.get()),
        (0, 0),
        "history refusal precedes physical callbacks"
    );
    let refusal = crate::density::source_refusal();
    let failed = owner.bdf2_navier_stokes_mixed_rheology_variable_step_flow_control(
        &force,
        &boundaries,
        &[&inlet],
        &[-0.4],
        &previous,
        Some((&older, 0.1)),
        0.1,
        controls,
        1e-12,
        |point, _| {
            queries += 1;
            if point.shear_rate_per_s > 0.25 {
                Err(refusal.clone())
            } else {
                Ok(reference)
            }
        },
    );
    assert!(matches!(failed, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert!(
        queries > mesh.tris.len(),
        "the original source refuses the updated shear"
    );
    let failed = owner.bdf2_navier_stokes_mixed_rheology_variable_step_flow_control(
        &force,
        &boundaries,
        &[&inlet],
        &[-0.4],
        &previous,
        Some((&older, 0.1)),
        0.1,
        controls,
        1e-12,
        |_, _| Ok(transport(3.0, 1.0)),
    );
    assert!(matches!(failed, Err(MaterialFluxError::Invalid(_))));
    let failed = owner.bdf2_navier_stokes_mixed_rheology_variable_step_flow_control(
        &force,
        &boundaries,
        &[&inlet],
        &[-0.4],
        &previous,
        Some((&older, 0.1)),
        0.1,
        NavierStokesControls {
            maximum_iterations: 1,
            ..controls
        },
        1e-12,
        |_, _| Ok(reference),
    );
    assert!(matches!(
        failed,
        Err(MaterialFluxError::Boundary(StokesBoundaryError::Rheology(
            fs_flux::FluxRheologyError::NonlinearConvergence { iterations: 1, .. }
        )))
    ));
    assert_eq!(previous, snapshots[0]);
    assert_eq!(older, snapshots[1]);
}
