//! Coupled inertia/strain closure with an independently derived exact flow.

use super::*;

#[test]
fn pressure_flow_closes_convection_and_spatial_shear_viscosity_at_the_same_endpoint() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let dt = 0.02;
    let old = |p: [f64; 2]| [0.1 * p[0] + 0.1 * p[1] + 0.2, -0.1 * p[1] + 0.05];
    let target = |p: [f64; 2]| [0.3 * p[0] + 0.2 * p[1] + 0.4, -0.3 * p[1] + 0.1];
    let previous = moments(&mesh, &old);
    let history = previous.clone();
    let initial_rates = system.shear_rates(&previous);
    // Synthetic shear-thinning law: gamma² = 4*.3² + .2² = .4 for the
    // target. nu varies linearly in x, so its stress divergence is nonzero.
    let law = |_: usize, p: [f64; 2], rate: f64| (0.5 + 0.1 * p[0]) / (1.0 + 0.25 * rate);
    let target_rate = 0.4_f64.sqrt();
    let viscosity_slope = 0.1 / (1.0 + 0.25 * target_rate);
    let force = |p: [f64; 2]| {
        let u = target(p);
        let v = old(p);
        [
            (u[0] - v[0]) / dt + 0.09 * p[0] + 0.14 - 0.6 * viscosity_slope + 2.0,
            (u[1] - v[1]) / dt + 0.09 * p[1] - 0.03 - 0.2 * viscosity_slope - 3.0,
        ]
    };
    for (policy, shift) in [
        (OpenBoundaryPolicy::NaturalTraction, 0.0),
        (OpenBoundaryPolicy::DissipativeBackflow, 17.0),
    ] {
        let traction = |p: [f64; 2], n: [f64; 2]| {
            let pressure = 5.0 + shift + 2.0 * p[0] - 3.0 * p[1];
            let nu = law(0, p, target_rate);
            let u = target(p);
            // The supplied load excludes the declared added resistance. For
            // the dissipative policy, total stress = t + min(u.n,0)*u.
            let resistance = match policy {
                OpenBoundaryPolicy::NaturalTraction => 0.0,
                OpenBoundaryPolicy::DissipativeBackflow => (u[0] * n[0] + u[1] * n[1]).min(0.0),
            };
            [
                (0.6 * nu - pressure) * n[0] + 0.2 * nu * n[1] - resistance * u[0],
                0.2 * nu * n[0] + (-0.6 * nu - pressure) * n[1] - resistance * u[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[0].abs() > 0.5 {
                    NavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: policy,
                    }
                } else {
                    NavierStokesBoundary::Velocity {
                        edge,
                        value: &target,
                    }
                }
            })
            .collect();
        let step = system
            .bdf1_navier_stokes_mixed_rheology_step(
                params,
                &law,
                &force,
                &boundaries,
                &previous,
                dt,
                NavierStokesControls {
                    relative_tolerance: 1e-10,
                    relaxation: 0.8,
                    ..NavierStokesControls::default()
                },
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&step.solution.field.x, &target);
        assert!(error < 2e-9 && divergence < 2e-9, "{error}, {divergence}");
        assert!(step.iterations > 1 && step.nonlinear_residual <= 1e-10);
        assert_eq!(
            step.solution.pressure_reference,
            StokesPressureReference::TractionDatum
        );
        assert!(
            step.open_boundary_policies
                .iter()
                .all(|(_, actual)| *actual == policy)
        );
        for (pressure, point) in step.solution.pressure.iter().zip(&mesh.centroids) {
            assert!((pressure - (5.0 + shift + 2.0 * point[0] - 3.0 * point[1])).abs() < 3e-8);
        }
        let rates = system.shear_rates(&step.solution.field.x);
        assert!(rates.iter().all(|rate| (rate - target_rate).abs() < 2e-8));
        // Independent frozen-operator assembly at the returned velocity must
        // pass. Freezing EITHER strain or convection at history must fail.
        for (rate_field, advection, accepted) in [
            (rates.as_slice(), &step.solution.field.x[..system.n_u], true),
            (
                initial_rates.as_slice(),
                &step.solution.field.x[..system.n_u],
                false,
            ),
            (rates.as_slice(), previous.as_slice(), false),
        ] {
            let assembled = system
                .assemble_navier_stokes_mixed_bdf1_variable_viscosity(
                    params,
                    &|cell, point| law(cell, point, rate_field[cell]),
                    &force,
                    &boundaries,
                    advection,
                    &previous,
                    dt,
                )
                .unwrap();
            let mut residual = vec![0.0; system.n];
            assembled
                .matrix()
                .spmv(&step.solution.field.x, &mut residual);
            let norm = residual
                .iter()
                .zip(assembled.rhs())
                .fold(0.0_f64, |r, (ax, rhs)| r.hypot(ax - rhs));
            if accepted {
                assert!(norm < 2e-8, "endpoint defect {norm}");
            } else {
                assert!(
                    norm > 1e-5,
                    "lagging a live nonlinear term gave defect {norm}"
                );
            }
        }
    }
    assert_eq!(previous, history);
}

#[test]
fn mixed_rheology_rejects_invalid_coefficients_and_exhausted_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let previous = vec![0.0; system.n_u];
    let traction = |p: [f64; 2], n: [f64; 2]| {
        let pressure = 3.0 - p[0];
        [-pressure * n[0], -pressure * n[1]]
    };
    let boundaries = open(&mesh, &traction, OpenBoundaryPolicy::NaturalTraction);
    for viscosity in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let result = system.bdf1_navier_stokes_mixed_rheology_step(
            FluxParams::default(),
            &|_, _, _| viscosity,
            &|_| [0.0; 2],
            &boundaries,
            &previous,
            0.02,
            NavierStokesControls::default(),
        );
        assert!(matches!(
            result,
            Err(StokesBoundaryError::InvalidViscosity { .. })
        ));
    }
    let result = system.bdf1_navier_stokes_mixed_rheology_step(
        FluxParams::default(),
        &|_, _, rate| 1.0 / (1.0 + rate),
        &|_| [0.0; 2],
        &boundaries,
        &previous,
        0.02,
        NavierStokesControls {
            maximum_iterations: 1,
            ..NavierStokesControls::default()
        },
    );
    assert!(matches!(
        result,
        Err(StokesBoundaryError::Rheology(
            fs_flux::FluxRheologyError::NonlinearConvergence { iterations: 1, .. }
        ))
    ));
    assert!(previous.iter().all(|value| *value == 0.0));
}

#[test]
fn bdf2_quadratic_time_flow_closes_spatial_rheology_and_actual_history_strain() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let dt = 0.08;
    let time = 3.0 * dt;
    let coefficients = |t: f64| {
        [
            0.2 + 0.3 * t + 0.1 * t * t,
            0.35 + 0.2 * t - 0.15 * t * t,
            0.4 + 0.1 * t + 0.4 * t * t,
            0.1 - 0.05 * t + 0.2 * t * t,
        ]
    };
    let velocity = |t: f64, p: [f64; 2]| {
        let [a, b, c, d] = coefficients(t);
        [a * p[0] + b * p[1] + c, -a * p[1] + d]
    };
    let previous = moments(&mesh, &|p| velocity(time - dt, p));
    let older = moments(&mesh, &|p| velocity(time - 2.0 * dt, p));
    let snapshots = [previous.clone(), older.clone()];
    let previous_rates = system.shear_rates(&previous);
    let first_rates = std::cell::RefCell::new(vec![None; mesh.tris.len()]);
    let law = |cell: usize, p: [f64; 2], rate: f64| {
        first_rates.borrow_mut()[cell].get_or_insert(rate);
        (0.7 + 0.12 * p[0] + 0.08 * p[1]) / (1.0 + 0.25 * rate * rate)
    };
    let [a, b, c, d] = coefficients(time);
    let rate = (4.0 * a * a + b * b).sqrt();
    let denominator = 1.0 + 0.25 * rate * rate;
    let target = |p| velocity(time, p);
    // Independent differential equation: affine spatial velocity, quadratic
    // time, linear spatial viscosity. BDF2 differentiates this time law exactly.
    let force = |p: [f64; 2]| {
        let acceleration = [
            (0.3 + 0.2 * time) * p[0] + (0.2 - 0.3 * time) * p[1] + 0.1 + 0.8 * time,
            -(0.3 + 0.2 * time) * p[1] - 0.05 + 0.4 * time,
        ];
        [
            acceleration[0] + a * a * p[0] + a * c + b * d
                - (2.0 * a * 0.12 + b * 0.08) / denominator
                + 2.0,
            acceleration[1] + a * a * p[1]
                - a * d
                - (b * 0.12 - 2.0 * a * 0.08) / denominator
                - 3.0,
        ]
    };
    for policy in [
        OpenBoundaryPolicy::NaturalTraction,
        OpenBoundaryPolicy::DissipativeBackflow,
    ] {
        let traction = |p: [f64; 2], n: [f64; 2]| {
            let pressure = 5.0 + 2.0 * p[0] - 3.0 * p[1];
            let nu = (0.7 + 0.12 * p[0] + 0.08 * p[1]) / denominator;
            let u = target(p);
            let resistance = match policy {
                OpenBoundaryPolicy::NaturalTraction => 0.0,
                OpenBoundaryPolicy::DissipativeBackflow => (u[0] * n[0] + u[1] * n[1]).min(0.0),
            };
            [
                (2.0 * a * nu - pressure) * n[0] + b * nu * n[1] - resistance * u[0],
                b * nu * n[0] + (-2.0 * a * nu - pressure) * n[1] - resistance * u[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[0].abs() > 0.5 {
                    NavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: policy,
                    }
                } else {
                    NavierStokesBoundary::Velocity {
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
        let check = |step: &fs_flux::NavierStokesStep| {
            let (error, divergence) = system.velocity_error(&step.solution.field.x, &target);
            assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
            assert!(step.nonlinear_residual <= 1e-10);
            for (pressure, p) in step.solution.pressure.iter().zip(&mesh.centroids) {
                assert!((pressure - (5.0 + 2.0 * p[0] - 3.0 * p[1])).abs() < 4e-8);
            }
        };
        let mut equal_signature = None;
        // Exercise the legacy entry, exact step doubling, equal-step parity,
        // and a step reduction against the same differentiated PDE.
        for (previous_dt, legacy) in [
            (dt, true),
            (0.5 * dt, false),
            (dt, false),
            (2.0 * dt, false),
        ] {
            let history = if previous_dt == dt {
                older.clone()
            } else {
                moments(&mesh, &|p| velocity(time - dt - previous_dt, p))
            };
            let history_snapshot = history.clone();
            first_rates.borrow_mut().fill(None);
            let step = if legacy {
                system.bdf2_navier_stokes_mixed_rheology_step(
                    FluxParams::default(),
                    &law,
                    &force,
                    &boundaries,
                    &previous,
                    Some(&history),
                    dt,
                    controls,
                )
            } else {
                system.bdf2_navier_stokes_mixed_rheology_variable_step(
                    FluxParams::default(),
                    &law,
                    &force,
                    &boundaries,
                    &previous,
                    Some((&history, previous_dt)),
                    dt,
                    controls,
                )
            }
            .unwrap();
            check(&step);
            for (initial, expected) in first_rates.borrow().iter().zip(&previous_rates) {
                assert_eq!(
                    *initial,
                    Some(*expected),
                    "mass history is not an initial strain state"
                );
            }
            let signature = (
                step.solution
                    .field
                    .x
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                step.iterations,
                step.relative_update.to_bits(),
                step.nonlinear_residual.to_bits(),
            );
            if legacy {
                equal_signature = Some(signature);
            } else if previous_dt == dt {
                assert_eq!(&signature, equal_signature.as_ref().unwrap());
            }
            if previous_dt == 0.5 * dt && policy == OpenBoundaryPolicy::NaturalTraction {
                let frozen =
                    |_: usize, p: [f64; 2]| (0.7 + 0.12 * p[0] + 0.08 * p[1]) / denominator;
                let step = system
                    .bdf2_navier_stokes_mixed_variable_viscosity_variable_step(
                        FluxParams::default(),
                        &frozen,
                        &force,
                        &boundaries,
                        &previous,
                        Some((&history, previous_dt)),
                        dt,
                        controls,
                    )
                    .unwrap();
                check(&step);
            }
            assert_eq!(history, history_snapshot);
        }
    }
    assert_eq!(previous, snapshots[0]);
    assert_eq!(older, snapshots[1]);
}

#[test]
fn bdf2_rheology_has_second_order_temporal_refinement_with_bdf1_startup() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let base = |p: [f64; 2]| [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1];
    let law = |_: usize, p: [f64; 2], rate: f64| (0.6 + 0.1 * p[0]) / (1.0 + 0.25 * rate * rate);
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let dt = 0.4 / f64::from(steps);
        let mut previous = moments(&mesh, &base);
        let mut older = None;
        for step in 1..=steps {
            let time = f64::from(step) * dt;
            let amplitude = (-time).exp();
            let rate = 0.5 * amplitude;
            let viscosity_slope = 0.1 / (1.0 + 0.25 * rate * rate);
            let force = |p: [f64; 2]| {
                let u = base(p).map(|value| amplitude * value);
                [
                    -u[0] + amplitude * amplitude * (0.04 * p[0] + 0.11)
                        - 0.4 * amplitude * viscosity_slope
                        + 2.0,
                    -u[1] + amplitude * amplitude * (0.04 * p[1] - 0.02)
                        - 0.3 * amplitude * viscosity_slope
                        - 3.0,
                ]
            };
            let traction = |p: [f64; 2], n: [f64; 2]| {
                let pressure = 5.0 + 2.0 * p[0] - 3.0 * p[1];
                let nu = law(0, p, rate);
                [
                    (0.4 * amplitude * nu - pressure) * n[0] + 0.3 * amplitude * nu * n[1],
                    0.3 * amplitude * nu * n[0] + (-0.4 * amplitude * nu - pressure) * n[1],
                ]
            };
            let boundaries = open(&mesh, &traction, OpenBoundaryPolicy::NaturalTraction);
            let endpoint = system
                .bdf2_navier_stokes_mixed_rheology_step(
                    FluxParams::default(),
                    &law,
                    &force,
                    &boundaries,
                    &previous,
                    older.as_deref(),
                    dt,
                    NavierStokesControls {
                        relative_tolerance: 1e-11,
                        ..NavierStokesControls::default()
                    },
                )
                .unwrap();
            older = Some(previous);
            previous = endpoint.solution.field.x[..system.n_u].to_vec();
        }
        let error = system
            .velocity_error(&previous, &|p| base(p).map(|u| (-0.4_f64).exp() * u))
            .0;
        errors.push(error);
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.2 && pair[0] / pair[1] < 4.8,
            "{errors:?}"
        );
    }
}

#[test]
fn bdf2_rheology_startup_is_bdf1_and_bad_histories_never_query_a_source() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let previous = vec![0.0; system.n_u];
    let traction = |p: [f64; 2], n: [f64; 2]| {
        let pressure = 3.0 - p[0];
        [-pressure * n[0], -pressure * n[1]]
    };
    let boundaries = open(&mesh, &traction, OpenBoundaryPolicy::NaturalTraction);
    let law = |_: usize, _: [f64; 2], rate: f64| 1.0 / (1.0 + 0.2 * rate);
    let force = |_: [f64; 2]| [0.0; 2];
    let bdf1 = system
        .bdf1_navier_stokes_mixed_rheology_step(
            FluxParams::default(),
            &law,
            &force,
            &boundaries,
            &previous,
            0.02,
            NavierStokesControls::default(),
        )
        .unwrap();
    let startup = system
        .bdf2_navier_stokes_mixed_rheology_step(
            FluxParams::default(),
            &law,
            &force,
            &boundaries,
            &previous,
            None,
            0.02,
            NavierStokesControls::default(),
        )
        .unwrap();
    assert_eq!(startup.solution.field.x, bdf1.solution.field.x);
    assert_eq!(startup.iterations, bdf1.iterations);
    assert_eq!(startup.nonlinear_residual, bdf1.nonlinear_residual);
    let variable_startup = system
        .bdf2_navier_stokes_mixed_rheology_variable_step(
            FluxParams::default(),
            &law,
            &force,
            &boundaries,
            &previous,
            None,
            0.02,
            NavierStokesControls::default(),
        )
        .unwrap();
    assert_eq!(
        variable_startup
            .solution
            .field
            .x
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>(),
        bdf1.solution
            .field
            .x
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(variable_startup.iterations, bdf1.iterations);
    assert_eq!(
        variable_startup.nonlinear_residual.to_bits(),
        bdf1.nonlinear_residual.to_bits()
    );
    for (recent, old, dt) in [
        (previous.clone(), vec![0.0; system.n_u - 1], 0.02),
        (previous.clone(), vec![f64::NAN; system.n_u], 0.02),
        (
            vec![f64::MAX; system.n_u],
            vec![-f64::MAX; system.n_u],
            0.02,
        ),
        (previous.clone(), previous.clone(), 1.25 / f64::MAX),
    ] {
        let result = system.bdf2_navier_stokes_mixed_rheology_step(
            FluxParams::default(),
            &|_, _, _| panic!("invalid history queried viscosity"),
            &|_| panic!("invalid history queried force"),
            &boundaries,
            &recent,
            Some(&old),
            dt,
            NavierStokesControls::default(),
        );
        assert!(matches!(result, Err(StokesBoundaryError::InvalidInput(_))));
    }
    let no_traction =
        |_: [f64; 2], _: [f64; 2]| -> [f64; 2] { panic!("invalid interval queried traction") };
    let unqueried = open(&mesh, &no_traction, OpenBoundaryPolicy::NaturalTraction);
    for previous_dt in [
        0.0,
        -0.02,
        f64::NAN,
        f64::INFINITY,
        f64::from_bits(1),
        0.009,
        1e308,
    ] {
        let result = system.bdf2_navier_stokes_mixed_rheology_variable_step(
            FluxParams::default(),
            &|_, _, _| panic!("invalid interval queried viscosity"),
            &|_| panic!("invalid interval queried force"),
            &unqueried,
            &previous,
            Some((&previous, previous_dt)),
            0.02,
            NavierStokesControls::default(),
        );
        assert!(matches!(result, Err(StokesBoundaryError::InvalidInput(_))));
    }
    assert!(previous.iter().all(|value| *value == 0.0));
}
