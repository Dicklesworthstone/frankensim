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
