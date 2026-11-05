//! G1/G3: weighted BDF2 time accuracy, exact history and two-state energy.

use super::{density_test_moments, traction_boundaries, velocity_boundaries};
use fs_flux::{FluxParams, FluxSystem, StokesBoundary, StokesBoundaryError, TriMesh};
use fs_solid::Mesh2;
use std::cell::Cell;

fn density(y: f64) -> f64 {
    1.0 + 2.0 * y + 3.0 * y * y
}

#[test]
fn quadratic_time_graded_shear_uses_both_histories_and_one_reference_density() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let dt = 0.13;
    let time = 3.0 * dt;
    let intercept = |t: f64| 0.3 + 1.1 * t + 0.4 * t * t;
    let slope = |t: f64| 0.4 + 0.6 * t - 0.2 * t * t;
    let previous = density_test_moments(&system, &|point| {
        [intercept(time - dt) + slope(time - dt) * point[1], 0.0]
    });
    let older = density_test_moments(&system, &|point| {
        [
            intercept(time - 2.0 * dt) + slope(time - 2.0 * dt) * point[1],
            0.0,
        ]
    });
    let snapshots = [previous.clone(), older.clone()];
    let velocity = |point: [f64; 2]| [intercept(time) + slope(time) * point[1], 0.0];
    let acceleration = |y: f64| 1.1 + 0.8 * time + (0.6 - 0.4 * time) * y;
    let edges: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .collect();
    // BDF2 differentiates quadratic time dependence exactly. rho is quadratic,
    // mu affine, velocity affine, so all spatial products fit native quadrature.
    let energy = 0.5
        * (3.0 * intercept(time).powi(2)
            + (23.0 / 6.0) * intercept(time) * slope(time)
            + (43.0 / 30.0) * slope(time).powi(2));
    for reference in [1.0, 5.0] {
        let samples = Cell::new(0);
        let rho = |_: usize, point: [f64; 2]| {
            samples.set(samples.get() + 1);
            density(point[1]) / reference
        };
        let force = |point: [f64; 2]| {
            [
                (density(point[1]) * acceleration(point[1]) - 1.2 * slope(time) - 2.0) / reference,
                0.0,
            ]
        };
        let traction = |point: [f64; 2], normal: [f64; 2]| {
            let pressure = 8.0 - 2.0 * point[0];
            let shear = (0.7 + 1.2 * point[1]) * slope(time);
            [
                (shear * normal[1] - pressure * normal[0]) / reference,
                (shear * normal[0] - pressure * normal[1]) / reference,
            ]
        };
        let boundaries: Vec<_> = edges
            .iter()
            .map(|&edge| {
                if mesh.edges[edge].normal[1].abs() > 0.5 {
                    StokesBoundary::Velocity {
                        edge,
                        value: &velocity,
                    }
                } else {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                }
            })
            .collect();
        let assembly = system
            .assemble_stokes_mixed_bdf2_variable_density(
                params,
                &rho,
                &|_, point| (0.7 + 1.2 * point[1]) / reference,
                &force,
                &boundaries,
                &previous,
                Some(&older),
                dt,
            )
            .unwrap();
        let solved = assembly.solve().unwrap();
        let (error, divergence) = system.velocity_error(&solved.field.x, &velocity);
        assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
        for (actual, point) in solved.pressure.iter().zip(&mesh.centroids) {
            assert!((reference * actual - (8.0 - 2.0 * point[0])).abs() < 3e-9);
        }
        assert!(
            (reference * assembly.kinetic_energy(&solved.field.x).unwrap() - energy).abs() < 3e-10
        );
        let load = solved
            .reactions
            .resultant(&edges, [0.0; 2], &[None])
            .unwrap();
        let stress_gradient = 2.0 + 1.2 * slope(time);
        assert!((reference * load.stress.force[0] + stress_gradient).abs() < 3e-9);
        assert!((reference * load.stress.force[1]).abs() < 3e-9);
        assert!((reference * load.stress.moment - 0.5 * stress_gradient).abs() < 3e-9);
        assert_eq!(
            samples.get(),
            6 * mesh.tris.len(),
            "both histories share one density sample set"
        );
        assert_eq!(previous, snapshots[0]);
        assert_eq!(older, snapshots[1]);
    }
}

#[test]
fn weighted_startup_is_original_bdf1_and_equal_step_g_energy_is_dissipative() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let rho = |_: usize, point: [f64; 2]| density(point[1]);
    let viscosity = |_: usize, _: [f64; 2]| 0.7;
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries = traction_boundaries(&mesh, &traction);
    let previous = density_test_moments(&system, &|_| [0.3, 0.0]);
    let older = vec![0.0; system.n_u];
    let bdf1 = system
        .assemble_stokes_mixed_bdf1_variable_density(
            params,
            &rho,
            &viscosity,
            &zero,
            &boundaries,
            &previous,
            0.2,
        )
        .unwrap();
    let startup = system
        .assemble_stokes_mixed_bdf2_variable_density(
            params,
            &rho,
            &viscosity,
            &zero,
            &boundaries,
            &previous,
            None,
            0.2,
        )
        .unwrap();
    assert_eq!(startup.matrix().to_dense(), bdf1.matrix().to_dense());
    assert_eq!(startup.rhs(), bdf1.rhs());
    let assembly = system
        .assemble_stokes_mixed_bdf2_variable_density(
            params,
            &rho,
            &viscosity,
            &zero,
            &boundaries,
            &previous,
            Some(&older),
            0.2,
        )
        .unwrap();
    let solved = assembly.solve().unwrap();
    let recent = &solved.field.x[..system.n_u];
    let energy = |v: &[f64]| assembly.kinetic_energy(v).unwrap();
    let extrapolate =
        |a: &[f64], b: &[f64]| -> Vec<f64> { a.iter().zip(b).map(|(a, b)| 2.0 * a - b).collect() };
    let second_difference: Vec<_> = recent
        .iter()
        .zip(&previous)
        .zip(&older)
        .map(|((a, b), c)| a - 2.0 * b + c)
        .collect();
    let old_g = 0.5 * (energy(&previous) + energy(&extrapolate(&previous, &older)));
    let new_g = 0.5 * (energy(recent) + energy(&extrapolate(recent, &previous)));
    let temporal_loss = 0.5 * energy(&second_difference);
    // With zero forcing and a uniform velocity, viscous work vanishes.
    // G_new-G_old + ||u_new-2u_n+u_(n-1)||_rho²/4 = 0 independently
    // of dt. Physical kinetic energy can rise while this G energy decreases.
    assert!((new_g - old_g + temporal_loss).abs() < 3e-12);
    assert!(new_g < old_g && energy(recent) > energy(&previous));
    assert!((energy(recent) - 1.5 * 0.4_f64.powi(2)).abs() < 3e-12);
    let freestream = system
        .assemble_stokes_mixed_bdf2_variable_density(
            params,
            &rho,
            &viscosity,
            &zero,
            &boundaries,
            &previous,
            Some(&previous),
            0.2,
        )
        .unwrap()
        .solve()
        .unwrap();
    let (error, divergence) = system.velocity_error(&freestream.field.x, &|_| [0.3, 0.0]);
    assert!(error < 3e-11 && divergence < 3e-11);
}

#[test]
fn original_bdf1_startup_gives_second_order_time_convergence() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let rho = |_: usize, point: [f64; 2]| density(point[1]);
    let traction = |point: [f64; 2], normal: [f64; 2]| normal.map(|n| -(4.0 - 2.0 * point[0]) * n);
    let boundaries = traction_boundaries(&mesh, &traction);
    let final_time = 0.4;
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let dt = final_time / steps as f64;
        let mut previous = vec![0.0; system.n_u];
        let mut older: Option<Vec<f64>> = None;
        for step in 1..=steps {
            let time = step as f64 * dt;
            // u'=t, exact u=t²/2; stress gradient supplies the uniform -2 term.
            let force = |point: [f64; 2]| [density(point[1]) * time - 2.0, 0.0];
            let solved = system
                .assemble_stokes_mixed_bdf2_variable_density(
                    params,
                    &rho,
                    &|_, _| 0.7,
                    &force,
                    &boundaries,
                    &previous,
                    older.as_deref(),
                    dt,
                )
                .unwrap()
                .solve()
                .unwrap();
            older = Some(previous);
            previous = solved.field.x[..system.n_u].to_vec();
        }
        let error = system.velocity_at(&previous, [0.4, 0.4])[0] - 0.5 * final_time * final_time;
        // The independently solved scalar error recurrence is
        // e_n=(3/4)*dt²*(1-3^(-n)), including the actual startup error.
        let expected_error = 0.75 * dt * dt * (1.0 - 3.0_f64.powi(-steps));
        assert!((error - expected_error).abs() < 3e-11);
        errors.push(error.abs());
    }
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.1);
    }
}

#[test]
fn weighted_g_energy_includes_original_viscous_and_sip_dissipation() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.7,
        sigma: 30.0,
        ..FluxParams::default()
    };
    let rho = |_: usize, point: [f64; 2]| density(point[1]);
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries = velocity_boundaries(&mesh, &zero);
    // curl of psi=x(1-x)y(1-y): zero boundary-normal velocity. Original
    // two-point edge Gauss integrates its cubic normal means exactly, so the
    // BDM projection remains exactly divergence-free despite tangential jumps.
    let previous = density_test_moments(&system, &|point| {
        let [x, y] = point;
        [
            x * (1.0 - x) * (1.0 - 2.0 * y),
            -(1.0 - 2.0 * x) * y * (1.0 - y),
        ]
    });
    let dt = 0.05;
    let assembly = system
        .assemble_stokes_mixed_bdf2_variable_density(
            params,
            &rho,
            &|_, _| params.nu,
            &zero,
            &boundaries,
            &previous,
            Some(&previous),
            dt,
        )
        .unwrap();
    let solved = assembly.solve().unwrap();
    let recent = &solved.field.x[..system.n_u];
    let extrapolated: Vec<_> = recent
        .iter()
        .zip(&previous)
        .map(|(a, b)| 2.0 * a - b)
        .collect();
    let second_difference: Vec<_> = recent.iter().zip(&previous).map(|(a, b)| a - b).collect();
    let energy = |v: &[f64]| assembly.kinetic_energy(v).unwrap();
    let old_g = energy(&previous);
    let new_g = 0.5 * (energy(recent) + energy(&extrapolated));
    let temporal_loss = 0.5 * energy(&second_difference);
    let spatial = system
        .assemble_stokes_mixed(params, &zero, &boundaries)
        .unwrap();
    let mut applied = vec![0.0; system.n];
    spatial.matrix().spmv(&solved.field.x, &mut applied);
    let viscous_power: f64 = solved.field.x.iter().zip(applied).map(|(x, y)| x * y).sum();
    assert!(viscous_power > 0.0 && new_g < old_g);
    assert!((new_g - old_g + temporal_loss + dt * viscous_power).abs() < 3e-11);
}

#[test]
fn bad_second_history_refuses_before_any_field_or_boundary_query() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let queries = Cell::new(0);
    let coefficient = |_: usize, _: [f64; 2]| {
        queries.set(queries.get() + 1);
        1.0
    };
    let value = |_: [f64; 2]| {
        queries.set(queries.get() + 1);
        [0.0; 2]
    };
    let boundaries = velocity_boundaries(&mesh, &value);
    for (recent, older, dt) in [
        (vec![0.0; system.n_u], Vec::new(), 0.1),
        (vec![0.0; system.n_u], vec![f64::NAN; system.n_u], 0.1),
        (vec![0.0; system.n_u], vec![0.0; system.n_u], 0.0),
        (vec![f64::MAX; system.n_u], vec![-f64::MAX; system.n_u], 0.1),
    ] {
        assert!(matches!(
            system.assemble_stokes_mixed_bdf2_variable_density(
                FluxParams::default(),
                &coefficient,
                &coefficient,
                &value,
                &boundaries,
                &recent,
                Some(&older),
                dt,
            ),
            Err(StokesBoundaryError::InvalidInput(_))
        ));
    }
    assert_eq!(queries.get(), 0);
}
