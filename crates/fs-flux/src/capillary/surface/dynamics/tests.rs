use super::*;

fn controls() -> SurfaceDynamicsControls {
    let mut solve = ImplicitSolveConfig::default();
    solve.max_newton_iterations = 32;
    solve.newton.absolute_tolerance = 1e-13;
    solve.newton.relative_tolerance = 1e-12;
    SurfaceDynamicsControls {
        solve,
        maximum_krylov_restart: 32,
        impulse_tolerance_n_s: 1e-11,
        kinematic_tolerance_m: 1e-12,
        energy_tolerance_j: 1e-11,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

fn positions(state: &SecondOrderState) -> Vec<Vector> {
    state
        .q
        .chunks_exact(3)
        .map(|v| [v[0], v[1], v[2]])
        .collect()
}

fn momentum(state: &SecondOrderState, masses: &[f64]) -> (Vector, Vector) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for ((q, v), mass) in state
        .q
        .chunks_exact(3)
        .zip(state.v.chunks_exact(3))
        .zip(masses)
    {
        let p = [mass * v[0], mass * v[1], mass * v[2]];
        let l = cross([q[0], q[1], q[2]], p);
        for d in 0..3 {
            linear[d] += p[d];
            angular[d] += l[d];
        }
    }
    (linear, angular)
}

fn tetrahedron() -> ([Vector; 4], [[usize; 3]; 4], [f64; 4], [Vector; 4]) {
    let q = [
        [1.1, 0.1, -0.1],
        [-0.4, 0.9, 0.2],
        [-0.5, -0.5, 0.4],
        [0.1, -0.2, 1.3],
    ];
    let faces = [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]];
    let masses = [0.7, 1.1, 1.3, 0.9];
    let v = core::array::from_fn(|i| {
        let spin = cross([0.7, -0.4, 1.8], q[i]);
        [
            spin[0] + 0.2 * q[i][1],
            spin[1] - 0.1 * q[i][2],
            spin[2] + 0.15 * q[i][0],
        ]
    });
    (q, faces, masses, v)
}

#[test]
fn g1_equilateral_breathing_matches_midpoint_solution_and_second_order_time_refinement() {
    let q = [
        [1.0, 0.0, 0.0],
        [-0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        [-0.5, -3.0_f64.sqrt() / 2.0, 0.0],
    ];
    let mass = [0.8; 3];
    let faces = [[0, 1, 2]];
    let sigma = 0.6;
    let field = SurfaceTensionField::Uniform(sigma);
    let surface = CapillarySurfaceDynamics::try_new(&mass, &faces, 3, 1).unwrap();
    let velocity = q.map(|p| p.map(|v| 0.2 * v));
    let force = [[0.0; 3]; 3];
    let initial = surface
        .initial_state(0.0, &q, &velocity, field, &force)
        .unwrap();
    // A=(3 sqrt(3)/4) r² and three masses give r''+omega² r=0.
    // The original structural midpoint step is its exact Cayley map.
    let square = 3.0_f64.sqrt() * sigma / (2.0 * mass[0]);
    let h = 0.08;
    let factor = 0.25 * h * h * square;
    let expected_r = ((1.0 - factor) + h * 0.2) / (1.0 + factor);
    let expected_v = ((1.0 - factor) * 0.2 - h * square) / (1.0 + factor);
    let mut state = initial.clone();
    let report = surface
        .advance(&mut state, field, field, &force, h, controls(), || false)
        .unwrap();
    for i in 0..3 {
        for d in 0..3 {
            close(state.q[3 * i + d], expected_r * q[i][d], 2e-13);
            close(state.v[3 * i + d], expected_v * q[i][d], 2e-12);
        }
    }
    close(
        report.kinetic_change_j + report.surface.energy_change_j,
        0.0,
        2e-13,
    );
    assert_eq!(state.steps, 1);
    assert_eq!(state.history.len(), 1);
    let time = 0.4;
    let omega = square.sqrt();
    let exact = (omega * time).cos() + 0.2 * (omega * time).sin() / omega;
    let errors: Vec<_> = [8, 16, 32]
        .map(|steps| {
            let mut state = initial.clone();
            for _ in 0..steps {
                surface
                    .advance(
                        &mut state,
                        field,
                        field,
                        &force,
                        time / steps as f64,
                        controls(),
                        || false,
                    )
                    .unwrap();
            }
            assert_eq!(state.history.len(), steps);
            (state.q[0] - exact).abs()
        })
        .into();
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.9 && pair[0] / pair[1] < 4.1,
            "{errors:?}"
        );
    }
}

#[test]
fn g1_rotating_nonuniform_surface_preserves_actual_energy_linear_and_angular_momentum() {
    let (q, faces, masses, v) = tetrahedron();
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let sigma = [0.2, 0.3, 0.5, 0.4];
    let field = SurfaceTensionField::PerFace(&sigma);
    let force = [[0.0; 3]; 4];
    let mut state = surface.initial_state(0.0, &q, &v, field, &force).unwrap();
    let (linear, angular) = momentum(&state, &masses);
    let initial_surface = assemble_surface_tension_step(&q, &q, &faces, field, field).unwrap();
    let initial_energy =
        kinetic(&masses, &state.v).unwrap() + initial_surface.energy.energy_before_j;
    for _ in 0..40 {
        let report = surface
            .advance(&mut state, field, field, &force, 0.005, controls(), || {
                false
            })
            .unwrap();
        close(
            report.kinetic_after_j + report.surface.energy_after_j,
            initial_energy,
            3e-11,
        );
        assert_eq!(report.surface.coefficient_supply_j, 0.0);
        assert_eq!(report.external_work_j, 0.0);
        let (p, l) = momentum(&state, &masses);
        for d in 0..3 {
            close(p[d], linear[d], 2e-11);
            close(l[d], angular[d], 3e-11);
        }
    }
    let edges = |p: &[Vector]| {
        (0..3)
            .map(|i| dot(sub(p[i + 1], p[0]), sub(p[i + 1], p[0])))
            .collect::<Vec<_>>()
    };
    let before = edges(&q);
    let after = edges(&positions(&state));
    assert!(
        before.iter().zip(after).any(|(a, b)| (a - b).abs() > 0.01),
        "the fixture must actually deform, not merely test rigid transport"
    );
}

#[test]
fn g1_prescribed_sigma_change_supplies_free_energy_and_held_force_equilibrium_is_stationary() {
    let q = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let faces = [[0, 1, 2]];
    let masses = [1.0; 3];
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 3, 1).unwrap();
    let zero = SurfaceTensionField::Uniform(0.0);
    let sigma = SurfaceTensionField::Uniform(0.6);
    let no_force = [[0.0; 3]; 3];
    let mut state = surface
        .initial_state(0.0, &q, &no_force, zero, &no_force)
        .unwrap();
    let report = surface
        .advance(&mut state, zero, sigma, &no_force, 0.1, controls(), || {
            false
        })
        .unwrap();
    assert!(report.kinetic_after_j > 0.0);
    assert!(report.surface.area_after_m2 < report.surface.area_before_m2);
    close(
        report.surface.coefficient_supply_j,
        0.3 * (report.surface.area_before_m2 + report.surface.area_after_m2),
        1e-14,
    );
    close(
        report.kinetic_change_j + report.surface.energy_change_j,
        report.surface.coefficient_supply_j,
        2e-12,
    );

    let initial_force = triangle_surface_tension(q, 0.6).unwrap().forces_n;
    let balanced = initial_force.map(|f| f.map(|v| -v));
    let mut equilibrium = surface
        .initial_state(0.0, &q, &no_force, sigma, &balanced)
        .unwrap();
    let initial = equilibrium.clone();
    let report = surface
        .advance(
            &mut equilibrium,
            sigma,
            sigma,
            &balanced,
            0.3,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(equilibrium.q, initial.q);
    assert_eq!(equilibrium.v, initial.v);
    assert_eq!(report.external_work_j, 0.0);
    assert_eq!(report.kinetic_change_j, 0.0);

    // The zero-tension limit is actual uniformly accelerated nodal inertia.
    // Nonzero exterior work must enter once with the correct sign and units.
    let load = [[0.3, -0.2, 0.1]; 3];
    let mut ballistic = surface
        .initial_state(0.0, &q, &no_force, zero, &load)
        .unwrap();
    let h = 0.2;
    let report = surface
        .advance(&mut ballistic, zero, zero, &load, h, controls(), || false)
        .unwrap();
    for i in 0..3 {
        for d in 0..3 {
            close(
                ballistic.q[3 * i + d],
                q[i][d] + 0.5 * h * h * load[i][d] / masses[i],
                1e-14,
            );
            close(ballistic.v[3 * i + d], h * load[i][d] / masses[i], 1e-13);
        }
    }
    assert!(report.external_work_j > 0.0);
    close(report.kinetic_change_j, report.external_work_j, 2e-14);
    assert_eq!(report.surface.energy_change_j, 0.0);
}

#[test]
fn g3_surface_motion_is_objective_under_finite_rotation_translation_and_checkpoint_replay() {
    let (q, faces, masses, v) = tetrahedron();
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let rotate = |p: Vector| [p[2], p[1], -p[0]];
    let transformed = q.map(|p| {
        let p = rotate(p);
        [p[0] + 2.0, p[1] - 1.0, p[2] + 3.0]
    });
    let sigma0 = [0.2, 0.3, 0.5, 0.4];
    let sigma1 = [0.25, 0.28, 0.6, 0.35];
    let a = SurfaceTensionField::PerFace(&sigma0);
    let b = SurfaceTensionField::PerFace(&sigma1);
    let force = [[0.0; 3]; 4];
    let mut original = surface.initial_state(0.0, &q, &v, a, &force).unwrap();
    let mut replay = original.clone();
    let mut rotated = surface
        .initial_state(0.0, &transformed, &v.map(rotate), a, &force)
        .unwrap();
    let report = surface
        .advance(&mut original, a, b, &force, 0.03, controls(), || false)
        .unwrap();
    let other = surface
        .advance(&mut rotated, a, b, &force, 0.03, controls(), || false)
        .unwrap();
    surface
        .advance(&mut replay, a, b, &force, 0.03, controls(), || false)
        .unwrap();
    assert_eq!(original, replay);
    for i in 0..4 {
        let p = rotate([
            original.q[3 * i],
            original.q[3 * i + 1],
            original.q[3 * i + 2],
        ]);
        let v = rotate([
            original.v[3 * i],
            original.v[3 * i + 1],
            original.v[3 * i + 2],
        ]);
        for d in 0..3 {
            close(rotated.q[3 * i + d], p[d] + [2.0, -1.0, 3.0][d], 2e-12);
            close(rotated.v[3 * i + d], v[d], 2e-11);
        }
    }
    close(
        report.surface.energy_after_j,
        other.surface.energy_after_j,
        2e-12,
    );
    close(
        report.surface.coefficient_supply_j,
        other.surface.coefficient_supply_j,
        2e-13,
    );
    close(report.energy_residual_j, 0.0, 2e-12);
}

#[test]
fn g4_cancelled_invalid_and_unconverged_steps_preserve_native_history_and_every_state_field() {
    let (q, faces, masses, v) = tetrahedron();
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let field = SurfaceTensionField::Uniform(0.4);
    let force = [[0.0; 3]; 4];
    let mut state = surface.initial_state(0.0, &q, &v, field, &force).unwrap();
    surface
        .advance(&mut state, field, field, &force, 0.01, controls(), || false)
        .unwrap();
    let checkpoint = state.clone();
    let mut calls = 0;
    let error = surface
        .advance(&mut state, field, field, &force, 0.01, controls(), || {
            calls += 1;
            calls >= 9
        })
        .unwrap_err();
    assert!(matches!(error, SurfaceDynamicsError::Cancelled));
    assert_eq!(state, checkpoint);
    let bad = SurfaceTensionField::Uniform(-0.1);
    assert!(matches!(
        surface.advance(&mut state, field, bad, &force, 0.01, controls(), || false),
        Err(SurfaceDynamicsError::Surface(
            CapillarySurfaceError::InvalidSurfaceTension
        ))
    ));
    assert_eq!(state, checkpoint);
    let mut budget = controls();
    budget.solve.max_newton_iterations = 1;
    assert!(matches!(
        surface.advance(&mut state, field, field, &force, 0.3, budget, || false),
        Err(SurfaceDynamicsError::Time(_))
    ));
    assert_eq!(state, checkpoint);
    let mut exact_gate = controls();
    exact_gate.impulse_tolerance_n_s = 0.0;
    exact_gate.kinematic_tolerance_m = 0.0;
    exact_gate.energy_tolerance_j = 0.0;
    assert!(matches!(
        surface.advance(&mut state, field, field, &force, 0.01, exact_gate, || false),
        Err(SurfaceDynamicsError::Balance { .. })
    ));
    assert_eq!(state, checkpoint);
    assert!(CapillarySurfaceDynamics::try_new(&[1.0, 0.0, 1.0], &[[0, 1, 2]], 3, 1).is_err());
    assert!(CapillarySurfaceDynamics::try_new(&masses, &faces, 3, 4).is_err());
}
