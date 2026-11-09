use super::*;
use crate::rans_card::RansCardDraft;

fn solver() -> RansSolver {
    let card = RansCardDraft::launder_sharma_channel("transport-regression")
        .freeze()
        .unwrap();
    let mut solver = RansSolver::new(card);
    solver.nu_laminar = 0.015;
    solver
}

fn grid() -> RansChannelGrid {
    RansChannelGrid {
        num_cells: 6,
        half_height_m: 1.0,
        y_nodes: vec![0.0, 0.03, 0.08, 0.18, 0.37, 0.66, 1.0],
    }
}

fn seeded(grid: &RansChannelGrid) -> RansFieldState {
    let shape: Vec<_> = grid
        .y_nodes
        .iter()
        .map(|&y| y * (2.0 * grid.half_height_m - y))
        .collect();
    RansFieldState {
        u: shape.iter().map(|&s| 1.1 * s).collect(),
        k: shape.iter().map(|&s| (0.32 * s).powi(2)).collect(),
        eps: shape.iter().map(|&s| 0.06 * s).collect(),
        nu_t: vec![0.0; shape.len()],
        temp_k: vec![300.0; shape.len()],
    }
}

fn bits(state: &RansFieldState) -> [Vec<u64>; 5] {
    [&state.u, &state.k, &state.eps, &state.nu_t, &state.temp_k]
        .map(|values| values.iter().map(|value| value.to_bits()).collect())
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

// Trapezoidal integration of nodal data gives the actual dual-cell storage,
// including the centerline half-cell. This oracle does not call volume().
fn integral(grid: &RansChannelGrid, field: &[f64]) -> f64 {
    grid.y_nodes
        .windows(2)
        .zip(field.windows(2))
        .map(|(y, f)| 0.5 * (y[1] - y[0]) * (f[0] + f[1]))
        .sum()
}

fn viscosities(solver: &RansSolver, state: &RansFieldState) -> Vec<f64> {
    state
        .k
        .iter()
        .zip(&state.eps)
        .map(|(&k, &eps)| {
            if k == 0.0 {
                0.0
            } else {
                let re_t = k * k / (solver.nu_laminar * eps);
                solver.coeffs.c_mu * (-3.4 / (1.0 + re_t / 50.0).powi(2)).exp() * k * k / eps
            }
        })
        .collect()
}

// An independent dense elimination oracle for the small momentum problem.
// Assemble the symmetric integrated face operator, without using any native
// transport row, diffusion, volume or tridiagonal-solve helper.
fn dense_momentum(
    grid: &RansChannelGrid,
    old: &[f64],
    nu: f64,
    acceleration: f64,
    dt: f64,
) -> Vec<f64> {
    let n = grid.y_nodes.len();
    let mut storage = vec![0.0; n];
    let mut matrix = vec![vec![0.0; n - 1]; n - 1];
    for i in 0..n - 1 {
        let width = grid.y_nodes[i + 1] - grid.y_nodes[i];
        storage[i] += 0.5 * width;
        storage[i + 1] += 0.5 * width;
        let conductance = nu / width;
        matrix[i][i] += conductance;
        if i > 0 {
            matrix[i - 1][i - 1] += conductance;
            matrix[i - 1][i] -= conductance;
            matrix[i][i - 1] -= conductance;
        }
    }
    let mut rhs = Vec::with_capacity(n - 1);
    for i in 1..n {
        matrix[i - 1][i - 1] += storage[i] / dt;
        rhs.push(storage[i] * (old[i] / dt + acceleration));
    }
    for col in 0..n - 1 {
        let pivot = (col..n - 1)
            .max_by(|&a, &b| matrix[a][col].abs().total_cmp(&matrix[b][col].abs()))
            .unwrap();
        matrix.swap(col, pivot);
        rhs.swap(col, pivot);
        let diagonal = matrix[col][col];
        for value in &mut matrix[col][col..] {
            *value /= diagonal;
        }
        rhs[col] /= diagonal;
        let pivot_row = matrix[col].clone();
        let pivot_rhs = rhs[col];
        for row in 0..n - 1 {
            if row != col {
                let factor = matrix[row][col];
                for (entry, &pivot_value) in matrix[row][col..].iter_mut().zip(&pivot_row[col..]) {
                    *entry -= factor * pivot_value;
                }
                rhs[row] -= factor * pivot_rhs;
            }
        }
    }
    let mut result = vec![0.0];
    result.extend(rhs);
    result
}

#[test]
fn laminar_invariant_and_momentum_match_independent_dense_solve() {
    let grid = grid();
    let solver = solver();
    let mut state = seeded(&grid);
    state.u.fill(0.0);
    state.k.fill(0.0);
    state.eps.fill(0.0);
    let temperature = state.temp_k.clone();
    let dt = 4.0;
    let pressure = -0.48;
    let expected = dense_momentum(
        &grid,
        &state.u,
        solver.nu_laminar,
        -pressure / solver.rho,
        dt,
    );
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let report = solver
        .advance_launder_sharma(&grid, &mut state, pressure, dt, &mut work, || false)
        .unwrap();
    for (actual, expected) in state.u.iter().zip(expected) {
        close(*actual, expected, 3e-14);
    }
    assert_eq!(state.u[0], 0.0);
    assert!(
        state.u[6] > state.u[5],
        "centerline must have its own solved half-cell"
    );
    assert!(
        state
            .k
            .iter()
            .chain(&state.eps)
            .chain(&state.nu_t)
            .all(|&v| v == 0.0)
    );
    assert_eq!(state.temp_k, temperature);
    assert_eq!(report.time_step_s, dt);
    assert_eq!(report.turbulence_production, 0.0);
    assert_eq!(report.turbulent_energy_after, 0.0);
    assert_eq!(report.turbulent_destruction, 0.0);
}

#[test]
fn large_implicit_steps_transport_positive_fields_and_recompute_viscosity() {
    let grid = grid();
    let solver = solver();
    let initial = seeded(&grid);
    let mut state = initial.clone();
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let addresses = [
        work.velocity.as_ptr(),
        work.kinetic.as_ptr(),
        work.dissipation.as_ptr(),
    ];
    // This exceeds the explicit near-wall molecular diffusion limit by
    // orders of magnitude. Positivity comes from the implicit operator.
    let dt = 20.0;
    assert!(dt * solver.nu_laminar / grid.y_nodes[1].powi(2) > 100.0);
    solver
        .advance_launder_sharma(&grid, &mut state, -0.48, dt, &mut work, || false)
        .unwrap();
    assert_ne!(state.k, initial.k);
    assert_ne!(state.eps, initial.eps);
    assert_ne!(state.u, initial.u);
    for field in [&state.k, &state.eps, &state.nu_t] {
        assert_eq!(field[0], 0.0);
        assert!(field[1..].iter().all(|v| v.is_finite() && *v > 0.0));
    }
    for (actual, expected) in state.nu_t.iter().zip(viscosities(&solver, &state)) {
        close(*actual, expected, 3e-14);
    }
    let first_viscosity = state.nu_t.clone();
    solver
        .advance_launder_sharma(&grid, &mut state, -0.48, 0.1, &mut work, || false)
        .unwrap();
    assert_ne!(state.nu_t, first_viscosity);
    assert_eq!(
        addresses,
        [
            work.velocity.as_ptr(),
            work.kinetic.as_ptr(),
            work.dissipation.as_ptr()
        ]
    );
    assert_eq!(state.temp_k, initial.temp_k);
}

#[test]
fn low_re_wall_terms_recover_quadratic_derivatives_on_stretched_grid() {
    let grid = grid();
    let solver = solver();
    let state = seeded(&grid);
    let nut = viscosities(&solver, &state);
    for i in 1..=grid.num_cells {
        let (d, e) = low_re_terms(&grid, &state.k, &state.u, solver.nu_laminar, &nut, i);
        let expected_d = 8.0
            * solver.nu_laminar
            * 0.32_f64.powi(2)
            * (grid.half_height_m - grid.y_nodes[i]).powi(2);
        let expected_e = 8.0 * solver.nu_laminar * nut[i] * 1.1_f64.powi(2);
        close(d, expected_d, 3e-13);
        close(e, expected_e, 3e-13);
        if i == grid.num_cells {
            assert_eq!(d, 0.0);
            assert!(
                e > 0.0,
                "centerline symmetry must not erase curvature production"
            );
        }
    }
    assert_eq!(
        low_re_terms(&grid, &state.k, &state.u, solver.nu_laminar, &nut, 0),
        (0.0, 0.0)
    );
}

#[test]
fn shared_face_production_closes_independent_mean_and_turbulence_budgets() {
    let grid = grid();
    let solver = solver();
    let before = seeded(&grid);
    let old_nut = viscosities(&solver, &before);
    let mut after = before.clone();
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let dt = 0.2;
    let pressure = -0.48;
    let report = solver
        .advance_launder_sharma(&grid, &mut after, pressure, dt, &mut work, || false)
        .unwrap();
    let mean_energy =
        |u: &[f64]| integral(&grid, &u.iter().map(|v| 0.5 * v * v).collect::<Vec<_>>());
    let mean_before = mean_energy(&before.u);
    let mean_after = mean_energy(&after.u);
    let k_before = integral(&grid, &before.k);
    let k_after = integral(&grid, &after.k);
    let pressure_power = -pressure / solver.rho * integral(&grid, &after.u);
    let mut molecular = 0.0;
    let mut production = 0.0;
    let mut volumes = vec![0.0; after.u.len()];
    let mut production_power = vec![0.0; after.u.len()];
    let mut epsilon_flux = vec![0.0; grid.num_cells];
    for i in 0..grid.num_cells {
        let dy = grid.y_nodes[i + 1] - grid.y_nodes[i];
        let gradient = (after.u[i + 1] - after.u[i]) / dy;
        let nut_face = 0.5 * (old_nut[i] + old_nut[i + 1]);
        let face_power = nut_face * gradient.powi(2) * dy;
        molecular += solver.nu_laminar * gradient.powi(2) * dy;
        production += face_power;
        volumes[i] += 0.5 * dy;
        volumes[i + 1] += 0.5 * dy;
        if i == 0 {
            production_power[1] += face_power;
        } else {
            production_power[i] += 0.5 * face_power;
            production_power[i + 1] += 0.5 * face_power;
        }
        epsilon_flux[i] = (solver.nu_laminar + nut_face / solver.coeffs.sigma_eps)
            * (after.eps[i + 1] - after.eps[i])
            / dy;
    }

    // Check the actual transported epsilon equation, including additive E.
    // Differentiating the local Lagrange interpolant supplies an independent
    // curvature formula; no production assembly or low_re_terms helper is used.
    let mut observes_extra_source = false;
    for i in 1..after.eps.len() {
        let a = grid.y_nodes[i] - grid.y_nodes[i - 1];
        let curvature = if i == grid.num_cells {
            2.0 * (after.u[i - 1] - after.u[i]) / a.powi(2)
        } else {
            let b = grid.y_nodes[i + 1] - grid.y_nodes[i];
            2.0 * (after.u[i - 1] / (a * (a + b)) - after.u[i] / (a * b)
                + after.u[i + 1] / (b * (a + b)))
        };
        let extra = 2.0 * solver.nu_laminar * old_nut[i] * curvature.powi(2);
        let ret = before.k[i].powi(2) / (solver.nu_laminar * before.eps[i]);
        let f2 = 1.0 - 0.3 * (-ret.powi(2)).exp();
        let rate = before.eps[i] / before.k[i];
        let source = solver.coeffs.c_eps_1 * rate * production_power[i] / volumes[i];
        let sink = solver.coeffs.c_eps_2 * f2 * rate * after.eps[i];
        let right_flux = if i == grid.num_cells {
            0.0
        } else {
            epsilon_flux[i]
        };
        let diffusion = (right_flux - epsilon_flux[i - 1]) / volumes[i];
        let change = (after.eps[i] - before.eps[i]) / dt;
        let scale = change.abs() + diffusion.abs() + source + sink + extra;
        assert!(
            (change - (diffusion + source - sink + extra)).abs() < 5e-12 * scale.max(1e-12),
            "epsilon row {i}"
        );
        observes_extra_source |= extra > 1e-7 * scale;
    }
    assert!(
        observes_extra_source,
        "the additive low-Re source must be observable"
    );
    let temporal = integral(
        &grid,
        &after
            .u
            .iter()
            .zip(&before.u)
            .map(|(new, old)| 0.5 * (new - old).powi(2) / dt)
            .collect::<Vec<_>>(),
    );
    let destruction = integral(
        &grid,
        &(0..after.k.len())
            .map(|i| {
                if i == 0 {
                    return 0.0;
                }
                // Analytic derivative of the authored quadratic sqrt(k) profile.
                let d = 8.0
                    * solver.nu_laminar
                    * 0.32_f64.powi(2)
                    * (grid.half_height_m - grid.y_nodes[i]).powi(2);
                (before.eps[i] + d) * after.k[i] / before.k[i]
            })
            .collect::<Vec<_>>(),
    );
    let wall_loss = (solver.nu_laminar + old_nut[1] / (2.0 * solver.coeffs.sigma_k)) * after.k[1]
        / grid.y_nodes[1];
    for (actual, expected) in [
        (report.mean_energy_before, mean_before),
        (report.mean_energy_after, mean_after),
        (report.turbulent_energy_before, k_before),
        (report.turbulent_energy_after, k_after),
        (report.pressure_power, pressure_power),
        (report.molecular_dissipation, molecular),
        (report.turbulence_production, production),
        (report.temporal_dissipation, temporal),
        (report.turbulent_destruction, destruction),
        (report.turbulent_wall_loss, wall_loss),
    ] {
        close(actual, expected, 5e-13);
    }
    close(
        (mean_after - mean_before) / dt,
        pressure_power - molecular - production - temporal,
        5e-13,
    );
    close(
        (k_after - k_before) / dt,
        production - destruction - wall_loss,
        5e-13,
    );
    let total_change = (mean_after + k_after - mean_before - k_before) / dt;
    close(
        total_change,
        pressure_power - molecular - temporal - destruction - wall_loss,
        8e-13,
    );
    assert!(production > 0.0 && destruction > 0.0 && wall_loss > 0.0);
    assert!(report.mean_balance_defect.abs() < 2e-14);
    assert!(report.turbulent_balance_defect.abs() < 2e-14);
}

#[test]
fn unforced_turbulence_loses_energy_without_floors_or_automatic_seeding() {
    let grid = grid();
    let solver = solver();
    let mut state = seeded(&grid);
    state.u.fill(0.0);
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    for _ in 0..20 {
        let old_k = integral(&grid, &state.k);
        let old_e = integral(&grid, &state.eps);
        let report = solver
            .advance_launder_sharma(&grid, &mut state, 0.0, 0.1, &mut work, || false)
            .unwrap();
        assert!(integral(&grid, &state.k) < old_k);
        assert!(integral(&grid, &state.eps) < old_e);
        assert!(state.k[1..].iter().chain(&state.eps[1..]).all(|&v| v > 0.0));
        assert!(state.u.iter().all(|&v| v == 0.0));
        assert_eq!(report.turbulence_production, 0.0);
        assert_eq!(report.mean_energy_after, 0.0);
        assert!(report.turbulent_destruction > 0.0);
    }
}

#[test]
fn flow_reversal_and_dimensional_scaling_preserve_the_same_physical_step() {
    let grid = grid();
    let solver = solver();
    let initial = seeded(&grid);
    let mut baseline = initial.clone();
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let report = solver
        .advance_launder_sharma(&grid, &mut baseline, -0.48, 0.2, &mut work, || false)
        .unwrap();
    let mut reverse = initial.clone();
    for u in &mut reverse.u {
        *u = -*u;
    }
    solver
        .advance_launder_sharma(&grid, &mut reverse, 0.48, 0.2, &mut work, || false)
        .unwrap();
    for (&forward, &reversed) in baseline.u.iter().zip(&reverse.u) {
        assert_eq!(forward, -reversed);
    }
    assert_eq!(baseline.k, reverse.k);
    assert_eq!(baseline.eps, reverse.eps);
    assert_eq!(baseline.nu_t, reverse.nu_t);

    for (length, velocity) in [(7.0, 0.3), (0.0005, 20.0)] {
        let scaled_grid = RansChannelGrid {
            num_cells: grid.num_cells,
            half_height_m: grid.half_height_m * length,
            y_nodes: grid.y_nodes.iter().map(|y| y * length).collect(),
        };
        let mut scaled_solver = self::solver();
        scaled_solver.nu_laminar *= length * velocity;
        let mut state = initial.clone();
        let factors = [
            velocity,
            velocity * velocity,
            velocity.powi(3) / length,
            length * velocity,
        ];
        for (field, factor) in [&mut state.u, &mut state.k, &mut state.eps, &mut state.nu_t]
            .into_iter()
            .zip(factors)
        {
            for value in field {
                *value *= factor;
            }
        }
        let mut scaled_work = RansTransportWorkspace::try_new(&scaled_grid, 7).unwrap();
        let scaled_report = scaled_solver
            .advance_launder_sharma(
                &scaled_grid,
                &mut state,
                -0.48 * velocity.powi(2) / length,
                0.2 * length / velocity,
                &mut scaled_work,
                || false,
            )
            .unwrap();
        for ((actual, expected), factor) in [&state.u, &state.k, &state.eps, &state.nu_t]
            .into_iter()
            .zip([&baseline.u, &baseline.k, &baseline.eps, &baseline.nu_t])
            .zip(factors)
        {
            for (a, b) in actual.iter().zip(expected) {
                close(*a, *b * factor, 2e-12);
            }
        }
        close(
            scaled_report.mean_energy_after,
            report.mean_energy_after * length * velocity.powi(2),
            2e-12,
        );
        close(
            scaled_report.turbulent_energy_after,
            report.turbulent_energy_after * length * velocity.powi(2),
            2e-12,
        );
        close(
            scaled_report.pressure_power,
            report.pressure_power * velocity.powi(3),
            2e-12,
        );
        close(
            scaled_report.turbulence_production,
            report.turbulence_production * velocity.powi(3),
            2e-12,
        );
        assert_eq!(state.temp_k, initial.temp_k);
    }
}

#[test]
fn budgets_invalid_inputs_and_late_arithmetic_failure_preserve_every_field() {
    let grid = grid();
    let solver = solver();
    assert!(matches!(
        RansTransportWorkspace::try_new(&grid, 6),
        Err(RansTransportError::NodeBudget {
            needed: 7,
            maximum: 6
        })
    ));
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    assert_eq!(work.numerical_bytes(), 7 * 9 * core::mem::size_of::<f64>());
    for case in 0..10 {
        let mut state = seeded(&grid);
        let mut bad_grid = grid.clone();
        let mut dt = 0.2;
        let mut pressure = -0.48;
        match case {
            0 => state.u[3] = f64::NAN,
            1 => {
                state.k[3] = 0.0;
                state.eps[3] = 0.0;
            }
            2 => state.eps[3] = 0.0,
            3 => state.k[0] = 1e-4,
            4 => state.k[2] = -1.0,
            5 => state.temp_k[2] = 0.0,
            6 => bad_grid.y_nodes[3] = bad_grid.y_nodes[2],
            7 => bad_grid.half_height_m = 1.1,
            8 => dt = 0.0,
            9 => pressure = f64::INFINITY,
            _ => unreachable!(),
        }
        let before = bits(&state);
        assert!(
            matches!(
                solver.advance_launder_sharma(
                    &bad_grid,
                    &mut state,
                    pressure,
                    dt,
                    &mut work,
                    || false,
                ),
                Err(RansTransportError::InvalidInput(_))
            ),
            "case {case}"
        );
        assert_eq!(bits(&state), before, "case {case} changed state");
    }
    let mut state = seeded(&grid);
    let before = bits(&state);
    // Momentum remains representable, but its subsequent squared face work
    // overflows. A computed velocity trial must not leak into accepted state.
    assert!(matches!(
        solver.advance_launder_sharma(&grid, &mut state, -1e200, 0.2, &mut work, || false,),
        Err(RansTransportError::Arithmetic { .. })
    ));
    assert_eq!(bits(&state), before);
}

#[test]
fn every_cancellation_poll_is_atomic_and_workspace_can_be_reused() {
    let grid = RansChannelGrid::new_stretched(132, 1.0, 1.015);
    let solver = solver();
    let initial = seeded(&grid);
    let n = grid.num_cells + 1;
    let mut work = RansTransportWorkspace::try_new(&grid, n).unwrap();
    let mut expected = initial.clone();
    let mut polls = 0;
    let expected_report = solver
        .advance_launder_sharma(&grid, &mut expected, -0.48, 0.2, &mut work, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(
        polls > 20,
        "exercise assembly, elimination, diagnostics and final publication"
    );
    for fail_at in 0..polls {
        let mut state = initial.clone();
        let mut seen = 0;
        let result =
            solver.advance_launder_sharma(&grid, &mut state, -0.48, 0.2, &mut work, || {
                let stop = seen == fail_at;
                seen += 1;
                stop
            });
        assert_eq!(result.unwrap_err(), RansTransportError::Cancelled);
        assert_eq!(
            bits(&state),
            bits(&initial),
            "poll {fail_at} published a trial"
        );
        let report = solver
            .advance_launder_sharma(&grid, &mut state, -0.48, 0.2, &mut work, || false)
            .unwrap();
        assert_eq!(bits(&state), bits(&expected));
        assert_eq!(report, expected_report);
    }
}
