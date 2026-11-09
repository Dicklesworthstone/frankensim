use super::*;
use crate::rans_card::RansCardDraft;

mod conjugate;
mod dissipative;

fn solver() -> RansSolver {
    let card = RansCardDraft::launder_sharma_channel("thermal-transport-regression")
        .freeze()
        .unwrap();
    let mut solver = RansSolver::new(card);
    solver.nu_laminar = 0.015;
    solver.rho = 2.0;
    solver.cp = 3.0;
    solver.k_thermal = 0.06;
    solver.pr_turbulent = 0.7;
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
        temp_k: shape.iter().map(|&s| 300.0 + 12.0 * s).collect(),
    }
}

fn laminar(grid: &RansChannelGrid, temperature: f64) -> RansFieldState {
    let mut state = seeded(grid);
    state.u.fill(0.0);
    state.k.fill(0.0);
    state.eps.fill(0.0);
    state.temp_k.fill(temperature);
    state
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

fn integral(grid: &RansChannelGrid, values: &[f64]) -> f64 {
    grid.y_nodes
        .windows(2)
        .zip(values.windows(2))
        .map(|(y, value)| 0.5 * (y[1] - y[0]) * (value[0] + value[1]))
        .sum()
}

fn old_viscosities(solver: &RansSolver, state: &RansFieldState) -> Vec<f64> {
    state
        .k
        .iter()
        .zip(&state.eps)
        .map(|(&k, &eps)| {
            if k == 0.0 {
                0.0
            } else {
                let ret = k * k / (solver.nu_laminar * eps);
                solver.coeffs.c_mu * (-3.4 / (1.0 + ret / 50.0).powi(2)).exp() * k * k / eps
            }
        })
        .collect()
}

// Independent heat equation in physical units: storage is rho*cp*volume
// and face conductivity is kappa+rho*cp*nu_t/Pr. Dense pivoted elimination
// shares no production assembly, volume or tridiagonal helper.
#[allow(clippy::too_many_arguments)]
fn dense_heat(
    solver: &RansSolver,
    grid: &RansChannelGrid,
    old: &[f64],
    nut: &[f64],
    wall_flux: f64,
    heating: Option<&[f64]>,
    dt: f64,
) -> Vec<f64> {
    let n = old.len();
    let capacity = solver.rho * solver.cp;
    let mut volumes = vec![0.0; n];
    let mut matrix = vec![vec![0.0; n]; n];
    for i in 0..n - 1 {
        let dy = grid.y_nodes[i + 1] - grid.y_nodes[i];
        volumes[i] += 0.5 * dy;
        volumes[i + 1] += 0.5 * dy;
        let conductance = (solver.k_thermal
            + capacity * (nut[i] + nut[i + 1]) / (2.0 * solver.pr_turbulent))
            / dy;
        matrix[i][i] += conductance;
        matrix[i + 1][i + 1] += conductance;
        matrix[i][i + 1] -= conductance;
        matrix[i + 1][i] -= conductance;
    }
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        matrix[i][i] += capacity * volumes[i] / dt;
        rhs[i] = capacity * volumes[i] * old[i] / dt
            + volumes[i] * heating.map_or(0.0, |values| values[i]);
    }
    rhs[0] += wall_flux;
    for col in 0..n {
        let pivot = (col..n)
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
        for row in 0..n {
            if row != col {
                let factor = matrix[row][col];
                for (value, pivot_value) in matrix[row][col..].iter_mut().zip(&pivot_row[col..]) {
                    *value -= factor * pivot_value;
                }
                rhs[row] -= factor * pivot_rhs;
            }
        }
    }
    rhs
}

#[test]
fn implicit_heat_matches_dense_oracle_and_preserves_the_mechanical_endpoint() {
    let solver = solver();
    let grid = grid();
    let mut initial = seeded(&grid);
    initial.temp_k = vec![305.0, 301.0, 309.0, 297.0, 302.0, 311.0, 300.0];
    let heating = [0.2, -0.1, 0.0, 0.8, 0.4, 0.2, 0.1];
    let mut mechanical = initial.clone();
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let mechanical_report = solver
        .advance_launder_sharma(&grid, &mut mechanical, -0.8, 0.3, &mut work, || false)
        .unwrap();
    let expected = dense_heat(
        &solver,
        &grid,
        &initial.temp_k,
        &mechanical.nu_t,
        0.7,
        Some(&heating),
        0.3,
    );
    let mut state = initial.clone();
    let report = solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut state,
            -0.8,
            0.3,
            RansHeatLoad {
                wall_heat_flux_w_m2: 0.7,
                volumetric_heating_w_m3: Some(&heating),
                ..RansHeatLoad::default()
            },
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(report.mechanical, mechanical_report);
    assert_eq!(&bits(&state)[..4], &bits(&mechanical)[..4]);
    for (&actual, expected) in state.temp_k.iter().zip(expected) {
        close(actual, expected, 3e-13);
    }

    let capacity = solver.rho * solver.cp;
    let before = capacity * integral(&grid, &initial.temp_k);
    let after = capacity * integral(&grid, &state.temp_k);
    let increments: Vec<_> = state
        .temp_k
        .iter()
        .zip(&initial.temp_k)
        .map(|(a, b)| a - b)
        .collect();
    let change = capacity * integral(&grid, &increments);
    let volume_power = integral(&grid, &heating);
    close(report.thermal_energy_before_j_m2, before, 3e-14);
    close(report.thermal_energy_after_j_m2, after, 3e-14);
    close(report.thermal_energy_change_j_m2, change, 3e-13);
    close(report.volumetric_heat_input_w_m2, volume_power, 3e-14);
    assert_eq!(report.wall_heat_flux_w_m2, 0.7);
    assert!((change / 0.3 - (0.7 + volume_power)).abs() < 2e-10);
    assert!(report.thermal_balance_defect_w_m2.abs() < 2e-10);
}

#[test]
fn wall_and_volume_heating_have_the_correct_inventory_and_storage_rates() {
    let grid = grid();
    let heating = vec![12.0; 7];
    for (density, heat_capacity) in [(2.0, 3.0), (4.0, 3.0), (2.0, 5.0)] {
        let mut solver = solver();
        solver.rho = density;
        solver.cp = heat_capacity;
        let mut state = laminar(&grid, 303.25);
        let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
        let report = solver
            .advance_launder_sharma_with_heat(
                &grid,
                &mut state,
                0.0,
                0.125,
                RansHeatLoad {
                    wall_heat_flux_w_m2: 0.0,
                    volumetric_heating_w_m3: Some(&heating),
                    ..RansHeatLoad::default()
                },
                &mut work,
                || false,
            )
            .unwrap();
        let expected = 303.25 + 12.0 * 0.125 / (density * heat_capacity);
        for temperature in state.temp_k {
            close(temperature, expected, 3e-14);
        }
        close(
            report.volumetric_heat_input_w_m2,
            12.0 * grid.half_height_m,
            3e-14,
        );
        close(
            report.thermal_energy_change_j_m2,
            12.0 * grid.half_height_m * 0.125,
            2e-12,
        );
    }
    let solver = solver();
    let mut state = laminar(&grid, 300.0);
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    for step in 1..=5 {
        let report = solver
            .advance_launder_sharma_with_heat(
                &grid,
                &mut state,
                0.0,
                0.3,
                RansHeatLoad {
                    wall_heat_flux_w_m2: 0.8,
                    volumetric_heating_w_m3: None,
                    ..RansHeatLoad::default()
                },
                &mut work,
                || false,
            )
            .unwrap();
        close(report.thermal_energy_change_j_m2, 0.3 * 0.8, 3e-11);
        let expected =
            300.0 * grid.half_height_m + f64::from(step) * 0.3 * 0.8 / (solver.rho * solver.cp);
        close(integral(&grid, &state.temp_k), expected, 3e-14);
        assert!(state.temp_k[0] > state.temp_k[1]);
        assert!(report.thermal_balance_defect_w_m2.abs() < 2e-10);
    }
}

#[test]
fn thermal_diffusion_uses_updated_turbulence_and_the_authored_prandtl_number() {
    let grid = grid();
    let mut initial = seeded(&grid);
    initial.temp_k = vec![300.0, 315.0, 294.0, 327.0, 302.0, 319.0, 298.0];
    let mut outcomes = Vec::new();
    for prandtl in [0.4, 1.6] {
        let mut solver = solver();
        solver.pr_turbulent = prandtl;
        let old_nut = old_viscosities(&solver, &initial);
        let frozen = dense_heat(&solver, &grid, &initial.temp_k, &old_nut, 0.0, None, 1.5);
        let mut state = initial.clone();
        let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
        let report = solver
            .advance_launder_sharma_with_heat(
                &grid,
                &mut state,
                -0.8,
                1.5,
                RansHeatLoad {
                    wall_heat_flux_w_m2: 0.0,
                    volumetric_heating_w_m3: None,
                    ..RansHeatLoad::default()
                },
                &mut work,
                || false,
            )
            .unwrap();
        let updated = dense_heat(&solver, &grid, &initial.temp_k, &state.nu_t, 0.0, None, 1.5);
        for (&actual, expected) in state.temp_k.iter().zip(updated) {
            close(actual, expected, 3e-13);
        }
        let frozen_error = state
            .temp_k
            .iter()
            .zip(frozen)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(
            frozen_error > 1e-5,
            "fixture must distinguish new and frozen eddy viscosity"
        );
        outcomes.push((state, report));
    }
    assert_eq!(&bits(&outcomes[0].0)[..4], &bits(&outcomes[1].0)[..4]);
    assert_eq!(outcomes[0].1.mechanical, outcomes[1].1.mechanical);
    let prandtl_effect = outcomes[0]
        .0
        .temp_k
        .iter()
        .zip(&outcomes[1].0.temp_k)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(prandtl_effect > 1e-4);
}

#[test]
fn no_source_heat_obeys_the_maximum_principle_and_preserves_uniform_temperature() {
    let grid = grid();
    let solver = solver();
    let mut state = seeded(&grid);
    state.temp_k = vec![295.0, 340.0, 300.0, 320.0, 297.0, 330.0, 305.0];
    let inventory = integral(&grid, &state.temp_k);
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut state,
            0.0,
            100.0,
            RansHeatLoad {
                wall_heat_flux_w_m2: 0.0,
                volumetric_heating_w_m3: None,
                ..RansHeatLoad::default()
            },
            &mut work,
            || false,
        )
        .unwrap();
    assert!(state.temp_k.iter().all(|&t| (295.0..=340.0).contains(&t)));
    close(integral(&grid, &state.temp_k), inventory, 3e-13);

    let mut uniform = seeded(&grid);
    uniform.temp_k.fill(300.125);
    let temperature_bits = bits(&uniform)[4].clone();
    let old_velocity = uniform.u.clone();
    let report = solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut uniform,
            -2.0,
            20.0,
            RansHeatLoad {
                wall_heat_flux_w_m2: 0.0,
                volumetric_heating_w_m3: None,
                ..RansHeatLoad::default()
            },
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(bits(&uniform)[4], temperature_bits);
    assert_ne!(uniform.u, old_velocity);
    assert!(report.mechanical.molecular_dissipation > 0.0);
    assert_eq!(report.thermal_energy_change_j_m2, 0.0);
    assert_eq!(report.thermal_balance_defect_w_m2, 0.0);
}

#[test]
fn signed_cooling_is_admitted_until_temperature_would_be_nonpositive() {
    let grid = grid();
    let solver = solver();
    let cooling = vec![-30.0; 7];
    let mut state = laminar(&grid, 3.0);
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let report = solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut state,
            0.0,
            0.2,
            RansHeatLoad {
                wall_heat_flux_w_m2: 0.0,
                volumetric_heating_w_m3: Some(&cooling),
                ..RansHeatLoad::default()
            },
            &mut work,
            || false,
        )
        .unwrap();
    for &temperature in &state.temp_k {
        close(temperature, 2.0, 3e-14);
    }
    close(report.thermal_energy_change_j_m2, -6.0, 3e-14);

    let mut state = seeded(&grid);
    state.temp_k.fill(3.0);
    let before = bits(&state);
    let severe_cooling = vec![-6000.0; 7];
    assert!(matches!(
        solver.advance_launder_sharma_with_heat(
            &grid,
            &mut state,
            -0.8,
            1.0,
            RansHeatLoad {
                wall_heat_flux_w_m2: 0.0,
                volumetric_heating_w_m3: Some(&severe_cooling),
                ..RansHeatLoad::default()
            },
            &mut work,
            || false,
        ),
        Err(RansTransportError::Arithmetic { .. })
    ));
    assert_eq!(bits(&state), before);
}

#[test]
fn thermal_admission_and_late_overflow_preserve_the_complete_accepted_state() {
    let grid = grid();
    let mut work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    for case in 0..7 {
        let mut solver = solver();
        let mut state = seeded(&grid);
        let mut wall_flux = 0.7;
        let mut sources = vec![0.2; 7];
        match case {
            0 => solver.cp = 0.0,
            1 => solver.k_thermal = f64::INFINITY,
            2 => solver.pr_turbulent = 0.0,
            3 => sources.truncate(6),
            4 => sources[3] = f64::NAN,
            5 => wall_flux = f64::INFINITY,
            6 => solver.rho = 0.0,
            _ => unreachable!(),
        }
        let before = bits(&state);
        assert!(
            matches!(
                solver.advance_launder_sharma_with_heat(
                    &grid,
                    &mut state,
                    -0.8,
                    0.3,
                    RansHeatLoad {
                        wall_heat_flux_w_m2: wall_flux,
                        volumetric_heating_w_m3: Some(&sources),
                        ..RansHeatLoad::default()
                    },
                    &mut work,
                    || false,
                ),
                Err(RansTransportError::InvalidInput(_))
            ),
            "case {case}"
        );
        assert_eq!(bits(&state), before);
    }
    let solver = solver();
    let mut state = seeded(&grid);
    let before = bits(&state);
    assert!(matches!(
        solver.advance_launder_sharma_with_heat(
            &grid,
            &mut state,
            -0.8,
            2.0,
            RansHeatLoad {
                wall_heat_flux_w_m2: f64::MAX,
                volumetric_heating_w_m3: None,
                ..RansHeatLoad::default()
            },
            &mut work,
            || false,
        ),
        Err(RansTransportError::Arithmetic { .. })
    ));
    assert_eq!(bits(&state), before);
    let mut expected = state.clone();
    let mut fresh = RansTransportWorkspace::try_new(&grid, 7).unwrap();
    let expected_report = solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut expected,
            -0.8,
            0.3,
            RansHeatLoad {
                wall_heat_flux_w_m2: 0.7,
                volumetric_heating_w_m3: None,
                ..RansHeatLoad::default()
            },
            &mut fresh,
            || false,
        )
        .unwrap();
    let report = solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut state,
            -0.8,
            0.3,
            RansHeatLoad {
                wall_heat_flux_w_m2: 0.7,
                volumetric_heating_w_m3: None,
                ..RansHeatLoad::default()
            },
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(bits(&state), bits(&expected));
    assert_eq!(report, expected_report);
}

#[test]
fn every_cancellation_poll_rolls_back_heat_and_mechanics_then_replays() {
    let grid = RansChannelGrid::new_stretched(132, 1.0, 1.015);
    let solver = solver();
    let initial = seeded(&grid);
    let n = grid.num_cells + 1;
    let heating = vec![0.8; n];
    let mut work = RansTransportWorkspace::try_new(&grid, n).unwrap();
    let pointers = [
        work.velocity.as_ptr(),
        work.kinetic.as_ptr(),
        work.dissipation.as_ptr(),
        work.viscosity.as_ptr(),
        work.production.as_ptr(),
        work.lower.as_ptr(),
        work.diagonal.as_ptr(),
        work.upper.as_ptr(),
        work.rhs.as_ptr(),
    ];
    let mut expected = initial.clone();
    let mut polls = 0;
    let expected_report = solver
        .advance_launder_sharma_with_heat(
            &grid,
            &mut expected,
            -0.8,
            0.3,
            RansHeatLoad {
                wall_heat_flux_w_m2: 0.7,
                volumetric_heating_w_m3: Some(&heating),
                ..RansHeatLoad::default()
            },
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 20);
    for fail_at in 0..polls {
        let mut state = initial.clone();
        let mut seen = 0;
        let error = solver
            .advance_launder_sharma_with_heat(
                &grid,
                &mut state,
                -0.8,
                0.3,
                RansHeatLoad {
                    wall_heat_flux_w_m2: 0.7,
                    volumetric_heating_w_m3: Some(&heating),
                    ..RansHeatLoad::default()
                },
                &mut work,
                || {
                    let stop = seen == fail_at;
                    seen += 1;
                    stop
                },
            )
            .unwrap_err();
        assert_eq!(error, RansTransportError::Cancelled);
        assert_eq!(bits(&state), bits(&initial), "cancel poll {fail_at}");
        let report = solver
            .advance_launder_sharma_with_heat(
                &grid,
                &mut state,
                -0.8,
                0.3,
                RansHeatLoad {
                    wall_heat_flux_w_m2: 0.7,
                    volumetric_heating_w_m3: Some(&heating),
                    ..RansHeatLoad::default()
                },
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(bits(&state), bits(&expected));
        assert_eq!(report, expected_report);
    }
    assert_eq!(work.numerical_bytes(), n * 9 * core::mem::size_of::<f64>());
    assert_eq!(
        pointers,
        [
            work.velocity.as_ptr(),
            work.kinetic.as_ptr(),
            work.dissipation.as_ptr(),
            work.viscosity.as_ptr(),
            work.production.as_ptr(),
            work.lower.as_ptr(),
            work.diagonal.as_ptr(),
            work.upper.as_ptr(),
            work.rhs.as_ptr()
        ]
    );
}
