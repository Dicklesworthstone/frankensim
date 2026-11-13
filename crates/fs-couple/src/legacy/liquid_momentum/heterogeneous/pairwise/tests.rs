use super::*;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn axes_z(angle: f64) -> [[f64; 3]; 3] {
    let (s, c) = angle.sin_cos();
    [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]
}
fn axes_x(angle: f64) -> [[f64; 3]; 3] {
    let (s, c) = angle.sin_cos();
    [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]]
}
const BODY_MASS: f64 = 2.4;
const INERTIA: [f64; 3] = [0.31, 0.53, 0.79];

fn cell(mass: f64, arm: [f64; 3], rates: [f64; 3], axes: [[f64; 3]; 3]) -> Cell {
    let principal = axes_z(0.37);
    let mut columns = [[0.0; 3]; 6];
    for i in 0..3 {
        columns[i][i] = 1.0 / BODY_MASS.sqrt();
        columns[i + 3] = cross(principal[i], arm).map(|v| v / INERTIA[i].sqrt());
    }
    Cell {
        mass,
        slip: [0.34, -0.27, 0.19],
        body_columns: columns,
        principal_rates_per_s: rates,
        principal_axes: axes,
    }
}
fn heterogeneous() -> (Vec<Cell>, Vec<[f64; 3]>) {
    let arms = vec![
        [0.43, -0.21, 0.17],
        [-0.22, 0.37, 0.11],
        [0.06, 0.12, -0.39],
    ];
    let mut cells = vec![
        cell(0.6, arms[0], [0.3, 1.4, 2.1], axes_x(0.41)),
        cell(0.9, arms[1], [2.7, 0.0, 0.8], axes_z(-0.23)),
        cell(0.4, arms[2], [0.9, 3.1, 0.2], axes_x(-0.52)),
    ];
    cells[1].slip = [-0.13, 0.22, 0.31];
    cells[2].slip = [0.27, 0.15, -0.24];
    (cells, arms)
}
fn run(cells: &[Cell], duration: f64, substeps: usize) -> Plan {
    solve(
        cells,
        duration,
        substeps,
        (2 * cells.len() - 1) * substeps,
        || false,
    )
    .unwrap()
}

/// Independent simultaneous ODE action in mass-whitened coordinates. This
/// never invokes the pair mobility, its eigensystem or a split exchange.
fn derivative(cells: &[Cell], state: &[f64]) -> Vec<f64> {
    let offset = cells.len() * 3;
    let mut result = vec![0.0; offset + 6];
    for (i, cell) in cells.iter().enumerate() {
        let slip: [f64; 3] = core::array::from_fn(|a| {
            state[3 * i + a] / cell.mass.sqrt()
                - (0..6)
                    .map(|k| cell.body_columns[k][a] * state[offset + k])
                    .sum::<f64>()
        });
        let mut resistance = [0.0; 3];
        for mode in 0..3 {
            let amplitude: f64 = (0..3).map(|a| cell.principal_axes[mode][a] * slip[a]).sum();
            for a in 0..3 {
                resistance[a] +=
                    cell.principal_axes[mode][a] * cell.principal_rates_per_s[mode] * amplitude;
            }
        }
        for a in 0..3 {
            result[3 * i + a] = -cell.mass.sqrt() * resistance[a];
        }
        for k in 0..6 {
            result[offset + k] += cell.mass
                * (0..3)
                    .map(|a| cell.body_columns[k][a] * resistance[a])
                    .sum::<f64>();
        }
    }
    result
}
fn initial(cells: &[Cell]) -> Vec<f64> {
    cells
        .iter()
        .flat_map(|cell| cell.slip.map(|v| v * cell.mass.sqrt()))
        .chain([0.0; 6])
        .collect()
}
fn state(cells: &[Cell], plan: &Plan) -> Vec<f64> {
    cells
        .iter()
        .zip(&plan.impulses)
        .flat_map(|(cell, j)| {
            core::array::from_fn::<_, 3, _>(|a| {
                cell.mass.sqrt() * cell.slip[a] + j[a] / cell.mass.sqrt()
            })
        })
        .chain(plan.body_increment)
        .collect()
}
fn reference(cells: &[Cell], duration: f64) -> Vec<f64> {
    // Full-system exponential series on short intervals; the reference cases
    // have bounded small rates, so no stiff series evaluation is attempted.
    let mut value = initial(cells);
    let h = duration / 32.0;
    for _ in 0..32 {
        let mut term = value.clone();
        let mut next = value.clone();
        for order in 1..=36 {
            term = derivative(cells, &term)
                .into_iter()
                .map(|v| v * h / f64::from(order))
                .collect();
            for (sum, add) in next.iter_mut().zip(&term) {
                *sum += add;
            }
        }
        value = next;
    }
    value
}
fn error(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}
fn energy(state: &[f64]) -> f64 {
    0.5 * state.iter().map(|x| x * x).sum::<f64>()
}

#[test]
fn one_cell_is_exact_with_finite_mass_recoil() {
    let cell = cell(0.6, [0.0; 3], [2.3; 3], axes_z(0.0));
    let duration = 0.73;
    let plan = run(core::slice::from_ref(&cell), duration, 1);
    let mobility = 1.0 + cell.mass / BODY_MASS;
    let decay = (-2.3 * mobility * duration).exp();
    for a in 0..3 {
        let impulse = -cell.mass * cell.slip[a] * (1.0 - decay) / mobility;
        near(plan.impulses[0][a], impulse, 2e-14);
        near(plan.body_increment[a], -impulse / BODY_MASS.sqrt(), 2e-14);
    }
    let loss = 0.5 * cell.mass / mobility
        * cell.slip.iter().map(|v| v * v).sum::<f64>()
        * (1.0 - decay * decay);
    near(plan.dissipation_j, loss, 2e-14);
    assert_eq!(plan.pair_exchanges, 1);
}

#[test]
fn noncommuting_anisotropic_pair_matches_simultaneous_exponential() {
    let cell = cell(0.8, [0.43, -0.21, 0.17], [0.4, 1.7, 3.2], axes_x(0.49));
    let mut gamma = [[0.0; 3]; 3];
    let mut mobility = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            gamma[i][j] = (0..3)
                .map(|k| {
                    cell.principal_axes[k][i]
                        * cell.principal_rates_per_s[k]
                        * cell.principal_axes[k][j]
                })
                .sum();
            mobility[i][j] = if i == j { 1.0 } else { 0.0 }
                + cell.mass
                    * (0..6)
                        .map(|k| cell.body_columns[k][i] * cell.body_columns[k][j])
                        .sum::<f64>();
        }
    }
    let commutator: f64 = (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .map(|(i, j)| {
            (0..3)
                .map(|k| mobility[i][k] * gamma[k][j] - gamma[i][k] * mobility[k][j])
                .sum::<f64>()
                .abs()
        })
        .sum();
    assert!(commutator > 0.01);
    let cells = [cell];
    let plan = run(&cells, 0.31, 1);
    assert!(error(&state(&cells, &plan), &reference(&cells, 0.31)) < 3e-13);
    near(
        plan.dissipation_j,
        energy(&initial(&cells)) - energy(&state(&cells, &plan)),
        3e-14,
    );
}

#[test]
fn heterogeneous_palindromic_composition_is_second_order() {
    let (cells, _) = heterogeneous();
    let expected = reference(&cells, 0.63);
    let errors: Vec<_> = [2, 4, 8, 16]
        .into_iter()
        .map(|steps| {
            let plan = run(&cells, 0.63, steps);
            assert_eq!(plan.pair_exchanges, (2 * cells.len() - 1) * steps);
            error(&state(&cells, &plan), &expected)
        })
        .collect();
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.28 * pair[0], "refinement errors: {errors:?}");
    }
    assert!(errors[3] < 2e-5, "refinement errors: {errors:?}");
}

#[test]
fn full_inertia_recoil_closes_linear_angular_and_physical_energy() {
    let (cells, arms) = heterogeneous();
    let plan = run(&cells, 0.63, 8);
    let principal = axes_z(0.37);
    let mut angular = [0.0; 3];
    for a in 0..3 {
        let liquid = plan.impulses.iter().map(|j| j[a]).sum::<f64>();
        near(
            liquid + BODY_MASS.sqrt() * plan.body_increment[a],
            0.0,
            3e-14,
        );
        angular[a] = (0..3)
            .map(|k| principal[k][a] * INERTIA[k].sqrt() * plan.body_increment[k + 3])
            .sum();
    }
    for (arm, impulse) in arms.iter().zip(&plan.impulses) {
        let moment = cross(*arm, *impulse);
        for a in 0..3 {
            angular[a] += moment[a];
        }
    }
    for component in angular {
        near(component, 0.0, 3e-14);
    }
    assert!(plan.dissipation_j > 0.0);
    near(
        energy(&state(&cells, &plan)) - energy(&initial(&cells)) + plan.dissipation_j,
        0.0,
        1e-13,
    );
    assert_eq!(run(&cells, 0.63, 8), plan);
}

#[test]
fn stiff_null_directions_stay_exact_with_off_center_mobility() {
    let cell = cell(0.8, [0.43, -0.21, 0.17], [0.0, 1e300, 0.0], axes_z(0.0));
    let expected = -cell.mass * cell.slip[1]
        / (1.0 + cell.mass * cell.body_columns.iter().map(|k| k[1] * k[1]).sum::<f64>());
    let cells = [cell];
    let plan = run(&cells, 1e10, 1);
    assert_eq!(plan.impulses[0][0], 0.0);
    assert_eq!(plan.impulses[0][2], 0.0);
    near(plan.impulses[0][1], expected, 2e-14);
    near(
        energy(&state(&cells, &plan)) - energy(&initial(&cells)) + plan.dissipation_j,
        0.0,
        3e-14,
    );
}

#[test]
fn weak_direction_survives_overflowing_common_exposure() {
    let cell = cell(0.6, [0.0; 3], [1e300, 1e-10, 0.0], axes_z(0.0));
    let mobility = 1.0 + cell.mass / BODY_MASS;
    let plan = run(core::slice::from_ref(&cell), 1e10, 1);
    near(
        plan.impulses[0][0],
        -cell.mass * cell.slip[0] / mobility,
        2e-14,
    );
    near(
        plan.impulses[0][1],
        -cell.mass * cell.slip[1] / mobility * -(-mobility).exp_m1(),
        2e-11,
    );
    assert_eq!(plan.impulses[0][2], 0.0);
}

#[test]
fn zero_duration_and_zero_resistance_preserve_all_increment_bits() {
    let mut cell = cell(0.8, [0.43, -0.21, 0.17], [0.4, 1.7, 3.2], axes_x(0.49));
    for duration in [0.0, 1.0] {
        if duration > 0.0 {
            cell.principal_rates_per_s = [0.0; 3];
        }
        let plan = run(core::slice::from_ref(&cell), duration, 3);
        assert_eq!(plan.pair_exchanges, 0);
        for value in plan.impulses[0]
            .into_iter()
            .chain(plan.body_increment)
            .chain([plan.dissipation_j])
        {
            assert_eq!(value.to_bits(), 0.0_f64.to_bits());
        }
    }
}

#[test]
fn invalid_or_unbounded_sources_and_mid_composition_cancellation_refuse() {
    let (cells, _) = heterogeneous();
    assert!(work_bound(0, 1, 10).is_err());
    assert!(work_bound(3, 0, 10).is_err());
    assert!(work_bound(usize::MAX, 1, usize::MAX).is_err());
    assert!(solve(&cells, 0.1, 2, 9, || false).is_err());
    assert!(solve(&cells, -0.1, 1, 5, || false).is_err());
    let mut invalid_cells = cells.clone();
    invalid_cells[0].principal_rates_per_s[1] = -1.0;
    assert!(solve(&invalid_cells, 0.1, 1, 5, || false).is_err());
    invalid_cells[0].principal_rates_per_s[1] = f64::NAN;
    assert!(solve(&invalid_cells, 0.1, 1, 5, || false).is_err());
    let mut polls = 0;
    assert_eq!(
        solve(&cells, 0.1, 4, 20, || {
            polls += 1;
            polls == 7
        }),
        Err(Error::Cancelled)
    );
}

#[test]
fn positive_exposure_lost_to_underflow_refuses() {
    let mut cell = cell(1e200, [0.0; 3], [1e-200; 3], axes_z(0.0));
    cell.body_columns = [[0.0; 3]; 6];
    for i in 0..3 {
        cell.body_columns[i][i] = 1e-100;
    }
    // Exposure is below binary64, although mass * rate * duration is not.
    assert!(matches!(
        solve(&[cell], 2e-130, 1, 1, || false),
        Err(Error::Invalid(reason)) if reason.contains("exposure")
    ));
}
