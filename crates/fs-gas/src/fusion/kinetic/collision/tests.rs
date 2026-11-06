use super::*;
use fs_qty::Amount;

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn two_bath_operator() -> (AlphaEnergyGrid, AlphaCollisionOperator) {
    let grid =
        AlphaEnergyGrid::new(vec![Energy::new(1e-16), Energy::new(3e-16)], || false).unwrap();
    let mut up = [[0.0; MAX_EDGES]; 2];
    let mut down = [[0.0; MAX_EDGES]; 2];
    up[0][0] = 0.6;
    up[1][0] = 1.4;
    down[0][0] = 1.2;
    down[1][0] = 1.8;
    let operator = AlphaCollisionOperator::from_rates(&grid, &up, &down, || false).unwrap();
    (grid, operator)
}

#[test]
fn g1_two_cell_relaxation_and_each_bath_reward_match_the_exact_solution() {
    let (grid, operator) = two_bath_operator();
    let mut initial = [0.0; MAX_ALPHA_BINS];
    initial[0] = 4e-6;
    initial[1] = 1e-6;
    let total = 5e-6;
    let equilibrium_high = total * 2.0 / 5.0;
    for duration in [1e-12, 0.37, 100.0] {
        let result = operator
            .advance(&initial, Time::new(duration), 100_000, || false)
            .unwrap();
        // Exact N_high'=2*(N_total-N_high)-3*N_high, and exact
        // integral of each bath's (down*N_high-up*N_low)*DeltaE.
        let high = equilibrium_high + (initial[1] - equilibrium_high) * (-5.0 * duration).exp();
        let high_integral = equilibrium_high * duration
            + (initial[1] - equilibrium_high) * (-(-5.0 * duration).exp_m1()) / 5.0;
        close(result.amounts[1], high, 2e-12, 1e-20);
        close(result.amounts[0], total - high, 2e-12, 1e-20);
        for (bath, (up, down)) in [(0.6, 1.2), (1.4, 1.8)].into_iter().enumerate() {
            let expected =
                (2e-16 * AVOGADRO) * ((up + down) * high_integral - up * total * duration);
            close(result.deposited_j[bath], expected, 3e-12, 2e-10);
        }
        close(result.amounts.iter().sum(), total, 2e-12, 1e-20);
        assert!(
            result
                .amounts
                .iter()
                .all(|value| value.is_finite() && *value >= 0.0)
        );
        assert!(result.numerical_tail_bound < 1e-15);
        assert!(result.applications > 0);
        assert!(result.energy_residual_j.abs() < 2e-8);
        close(
            energy(&result.amounts, &grid).unwrap() + result.deposited_j.iter().sum::<f64>(),
            energy(&initial, &grid).unwrap(),
            3e-12,
            2e-8,
        );
    }
}

#[test]
fn g1_maxwellian_detailed_balance_is_stationary_and_refines_to_three_halves_kt() {
    let temperature = Temperature::new(1e7);
    let kt = BOLTZMANN_J_PER_K * temperature.value();
    let grid =
        AlphaEnergyGrid::uniform_speed(Energy::new(0.01 * kt), Energy::new(30.0 * kt), 64, || {
            false
        })
        .unwrap();
    // D_b=nu_b*kT*E has the Maxwellian energy equilibrium. This direct
    // diffusion fixture tests the energy discretization independently of the
    // Coulomb source's coefficient implementation.
    let operator = AlphaCollisionOperator::from_diffusion(
        &grid,
        [temperature; 2],
        |e| Ok([0.4 * kt * e.value(), 0.6 * kt * e.value()]),
        || false,
    )
    .unwrap();
    let weights: Vec<_> = grid
        .energies()
        .iter()
        .zip(grid.widths_j())
        .map(|(e, width)| width * e.value().sqrt() * (-e.value() / kt).exp())
        .collect();
    let normalization: f64 = weights.iter().sum();
    let mut equilibrium = [0.0; MAX_ALPHA_BINS];
    for (n, weight) in equilibrium.iter_mut().zip(weights) {
        *n = 7e-6 * weight / normalization;
    }
    for edge in 0..grid.len() - 1 {
        close(
            operator.upward[edge] * equilibrium[edge],
            operator.downward[edge + 1] * equilibrium[edge + 1],
            5e-13,
            1e-30,
        );
    }
    let relaxed = operator
        .advance(&equilibrium, Time::new(0.2), 100_000, || false)
        .unwrap();
    for (actual, expected) in relaxed.amounts.iter().zip(equilibrium) {
        close(*actual, expected, 2e-12, 1e-25);
    }
    for heat in relaxed.deposited_j {
        assert!(heat.abs() < 1e-8);
    }

    // Midpoints in speed approach zero with refinement. Merely refining at a
    // fixed positive first energy would retain its first-cell quadrature error.
    let mut previous_error = f64::INFINITY;
    for count in [16, 32, 64, 128] {
        let spacing = (30.0_f64 * kt).sqrt() / count as f64;
        let grid = AlphaEnergyGrid::new(
            (0..count)
                .map(|i| Energy::new(((i as f64 + 0.5) * spacing).powi(2)))
                .collect(),
            || false,
        )
        .unwrap();
        let mut population = 0.0;
        let mut energy = 0.0;
        for (e, width) in grid.energies().iter().zip(grid.widths_j()) {
            let weight = width * e.value().sqrt() * (-e.value() / kt).exp();
            population += weight;
            energy += e.value() * weight;
        }
        let error = (energy / population / kt - 1.5).abs();
        assert!(
            error < previous_error,
            "count={count}, error={error}, previous={previous_error}"
        );
        previous_error = error;
    }
    assert!(
        previous_error < 3e-5,
        "final equilibrium mean error={previous_error}"
    );
}

#[test]
fn g1_cold_kinetic_alphas_gain_energy_without_a_terminal_handoff() {
    let temperature = Temperature::new(1e7);
    let kt = BOLTZMANN_J_PER_K * temperature.value();
    let grid =
        AlphaEnergyGrid::uniform_speed(Energy::new(0.01 * kt), Energy::new(20.0 * kt), 48, || {
            false
        })
        .unwrap();
    let operator = AlphaCollisionOperator::from_diffusion(
        &grid,
        [temperature; 2],
        |e| Ok([0.0, kt * e.value()]),
        || false,
    )
    .unwrap();
    let mut initial = [0.0; MAX_ALPHA_BINS];
    grid.inject(&mut initial, Amount::new(2e-6), grid.energies()[0])
        .unwrap();
    let result = operator
        .advance(&initial, Time::new(0.5), 100_000, || false)
        .unwrap();
    close(result.amounts.iter().sum(), 2e-6, 2e-12, 1e-20);
    assert!(result.amounts[1..grid.len()].iter().any(|&n| n > 0.0));
    assert!(energy(&result.amounts, &grid).unwrap() > energy(&initial, &grid).unwrap());
    assert!(result.deposited_j[1] < 0.0);
    assert_eq!(result.deposited_j[0], 0.0);
    assert!(result.energy_residual_j.abs() < 1e-8);
}

#[test]
fn g4_collision_budgets_domains_and_late_cancellation_publish_no_partial_state() {
    let (_, operator) = two_bath_operator();
    let mut amounts = [0.0; MAX_ALPHA_BINS];
    amounts[0] = 1e-6;
    let original = amounts;
    let noop = operator
        .advance(&amounts, Time::new(0.0), 0, || false)
        .unwrap();
    assert_eq!(noop.amounts, amounts);
    assert_eq!(noop.deposited_j, [0.0; 2]);
    assert_eq!(noop.applications, 0);
    assert_eq!(
        operator
            .advance(&amounts, Time::new(1.0), 1, || false)
            .unwrap_err(),
        DtFusionError::StepBudget
    );
    assert!(matches!(
        operator.advance(&amounts, Time::new(-1.0), 1000, || false),
        Err(DtFusionError::Invalid(_))
    ));
    assert_eq!(
        operator
            .advance(&amounts, Time::new(1.0), 1000, || true)
            .unwrap_err(),
        DtFusionError::Cancelled
    );
    let mut polls = 0;
    assert_eq!(
        operator
            .advance(&amounts, Time::new(100.0), 100_000, || {
                polls += 1;
                polls > 30
            })
            .unwrap_err(),
        DtFusionError::Cancelled
    );
    assert_eq!(amounts, original);
    amounts[100] = 1e-6;
    assert!(matches!(
        operator.advance(&amounts, Time::new(1.0), 1000, || false),
        Err(DtFusionError::Invalid(_))
    ));
}
