use super::*;
use fs_thermochem::coulomb::{CoulombLimitsV1, CoulombScreeningV1};

fn plasma(ne: f64, neutral: f64, te: f64, th: f64, formation: f64) -> FrozenIonizedGasStateV1 {
    let source = state(te, th, 1.0, formation);
    FrozenIonizedGasStateV1::new(
        source.thermodynamics(),
        &[
            Amount::new(neutral / AVOGADRO_PER_MOL),
            Amount::new(ne / AVOGADRO_PER_MOL),
            Amount::new(ne / AVOGADRO_PER_MOL),
        ],
        Volume::new(1.0),
        Temperature::new(th),
        Temperature::new(te),
        1e-12,
        || false,
    )
    .unwrap()
}
fn coulomb(source: &FrozenIonizedGasStateV1) -> CoulombModelV1 {
    CoulombModelV1::new(
        source,
        CoulombScreeningV1::StaticAllChargedDebye,
        CoulombLimitsV1 {
            minimum_coulomb_logarithm: 4.0,
            minimum_particles_per_debye_sphere: 10.0,
            maximum_coulomb_coupling: 0.1,
            maximum_thermal_phase_space_density: 0.1,
            maximum_thermal_energy_rest_energy_ratio: 0.001,
            maximum_quantum_to_classical_ratio: 0.5,
            maximum_relative_drift_thermal_ratio: 0.1,
        },
        || false,
    )
    .unwrap()
}
fn plasma_controls() -> ElectronHeavyEnergySolve {
    let mut c = controls();
    c.maximum_charged_mole_fraction = 1.0;
    c.energy_tolerance_j = 1e-9;
    c.initial_step_s = 1e-9;
    c.maximum_step_s = 1e-6;
    c.minimum_step_s = 1e-22;
    c
}

#[test]
fn coulomb_only_thermal_transient_converges_to_independent_landau_moment_integral() {
    let initial = plasma(1e20, 0.0, 20000.0, 10000.0, 0.0);
    let model = coulomb(&initial);
    // Independent coherent-SI calculation of the declared static-Debye,
    // classical-cutoff Maxwellian moment. QUADPACK inversion of
    // t=integral[Cv dTe / (G(Te,30000-Te)*(2Te-30000))]
    // and DOP853 at rtol1e-12 agree within5e-9K. Neither uses Rust output.
    // See Fitzpatrick Eqs3.179/3.183:
    // https://farside.ph.utexas.edu/teaching/plasma/lectures/node42.html
    let reference_te = 16_755.611_400_759_084;
    let mut errors = [0.0; 3];
    // The endpoint-consistent backward-Euler method is first order; its
    // configured local estimate is not a global temperature-error bound.
    // Resolve the tightest run sufficiently for the independent 0.5 K target.
    for (i, relative) in [1e-3, 1e-5, 1e-8].into_iter().enumerate() {
        let mut c = plasma_controls();
        c.relative_tolerance = relative;
        let result = advance_frozen_energy_with_coulomb(
            &initial,
            &[],
            &model,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(5e-7),
            c,
            || false,
        )
        .unwrap();
        let final_state = result.state();
        let te = final_state.electron_temperature().value();
        let th = final_state.heavy_temperature().value();
        errors[i] = (te - reference_te).abs();
        assert!(te < 20000.0 && th > 10000.0 && te > th);
        assert!((te + th - 30000.0).abs() < 1e-7);
        assert_eq!(final_state.amounts(), initial.amounts());
        assert_eq!(result.coulomb_model(), Some(&model));
        let b = result.balance();
        assert!(b.electron_energy_residual_j.abs() <= c.energy_tolerance_j);
        assert!(b.heavy_energy_residual_j.abs() <= c.energy_tolerance_j);
        assert_eq!(b.electron_to_heavy_neutral_energy_j, 0.0);
        assert_eq!(
            b.electron_to_heavy_energy_j,
            b.electron_to_heavy_coulomb_energy_j
        );
        assert!(b.electron_to_heavy_coulomb_energy_j > 0.0);
        let cv = 1.5 * BOLTZMANN_J_PER_K * 1e20;
        let expected_entropy = cv * ((te / 20000.0).ln() + (th / 10000.0).ln());
        assert!(expected_entropy > 0.0);
        let entropy = result.entropy_balance();
        assert!((entropy.total_entropy_change_j_k() - expected_entropy).abs() < 2e-17);
        assert!(entropy.collision_production_j_k > 0.0);
        assert_eq!(entropy.external_power_contribution_j_k, 0.0);
        assert!(entropy.numerical_excess_j_k >= -entropy.tolerance_j_k);
        assert!(entropy.residual_j_k.abs() <= entropy.tolerance_j_k);
        assert!(b.maximum_normalized_local_error <= 1.0);
    }
    eprintln!("Coulomb transient temperature errors [K]: {errors:?}");
    assert!(
        errors[1] < 0.2 * errors[0] && errors[2] < 0.2 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 0.5, "{errors:?}");
}

#[test]
fn neutral_and_coulomb_exchange_share_one_conservative_energy_step() {
    let initial = plasma(1e18, 1e22, 20000.0, 10000.0, 0.0);
    let model = coulomb(&initial);
    let sources = [collision(1e-23, 500.0)];
    let c = plasma_controls();
    let combined = advance_frozen_energy_with_coulomb(
        &initial,
        &sources,
        &model,
        Power::new(0.0),
        Power::new(0.0),
        Time::new(2e-5),
        c,
        || false,
    )
    .unwrap();
    let mut neutral_controls = c;
    neutral_controls.maximum_charged_mole_fraction = 0.01;
    let neutral = advance_frozen_energy(
        &initial,
        &sources,
        Power::new(0.0),
        Power::new(0.0),
        Time::new(2e-5),
        neutral_controls,
        || false,
    )
    .unwrap();
    let b = combined.balance();
    assert!(b.electron_to_heavy_coulomb_energy_j > 0.0);
    assert!(b.electron_to_heavy_neutral_energy_j > 0.0);
    assert!(
        (b.electron_to_heavy_coulomb_energy_j + b.electron_to_heavy_neutral_energy_j
            - b.electron_to_heavy_energy_j)
            .abs()
            < 1e-15
    );
    assert!(
        combined.state().electron_temperature().value()
            < neutral.state().electron_temperature().value()
    );
    assert!(
        (combined.state().internal_energy().value() - initial.internal_energy().value()).abs()
            < 1e-8
    );
    assert_eq!(combined.state().amounts(), initial.amounts());
}

#[test]
fn coulomb_exchange_reverses_and_closes_held_work_with_large_formation_reference() {
    for shift in [0.0, 1e12] {
        let initial = plasma(1e20, 0.0, 10000.0, 20000.0, shift);
        let model = coulomb(&initial);
        let c = plasma_controls();
        let duration = 3e-7;
        let pe = 1e6;
        let ph = -2e5;
        let result = advance_frozen_energy_with_coulomb(
            &initial,
            &[],
            &model,
            Power::new(pe),
            Power::new(ph),
            Time::new(duration),
            c,
            || false,
        )
        .unwrap();
        let cv = 1.5 * BOLTZMANN_J_PER_K * 1e20;
        let de = cv * (result.state().electron_temperature().value() - 10000.0);
        let dh = cv * (result.state().heavy_temperature().value() - 20000.0);
        let b = result.balance();
        assert!(b.electron_to_heavy_coulomb_energy_j < 0.0);
        assert!((de - pe * duration + b.electron_to_heavy_coulomb_energy_j).abs() < 1e-8);
        assert!((dh - ph * duration - b.electron_to_heavy_coulomb_energy_j).abs() < 1e-8);
        assert!(b.electron_energy_residual_j.abs() <= c.energy_tolerance_j);
        assert!(b.heavy_energy_residual_j.abs() <= c.energy_tolerance_j);
    }
}

#[test]
fn accepted_mixed_steps_use_endpoint_rates_for_channels_and_entropy() {
    let initial = plasma(1e18, 1e22, 20000.0, 10000.0, 0.0);
    let model = coulomb(&initial);
    let coulomb_at = |temperatures: [f64; 2]| {
        model
            .evaluate(
                Temperature::new(temperatures[0]),
                Temperature::new(temperatures[1]),
                || false,
            )
            .unwrap()
            .electron_heavy_thermal_conductance_w_m3_k()
            * initial.volume().value()
    };
    let start = initial_point(&initial);
    let initial_coulomb = coulomb_at(start.temperatures);
    // Make the two physically different channels comparable at the start,
    // then allow their relative strength to evolve with the temperatures.
    let sigma = 1e-23 * initial_coulomb / constant_sigma_conductance(&initial, 1e-23, 20000.0);
    let sources = [collision(sigma, 500.0)];
    let dt = 0.5 * start.capacities[1] / initial_coulomb;
    let mut c = plasma_controls();
    c.relative_tolerance = 0.5;
    c.initial_step_s = dt;
    c.maximum_step_s = dt;
    c.maximum_steps = 1;
    let problem =
        Problem::new(&initial, &sources, Some(&model), [0.0; 2], c, &mut || false).unwrap();
    let initial_coupling = problem.coupling(&start, &mut || false).unwrap();
    let (half, half_coupling) = problem
        .step(&start, 0.5 * dt, initial_coupling, &mut || false)
        .unwrap();
    let (two, _) = problem
        .step(&half, dt, half_coupling, &mut || false)
        .unwrap();

    let mut neutral_energy = 0.0;
    let mut coulomb_energy = 0.0;
    let mut collision_entropy = 0.0;
    for point in [half, two] {
        let [th, te] = point.temperatures;
        let neutral = constant_sigma_conductance(&initial, sigma, te);
        let charged = coulomb_at(point.temperatures);
        neutral_energy += 0.5 * dt * neutral * (te - th);
        coulomb_energy += 0.5 * dt * charged * (te - th);
        collision_entropy += 0.5 * dt * (neutral + charged) * (te - th).powi(2) / (te * th);
    }
    let end_coulomb = coulomb_at(two.temperatures);
    let end_neutral = constant_sigma_conductance(&initial, sigma, two.temperatures[1]);
    assert!(end_coulomb / (end_coulomb + end_neutral) > 0.55);

    let result = advance_frozen_energy_with_coulomb(
        &initial,
        &sources,
        &model,
        Power::new(0.0),
        Power::new(0.0),
        Time::new(dt),
        c,
        || false,
    )
    .unwrap();
    let balance = result.balance();
    assert_eq!(balance.accepted_steps, 1);
    assert_eq!(balance.rejected_steps, 0);
    assert_eq!(balance.electron_to_heavy_energy_j, two.exchange);
    assert_eq!(
        balance.electron_to_heavy_coulomb_energy_j,
        two.coulomb_exchange
    );
    assert!(
        (balance.electron_to_heavy_neutral_energy_j - neutral_energy).abs()
            < 4.0 * c.energy_tolerance_j
    );
    assert!(
        (balance.electron_to_heavy_coulomb_energy_j - coulomb_energy).abs()
            < 4.0 * c.energy_tolerance_j
    );
    let entropy = result.entropy_balance();
    assert!(
        (entropy.collision_production_j_k - collision_entropy).abs() < 1e-9 * collision_entropy
    );
    assert!(entropy.numerical_excess_j_k >= -entropy.tolerance_j_k);
    assert!(entropy.residual_j_k.abs() <= entropy.tolerance_j_k);
}

#[test]
fn sourced_coulomb_exchange_brackets_a_nonmonotone_constitutive_residual() {
    let initial = plasma(1e18, 0.0, 30000.0, 1500.0, 0.0);
    let model = coulomb(&initial);
    let start = initial_point(&initial);
    let conductance = |th, te| {
        model
            .evaluate(Temperature::new(th), Temperature::new(te), || false)
            .unwrap()
            .electron_heavy_thermal_conductance_w_m3_k()
            * initial.volume().value()
    };
    let dt = 10.0 * start.capacities[1] / conductance(1500.0, 30000.0);
    let residual = |q: f64| {
        let th = 1500.0 + q / start.capacities[0];
        let te = 30000.0 - q / start.capacities[1];
        q - dt * conductance(th, te) * (te - th)
    };
    // This is an actual admitted source law, not a manufactured solver
    // callback. Static Debye screening and the thermal moment make the
    // endpoint residual initially DECREASE as energy flows to the heavy pool.
    // A derivative that assumes frozen G is not the nonlinear derivative.
    assert!(residual(10.0 * start.capacities[1]) < residual(0.0));
    let c = plasma_controls();
    let problem = Problem::new(&initial, &[], Some(&model), [0.0; 2], c, &mut || false).unwrap();
    let coupling = problem.coupling(&start, &mut || false).unwrap();
    let (end, endpoint_coupling) = problem.step(&start, dt, coupling, &mut || false).unwrap();
    assert!(end.temperatures[1] > end.temperatures[0]);
    assert!(end.temperatures[0] > 1500.0);
    assert!(end.temperatures[1] < 30000.0);
    assert!(residual(end.exchange).abs() < 4.0 * c.energy_tolerance_j);
    assert_eq!(end.coulomb_exchange, end.exchange);
    let entropy = problem
        .entropy_step(&start, &end, endpoint_coupling, &mut || false)
        .unwrap();
    assert!(entropy.collision > 0.0);
    assert!(entropy.numerical >= -entropy.tolerance);
}

#[test]
fn coulomb_population_binding_domain_and_cancellation_refuse_atomically() {
    let initial = plasma(1e20, 0.0, 20000.0, 10000.0, 0.0);
    let model = coulomb(&initial);
    let changed = plasma(2e20, 0.0, 20000.0, 10000.0, 0.0);
    let before = changed.clone();
    assert!(matches!(
        advance_frozen_energy_with_coulomb(
            &changed,
            &[],
            &model,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-6),
            plasma_controls(),
            || false
        ),
        Err(ElectronHeavyEnergyError::Invalid(_))
    ));
    assert_eq!(changed, before);
    assert_eq!(
        advance_frozen_energy_with_coulomb(
            &initial,
            &[],
            &model,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-6),
            plasma_controls(),
            || true
        ),
        Err(ElectronHeavyEnergyError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        advance_frozen_energy_with_coulomb(
            &initial,
            &[],
            &model,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-6),
            plasma_controls(),
            || {
                polls += 1;
                polls > 100
            }
        ),
        Err(ElectronHeavyEnergyError::Cancelled)
    );
    // The old neutral-only closure still refuses a fully charged fraction.
    assert!(matches!(
        advance_frozen_energy(
            &initial,
            &[],
            Power::new(0.0),
            Power::new(0.0),
            Time::new(0.0),
            plasma_controls(),
            || false
        ),
        Err(ElectronHeavyEnergyError::Invalid(_))
    ));
    let zero = advance_frozen_energy_with_coulomb(
        &initial,
        &[],
        &model,
        Power::new(0.0),
        Power::new(0.0),
        Time::new(0.0),
        plasma_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(zero.state(), &initial);
    assert_eq!(zero.balance().electron_to_heavy_coulomb_energy_j, 0.0);
}
