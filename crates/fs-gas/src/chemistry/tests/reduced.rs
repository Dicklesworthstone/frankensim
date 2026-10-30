use super::*;
use crate::chemistry::solve::NewtonWorkspace;

#[test]
fn reduced_newton_matches_full_extents_with_dependent_stoichiometry_and_heat() {
    let ns = 3;
    let nr = 7;
    // Three cyclic reaction directions, their duplicates and one doubled
    // direction. N has rank two; neither factorization may require full rank.
    let nu = [
        -1.0, 0.0, 1.0, -1.0, 0.0, 1.0, -2.0, 1.0, -1.0, 0.0, 1.0, -1.0, 0.0, 2.0, 0.0, 1.0, -1.0,
        0.0, 1.0, -1.0, 0.0,
    ];
    let jac = [
        3.0, -1.0, 0.0, 0.0, 7.0, -2.0, -5.0, 0.0, 11.0, 1.5, -0.5, 0.0, 0.0, 4.0, -1.0, -2.0, 0.0,
        6.0, 8.0, -3.0, 0.0,
    ];
    let temperature_rates = [0.003, 0.014, -0.007, 0.002, 0.008, -0.004, 0.006];
    let u = [-6000.0, -12000.0, -23000.0];
    let temperature_extents: Vec<f64> = (0..nr)
        .map(|j| -(0..ns).map(|i| u[i] * nu[i * nr + j]).sum::<f64>() / 900.0)
        .collect();
    let residual = [1.0, -0.4, 0.3, -0.2, 0.8, -1.1, 0.17];
    let mut reduced = NewtonWorkspace::new(ns, nr, &nu).unwrap();
    let mut full = NewtonWorkspace::full_extent(nr).unwrap();
    assert_eq!(reduced.dimension(), 4);
    assert_eq!(full.dimension(), 7);
    for step in [1e-8, 0.2, 1000.0] {
        let mut actual = vec![0.0; nr];
        let mut expected = vec![0.0; nr];
        for (solver, output) in [(&mut reduced, &mut actual), (&mut full, &mut expected)] {
            solver
                .solve(
                    &nu,
                    &jac,
                    &temperature_rates,
                    &temperature_extents,
                    850.0,
                    step,
                    &residual,
                    output,
                    &mut || false,
                )
                .unwrap();
        }
        for j in 0..nr {
            assert!((actual[j] - expected[j]).abs() < 2e-9 * (1.0 + expected[j].abs()));
            // Independently assemble and apply the full extent matrix. This
            // check does not use the Schur projection or its residual helper.
            let mut applied = 0.0;
            let mut magnitude = residual[j].abs();
            for k in 0..nr {
                let derivative = (0..ns)
                    .map(|i| jac[j * ns + i] * nu[i * nr + k])
                    .sum::<f64>()
                    + temperature_rates[j] * temperature_extents[k];
                let a = if j == k { 1.0 } else { 0.0 } - step * derivative;
                applied += a * actual[k];
                magnitude += (a * actual[k]).abs();
            }
            assert!((applied + residual[j]).abs() <= 3e-12 * (1.0 + magnitude));
        }
    }
    assert_eq!(
        NewtonWorkspace::new(10, 29, &[0.0; 290])
            .unwrap()
            .dimension(),
        11
    );
    assert_eq!(
        NewtonWorkspace::new(3, 3, &[0.0; 9]).unwrap().dimension(),
        3
    );
    assert_eq!(
        NewtonWorkspace::new(3, 4, &[0.0; 12]).unwrap().dimension(),
        4
    );
}

#[test]
fn stiff_schur_recovery_closes_the_original_extent_equation() {
    let nr = 8;
    let rate = 1e8;
    let step = 1.0;
    let nu: Vec<_> = [-1.0, 1.0]
        .into_iter()
        .flat_map(|v| std::iter::repeat_n(v, nr))
        .collect();
    let jac: Vec<_> = (0..nr).flat_map(|_| [rate, -rate]).collect();
    let residual = vec![-step * rate * 100.0; nr];
    let mut delta = vec![0.0; nr];
    let mut solver = NewtonWorkspace::new(2, nr, &nu).unwrap();
    solver
        .solve(
            &nu,
            &jac,
            &vec![0.0; nr],
            &vec![0.0; nr],
            500.0,
            step,
            &residual,
            &mut delta,
            &mut || false,
        )
        .unwrap();
    // The eight identical A <-> B rates have one nonzero eigenmode. Each
    // extent is independently known, including the identity in the stiff J.
    let expected = step * rate * 100.0 / (1.0 + 2.0 * nr as f64 * step * rate);
    for correction in delta {
        assert!(
            (correction - expected).abs() < 2e-11,
            "{correction} vs {expected}"
        );
    }
}

fn duplicate_isomers(q: f64, reference: f64, rate: f64, count: usize) -> IdealGasKineticsV1 {
    mechanism(
        vec![
            source("A", 0.01, 4.0, reference, 0.0),
            source("B", 0.01, 4.0, reference - q, 0.0),
        ],
        vec![1, 1],
        (0..count)
            .map(|j| {
                reaction(
                    &format!("A-B-{j}"),
                    &[("A", 1)],
                    &[("B", 1)],
                    rate / count as f64,
                )
            })
            .collect(),
    )
}

#[test]
fn reduced_stiff_reactor_preserves_elements_energy_and_analytic_relaxation() {
    let mechanism = duplicate_isomers(0.0, 0.0, 1e8, 8);
    let result = run(&mechanism, 500.0, 0.1, controls(1e-5));
    for amount in result.amounts() {
        assert!((amount.value() - 0.5).abs() < 2e-6);
        assert!(amount.value() >= 0.0);
    }
    assert!(
        result.balance().accepted_steps < 100,
        "{:?}",
        result.balance()
    );
    assert!(result.balance().mass_residual_kg.abs() < 1e-15);
    assert!(result.balance().maximum_relative_element_residual < 1e-14);
    assert!(result.balance().energy_residual_j.abs() <= controls(1e-5).energy_tolerance_j);
    assert!(result.balance().entropy_change_j_k > 0.0);
}

#[test]
fn reduced_exothermic_reactor_matches_lumped_rates_and_elemental_reference_shift() {
    let mut control = controls(2e-7);
    // At a 1e7 J/mol reference offset, this absolute energy gate still resolves
    // temperature far more tightly than the independent time-error tolerance.
    control.energy_tolerance_j = 1e-6;
    let lumped = run(
        &duplicate_isomers(12000.0, 0.0, 20.0, 1),
        300.0,
        0.15,
        control,
    );
    let duplicate = run(
        &duplicate_isomers(12000.0, 0.0, 20.0, 8),
        300.0,
        0.15,
        control,
    );
    let shifted = run(
        &duplicate_isomers(12000.0, 1e7, 20.0, 8),
        300.0,
        0.15,
        control,
    );
    for result in [&duplicate, &shifted] {
        for (a, b) in result.amounts().iter().zip(lumped.amounts()) {
            assert!((a.value() - b.value()).abs() < 2e-6);
        }
        let product = result.amounts()[1].value();
        let expected_temperature = 300.0 + 12000.0 * product / (3.0 * R);
        assert!((result.temperature().value() - expected_temperature).abs() < 2e-6);
        assert!(result.balance().mass_residual_kg.abs() < 1e-15);
        assert!(result.balance().maximum_relative_element_residual < 1e-14);
        assert!(result.balance().energy_residual_j.abs() <= control.energy_tolerance_j);
        assert!(result.balance().entropy_change_j_k > 0.0);
    }
}

#[test]
fn reduced_workspace_and_reactor_replay_after_cancellation() {
    let nu = [-1.0; 4].into_iter().chain([1.0; 4]).collect::<Vec<_>>();
    let mut workspace = NewtonWorkspace::new(2, 4, &nu).unwrap();
    let jac = [3.0, -3.0].repeat(4);
    let residual = [1.0, 2.0, 3.0, 4.0];
    let mut delta = [0.0; 4];
    let mut polls = 0;
    let result = workspace.solve(
        &nu,
        &jac,
        &[0.0; 4],
        &[0.0; 4],
        500.0,
        0.01,
        &residual,
        &mut delta,
        &mut || {
            polls += 1;
            polls == 3
        },
    );
    assert_eq!(result, Err(IsochoricKineticsError::Cancelled));
    workspace
        .solve(
            &nu,
            &jac,
            &[0.0; 4],
            &[0.0; 4],
            500.0,
            0.01,
            &residual,
            &mut delta,
            &mut || false,
        )
        .unwrap();
    let mut direct = [0.0; 4];
    NewtonWorkspace::new(2, 4, &nu)
        .unwrap()
        .solve(
            &nu,
            &jac,
            &[0.0; 4],
            &[0.0; 4],
            500.0,
            0.01,
            &residual,
            &mut direct,
            &mut || false,
        )
        .unwrap();
    assert_eq!(delta, direct);

    let mechanism = duplicate_isomers(7000.0, 0.0, 8.0, 6);
    let control = controls(1e-5);
    let mut prepared = IsochoricReactor::new(mechanism.clone(), control, || false).unwrap();
    let input = [Amount::new(1.0), Amount::new(0.0)];
    let storage = prepared.amounts.as_ptr();
    let mut polls = 0;
    let result = prepared.advance_adiabatic(
        &input,
        Volume::new(0.01),
        Temperature::new(500.0),
        Time::new(0.1),
        || {
            polls += 1;
            polls > 160
        },
    );
    assert!(matches!(result, Err(IsochoricKineticsError::Cancelled)));
    let recovered = prepared
        .advance_adiabatic(
            &input,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(0.1),
            || false,
        )
        .unwrap();
    assert_eq!(recovered.amounts().as_ptr(), storage);
    assert_eq!(
        recovered.to_owned(|| false).unwrap(),
        run(&mechanism, 500.0, 0.1, control)
    );
    assert_eq!(input, [Amount::new(1.0), Amount::new(0.0)]);
}

#[test]
fn sparse_gri_schur_matches_dense_original_equation_and_prepared_reuse() {
    use fs_thermochem::kinetics::gri30::{Gri30Calorics, cantera_gri30};
    let mechanism = cantera_gri30(Gri30Calorics::Original).unwrap();
    let ns = mechanism.thermodynamics().models().len();
    let nr = mechanism.reactions().len();
    let mut nu = vec![0.0; ns * nr];
    for i in 0..ns {
        for j in 0..nr {
            nu[i * nr + j] = mechanism.stoichiometry().get(i, j).unwrap() as f64;
        }
    }
    let mut sparse = NewtonWorkspace::new(ns, nr, &nu).unwrap();
    let mut dense = NewtonWorkspace::dense_reduced(ns, nr).unwrap();
    assert_eq!(sparse.dimension(), 54);
    assert_eq!(dense.dimension(), 54);
    let nonzero = nu.iter().filter(|&&x| x != 0.0).count();
    assert_eq!(sparse.nonzero_stoichiometry(), nonzero);
    assert!(nonzero < ns * nr / 8);
    let concentrations: Vec<_> = (0..ns).map(|i| 0.01 * (i + 1) as f64).collect();
    let mut kinetics = KineticWorkspaceV1::new(&mechanism).unwrap();
    let state = mechanism
        .evaluate(
            Temperature::new(1800.0),
            &concentrations,
            &mut kinetics,
            || false,
        )
        .unwrap();
    let mut dt = vec![0.0; nr];
    for j in 0..nr {
        let mut energy = Sum::default();
        for i in 0..ns {
            energy.add(kinetics.species_internal_energies()[i] * nu[i * nr + j]);
        }
        dt[j] = -energy.value() / state.thermodynamics.heat_capacity_v_density_j_m3_k;
    }
    let residual: Vec<_> = (0..nr).map(|j| 1e-7 * ((j % 17) as f64 - 8.0)).collect();
    let mut a = vec![0.0; nr];
    let mut b = vec![0.0; nr];
    for step in [1e-10, 1e-8, 1e-6] {
        for (workspace, output) in [(&mut sparse, &mut a), (&mut dense, &mut b)] {
            workspace
                .solve(
                    &nu,
                    kinetics.rate_jacobian(),
                    kinetics.temperature_derivatives(),
                    &dt,
                    1800.0,
                    step,
                    &residual,
                    output,
                    &mut || false,
                )
                .unwrap();
        }
        // Exact zeros contribute nothing to the compensated reductions; source
        // reaction order and the dense thermal coordinate remain unchanged.
        assert_eq!(a, b);
    }
    // A small source-specific release timing, without a timing assertion or a
    // benchmark framework. Both owners reuse all storage and solve identical
    // original equations; it measures this operation, not reactor throughput.
    #[cfg(not(debug_assertions))]
    {
        const REPEATS: usize = 40;
        let mut timings = [std::time::Duration::ZERO; 2];
        for (index, workspace) in [&mut sparse, &mut dense].into_iter().enumerate() {
            let start = std::time::Instant::now();
            for _ in 0..REPEATS {
                workspace
                    .solve(
                        &nu,
                        kinetics.rate_jacobian(),
                        kinetics.temperature_derivatives(),
                        &dt,
                        1800.0,
                        1e-8,
                        &residual,
                        &mut a,
                        &mut || false,
                    )
                    .unwrap();
                std::hint::black_box(&a);
            }
            timings[index] = start.elapsed();
        }
        eprintln!(
            "GRI Schur {REPEATS} prepared 54x54 solves: sparse={:?}, dense={:?}, N nonzeros={nonzero}/{}",
            timings[0],
            timings[1],
            ns * nr
        );
    }
}
