//! G0/G1: independent Norton/Arrhenius references and physical work accounting.
use fs_material::creep::{CreepState, CreepTolerances, NortonCreep, NortonCreepParameters};
use fs_material::elastic::IsotropicElastic;
use fs_material::tensor::{contract, rotate, rotation, trace, von_mises};
use fs_material::{SmallStrainLaw, Voigt};

fn parameters(n: f64) -> NortonCreepParameters {
    NortonCreepParameters {
        reference_stress_pa: 3.0e6,
        reference_rate_s: 2.0e-3,
        stress_exponent: n,
        activation_energy_j_mol: 0.0,
        reference_temperature_k: 500.0,
        temperature_range_k: [300.0, 900.0],
    }
}
fn elastic() -> IsotropicElastic {
    IsotropicElastic::new(1.2e9, 0.2, 0.05).unwrap()
}
fn law(n: f64, dt: f64) -> NortonCreep {
    NortonCreep::new(elastic(), parameters(n), 500.0, dt).unwrap()
}
fn tolerances() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 2.0e-6,
        energy_j_m3: 1.0e-7,
    }
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute,
        "{a} != {b}, tolerance {absolute}"
    );
}
fn deviatoric_strain(e: f64) -> Voigt {
    [e, -0.5 * e, -0.5 * e, 0.0, 0.0, 0.0]
}

#[test]
fn g1_linear_creep_has_no_yield_threshold_and_separates_heat_from_numerical_loss() {
    let law = law(1.0, 0.25);
    let (_, g) = law.elastic().lame();
    let strain = deviatoric_strain(0.004);
    let old = law.initial_state();
    let next = law
        .advance_with_energy(&[0.0; 6], &old, &strain, tolerances())
        .unwrap();
    let eta = parameters(1.0).reference_stress_pa / parameters(1.0).reference_rate_s;
    let delta = 3.0 * g * strain[0] / (3.0 * g + eta / law.duration_s());
    let q = eta * delta / law.duration_s();
    close(next.state.equivalent_creep_strain, delta, 2e-17);
    close(next.state.creep_strain[0], delta, 2e-17);
    close(next.state.creep_strain[1], -0.5 * delta, 2e-17);
    close(von_mises(&next.stress_pa), q, 2e-8);
    close(next.energy.physical_dissipation_j_m3, q * delta, 2e-10);
    let stored = 1.5 * g * (strain[0] - delta).powi(2);
    close(next.energy.final_stored_j_m3, stored, 2e-10);
    close(next.energy.numerical_dissipation_j_m3, stored, 2e-10);
    assert!(next.energy.residual_j_m3.abs() < 1e-8);
    assert_eq!(old, law.initial_state());
}

#[test]
fn g1_constant_uniaxial_stress_gives_linear_secondary_creep_and_permanent_set() {
    let law = law(3.0, 2.0);
    let sigma: f64 = 1.0e6;
    let elastic = [
        sigma / law.elastic().youngs,
        -law.elastic().poisson * sigma / law.elastic().youngs,
        -law.elastic().poisson * sigma / law.elastic().youngs,
        0.0,
        0.0,
        0.0,
    ];
    let rate =
        parameters(3.0).reference_rate_s * (sigma / parameters(3.0).reference_stress_pa).powi(3);
    let mut state = law.initial_state();
    let mut old_strain = elastic;
    for step in 1..=8 {
        let alpha = rate * step as f64 * law.duration_s();
        let creep = deviatoric_strain(alpha);
        let strain = core::array::from_fn(|i| elastic[i] + creep[i]);
        let next = law
            .advance_with_energy(&old_strain, &state, &strain, tolerances())
            .unwrap();
        close(next.stress_pa[0], sigma, 2e-8);
        for s in &next.stress_pa[1..] {
            close(*s, 0.0, 2e-8);
        }
        close(next.state.equivalent_creep_strain, alpha, 2e-17);
        close(
            next.energy.physical_dissipation_j_m3,
            sigma * rate * law.duration_s(),
            1e-9,
        );
        state = next.state;
        old_strain = strain;
    }
    // Stress-free unloading ends at the retained creep strain. It does not
    // recover like a finite-equilibrium-modulus generalized Maxwell solid.
    let permanent = law
        .advance_with_energy(&old_strain, &state, &state.creep_strain, tolerances())
        .unwrap();
    assert_eq!(permanent.state, state);
    assert!(permanent.stress_pa.iter().all(|s| s.abs() < 1e-12));
    close(permanent.energy.physical_dissipation_j_m3, 0.0, 0.0);
    assert!(state.creep_strain[0] > 0.001);
}

#[test]
fn g1_fixed_tensor_strain_relaxation_converges_to_the_analytic_norton_solution() {
    let (_, g) = elastic().lame();
    let strain = deviatoric_strain(0.005);
    let q0 = 3.0 * g * strain[0];
    let total_time = 0.4;
    for exponent in [1.0, 3.0] {
        let p = parameters(exponent);
        let rate = 3.0 * g * p.reference_rate_s / p.reference_stress_pa;
        let exact = if exponent == 1.0 {
            q0 * (-rate * total_time).exp()
        } else {
            q0 / (1.0
                + (exponent - 1.0)
                    * rate
                    * total_time
                    * (q0 / p.reference_stress_pa).powf(exponent - 1.0))
            .powf(1.0 / (exponent - 1.0))
        };
        let mut errors = Vec::new();
        for steps in [16, 32, 64] {
            let law = law(exponent, total_time / steps as f64);
            let mut state = law.initial_state();
            let mut heat = 0.0;
            let mut numerical = 0.0;
            let mut final_energy = 0.0;
            let mut q = 0.0;
            for _ in 0..steps {
                let next = law
                    .advance_with_energy(&strain, &state, &strain, tolerances())
                    .unwrap();
                assert_eq!(next.energy.mechanical_work_j_m3, 0.0);
                assert!(next.energy.physical_dissipation_j_m3 > 0.0);
                heat += next.energy.physical_dissipation_j_m3;
                numerical += next.energy.numerical_dissipation_j_m3;
                final_energy = next.energy.final_stored_j_m3;
                q = von_mises(&next.stress_pa);
                state = next.state;
            }
            let initial = 1.5 * g * strain[0].powi(2);
            close(initial - final_energy, heat + numerical, 1e-7);
            errors.push((q - exact).abs());
        }
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.2,
                "{errors:?}"
            );
        }
    }
}

#[test]
fn g0_tensor_shear_flow_and_zero_stress_linear_tangent_have_correct_factors() {
    let law = law(1.0, 0.5);
    let (_, g) = law.elastic().lame();
    let k = 3.0 * g * law.duration_s() * parameters(1.0).reference_rate_s
        / parameters(1.0).reference_stress_pa;
    let strain = [0.0, 0.0, 0.0, 0.003, 0.0, 0.0];
    let state = law.initial_state();
    let next = law
        .advance_with_energy(&[0.0; 6], &state, &strain, tolerances())
        .unwrap();
    close(next.stress_pa[3], 2.0 * g * strain[3] / (1.0 + k), 1e-8);
    close(next.state.creep_strain[3], k * strain[3] / (1.0 + k), 2e-17);
    close(
        next.state.equivalent_creep_strain,
        2.0 / 3.0_f64.sqrt() * next.state.creep_strain[3],
        2e-17,
    );
    close(
        law.tangent(&[0.0; 6], &state)[3][3],
        2.0 * g / (1.0 + k),
        1e-6,
    );
}

#[test]
fn g0_algorithmic_tangent_matches_all_tensor_columns_including_reversal() {
    for exponent in [1.0, 3.0] {
        let law = law(exponent, 0.2);
        let prestrain = [0.003, -0.001, -0.001, 0.0007, 0.0003, -0.0002];
        let state = law.update_state(&prestrain, &law.initial_state());
        for strain in [[0.0; 6], [-0.003, 0.001, 0.0, -0.001, 0.0008, 0.0005]] {
            let tangent = law.tangent(&strain, &state);
            for j in 0..6 {
                let mut plus = strain;
                let mut minus = strain;
                let h = 1e-8;
                plus[j] += h;
                minus[j] -= h;
                let sp = law.stress(&plus, &state);
                let sm = law.stress(&minus, &state);
                for i in 0..6 {
                    close((sp[i] - sm[i]) / (2.0 * h), tangent[i][j], 0.2);
                }
            }
        }
    }
}

#[test]
fn g0_arrhenius_uses_absolute_temperature_and_the_declared_reference_rate() {
    let mut p = parameters(3.0);
    p.activation_energy_j_mol = 60_000.0;
    let reference = NortonCreep::new(elastic(), p.clone(), 500.0, 0.1).unwrap();
    let warm = NortonCreep::new(elastic(), p.clone(), 700.0, 0.1).unwrap();
    let expected =
        (p.activation_energy_j_mol / 8.314_462_618_153_24 * (1.0 / 500.0 - 1.0 / 700.0)).exp();
    close(
        reference
            .equivalent_creep_rate_s(p.reference_stress_pa)
            .unwrap(),
        p.reference_rate_s,
        1e-17,
    );
    close(
        warm.equivalent_creep_rate_s(p.reference_stress_pa).unwrap()
            / reference
                .equivalent_creep_rate_s(p.reference_stress_pa)
                .unwrap(),
        expected,
        2e-12,
    );
    let strain = deviatoric_strain(0.002);
    let hot = warm
        .advance_with_energy(&strain, &warm.initial_state(), &strain, tolerances())
        .unwrap();
    let cold = reference
        .advance_with_energy(&strain, &reference.initial_state(), &strain, tolerances())
        .unwrap();
    assert!(hot.state.equivalent_creep_strain > cold.state.equivalent_creep_strain);
    assert!(von_mises(&hot.stress_pa) < von_mises(&cold.stress_pa));
}

#[test]
fn g0_hydrostatic_loading_does_not_creep_and_rotations_preserve_response() {
    let law = law(3.0, 0.2);
    let state = law.initial_state();
    let hydro = [0.001_234_567; 3]
        .into_iter()
        .chain([0.0; 3])
        .collect::<Vec<_>>();
    let hydro: Voigt = hydro.try_into().unwrap();
    let next = law
        .advance_with_energy(&[0.0; 6], &state, &hydro, tolerances())
        .unwrap();
    assert_eq!(next.state, state);
    assert_eq!(next.energy.physical_dissipation_j_m3, 0.0);
    let strain = [0.003, -0.001, 0.0005, 0.001, -0.0003, 0.0002];
    let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.7);
    let first = law
        .advance_with_energy(&[0.0; 6], &state, &strain, tolerances())
        .unwrap();
    let second = law
        .advance_with_energy(&[0.0; 6], &state, &rotate(&strain, &q), tolerances())
        .unwrap();
    for (a, b) in rotate(&first.stress_pa, &q).iter().zip(second.stress_pa) {
        close(*a, b, 3e-8);
    }
    close(
        first.state.equivalent_creep_strain,
        second.state.equivalent_creep_strain,
        1e-17,
    );
    close(
        first.energy.physical_dissipation_j_m3,
        second.energy.physical_dissipation_j_m3,
        1e-9,
    );
}

#[test]
fn g0_reversal_retains_signed_history_and_monotone_accumulation() {
    let law = law(3.0, 0.2);
    let initial = law.initial_state();
    let strain = deviatoric_strain(0.004);
    let first = law
        .advance_with_energy(&[0.0; 6], &initial, &strain, tolerances())
        .unwrap();
    let reverse = deviatoric_strain(-0.004);
    let second = law
        .advance_with_energy(&strain, &first.state, &reverse, tolerances())
        .unwrap();
    assert!(second.state.creep_strain[0] < first.state.creep_strain[0]);
    assert!(second.state.equivalent_creep_strain > first.state.equivalent_creep_strain);
    assert!(second.energy.physical_dissipation_j_m3 > 0.0);
    close(trace(&second.state.creep_strain), 0.0, 1e-18);
    assert!(
        second.state.equivalent_creep_strain
            >= (2.0 / 3.0 * contract(&second.state.creep_strain, &second.state.creep_strain))
                .sqrt()
    );
}

#[test]
fn g0_nearly_cancelled_reversal_preserves_exact_deviatoric_history() {
    let law = law(1.0, 0.2);
    let strain = [0.004, -0.001, -0.003, 0.0008, -0.0003, 0.0005];
    let first = law
        .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, tolerances())
        .unwrap();
    let (_, g) = law.elastic().lame();
    let k = 3.0 * g * law.duration_s() * parameters(1.0).reference_rate_s
        / parameters(1.0).reference_stress_pa;
    let reversed = core::array::from_fn(|i| -first.state.creep_strain[i] / k);
    let second = law
        .advance_with_energy(&strain, &first.state, &reversed, tolerances())
        .unwrap();
    assert_eq!(trace(&second.state.creep_strain), 0.0);
    assert!(second.state.creep_strain.iter().all(|e| e.abs() < 1e-17));
    assert!(second.state.equivalent_creep_strain > first.state.equivalent_creep_strain);
}

#[test]
fn g0_near_zero_load_newton_trials_admit_history_roundoff() {
    let law = law(3.0, 0.2);
    let loaded = law.update_state(&deviatoric_strain(0.004), &law.initial_state());
    let perturbation = deviatoric_strain(1e-17);
    let unloaded = core::array::from_fn(|i| loaded.creep_strain[i] + perturbation[i]);
    let next = law
        .advance_with_energy(&unloaded, &loaded, &unloaded, tolerances())
        .unwrap();
    assert_eq!(next.state, loaded);
    assert!(
        next.stress_pa
            .iter()
            .all(|x| x.is_finite() && x.abs() < 1e-7)
    );
    assert!(next.energy.physical_dissipation_j_m3 >= 0.0);
    assert!(next.energy.residual_j_m3.abs() < 1e-30);
}

#[test]
fn g0_invalid_inputs_refuse_without_touching_accepted_history() {
    for (temperature, dt) in [(299.0, 1.0), (901.0, 1.0), (500.0, 0.0), (f64::NAN, 1.0)] {
        assert!(NortonCreep::new(elastic(), parameters(3.0), temperature, dt).is_err());
    }
    assert!(NortonCreep::new(elastic(), parameters(0.5), 500.0, 1.0).is_err());
    let law = law(3.0, 0.2);
    let state = law.initial_state();
    let preserved = state.clone();
    for strain in [[0.1; 6], [f64::NAN; 6]] {
        assert!(
            law.advance_with_energy(&[0.0; 6], &state, &strain, tolerances())
                .is_err()
        );
        assert_eq!(state, preserved);
    }
    let invalid = CreepState {
        creep_strain: [0.001; 6],
        equivalent_creep_strain: 0.0,
    };
    assert!(law.validate_history(&[0.0; 6], &invalid).is_err());
    assert!(law.stress(&[0.0; 6], &invalid).iter().all(|x| x.is_nan()));
    assert_eq!(state, preserved);
}

#[test]
fn g0_logarithmic_return_avoids_overflowing_raw_stress_powers() {
    let mut p = parameters(100.0);
    p.reference_stress_pa = 1.0e6;
    let law = NortonCreep::new(elastic(), p, 500.0, 0.1).unwrap();
    let strain = deviatoric_strain(0.004);
    let next = law
        .advance_with_energy(
            &[0.0; 6],
            &law.initial_state(),
            &strain,
            CreepTolerances {
                stress_pa: 2e-5,
                energy_j_m3: 1e-6,
            },
        )
        .unwrap();
    assert!(next.stress_pa.iter().all(|x| x.is_finite()));
    assert!(next.state.equivalent_creep_strain > 0.0);
    assert!(von_mises(&next.stress_pa) < 1.1e6);
}
