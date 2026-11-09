//! G0/G1: common series-creep mechanics and explicit state-rounding accounts.
use fs_material::fiber::WoolFelt;
use fs_material::fiber::felt_rate::{CreepBranch, CreepSpectrum, MAX_FELT_RATE_BRANCHES, RateFelt};

fn law() -> CreepSpectrum {
    CreepSpectrum::new(vec![
        CreepBranch {
            modulus_pa: 1e5,
            viscosity_pa_s: 2000.0,
        },
        CreepBranch {
            modulus_pa: 7e4,
            viscosity_pa_s: 14000.0,
        },
    ])
    .unwrap()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-300),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_original_held_stress_and_recovery_with_independent_heat() {
    let spectrum = law();
    for dt in [1e-8, 1e-4, 0.01, 0.1, 2.0] {
        let p = spectrum.prepare(dt).unwrap();
        let stress = 1500.0;
        let mut state = [0.0; MAX_FELT_RATE_BRANCHES];
        state[0] = 0.013;
        state[1] = 0.004;
        let trial = p.advance(&state, stress).unwrap();
        let mut independent_heat = 0.0;
        for (i, b) in spectrum.branches().iter().enumerate() {
            let exposure = dt * b.modulus_pa / b.viscosity_pa_s;
            let a = (-exposure).exp();
            let target = stress / b.modulus_pa;
            let expected = if exposure < 0.5 {
                state[i] + (-(-exposure).exp_m1()) * (target - state[i])
            } else {
                target + a * (state[i] - target)
            };
            close(trial.strain[i], expected, 1e-14);
            independent_heat +=
                0.5 * b.modulus_pa * (-(-2.0 * exposure).exp_m1()) * (state[i] - target).powi(2);
        }
        close(trial.dissipated_energy_j_m3, independent_heat, 1e-13);
        assert!(trial.balance_residual_j_m3.abs() <= trial.roundoff_bound_j_m3);
        let released = p.advance(&trial.strain, 0.0).unwrap();
        assert!(released.dissipated_energy_j_m3 >= 0.0);
        for (i, b) in spectrum.branches().iter().enumerate() {
            close(
                released.strain[i],
                trial.strain[i] * (-dt * b.modulus_pa / b.viscosity_pa_s).exp(),
                2e-14,
            );
        }
    }
}

#[test]
fn g0_large_exposure_retains_representable_recovery_tail() {
    let spectrum = CreepSpectrum::new(vec![CreepBranch {
        modulus_pa: 3.0,
        viscosity_pa_s: 3.0,
    }])
    .unwrap();
    let p = spectrum.prepare(100.0).unwrap();
    let mut state = [0.0; MAX_FELT_RATE_BRANCHES];
    state[0] = 0.2;
    let expected = 0.2 * (-100.0_f64).exp();
    let free = p.free_strain(&state).unwrap();
    let trial = p.advance(&state, 0.0).unwrap();
    assert!(free > 0.0 && trial.strain[0] > 0.0);
    close(free, expected, 1e-14);
    close(trial.strain[0], expected, 1e-14);
    close(trial.dissipated_energy_j_m3, 0.5 * 3.0 * 0.2 * 0.2, 1e-14);
}

#[test]
fn g0_equilibrated_long_holds_expose_only_bounded_roundoff() {
    let p = law().prepare(0.001).unwrap();
    let mut state = [0.0; MAX_FELT_RATE_BRANCHES];
    let mut work = 0.0;
    let mut heat = 0.0;
    let mut allowed = 0.0;
    for _ in 0..25000 {
        let trial = p.advance(&state, 1234.56789).unwrap();
        assert!(trial.dissipated_energy_j_m3 >= 0.0);
        assert!(trial.balance_residual_j_m3.abs() <= trial.roundoff_bound_j_m3);
        work += trial.work_density_j_m3;
        heat += trial.dissipated_energy_j_m3;
        allowed += trial.roundoff_bound_j_m3;
        state = trial.strain;
    }
    let stored = p.stored_energy_density(&state).unwrap();
    assert!((work - stored - heat).abs() <= allowed + 1e-11 * stored);
    close(state[0], 1234.56789 / 1e5, 2e-14);
    close(state[1], 1234.56789 / 7e4, 3e-14);
}

#[test]
fn g1_shared_kernel_matches_nonlinear_material_creep_trial() {
    let branches = law().branches().to_vec();
    let full = RateFelt::new(
        WoolFelt::new(30000.0, 0.2, 2.2, 3.0, 0.15, 0.7).unwrap(),
        branches.clone(),
    )
    .unwrap();
    let p = full.prepare(0.0007).unwrap();
    let creep = CreepSpectrum::new(branches)
        .unwrap()
        .prepare(0.0007)
        .unwrap();
    let mut state = full.initial_state();
    for i in 0..300 {
        let stress = if i < 180 {
            1000.0 + 200.0 * (i as f64 * 0.07).sin()
        } else {
            0.0
        };
        let expected = creep.advance(&state.creep_strain, stress).unwrap();
        let trial = p.step_stress_driven(&state, stress).unwrap();
        assert_eq!(expected.strain, trial.state.creep_strain);
        assert_eq!(
            expected.dissipated_energy_j_m3.to_bits(),
            trial.creep_dissipation_j_m3.to_bits()
        );
        state = trial.state;
    }
}

#[test]
fn g0_empty_zero_duration_and_invalid_histories_are_explicit() {
    let empty = CreepSpectrum::new(vec![]).unwrap().prepare(0.0).unwrap();
    let virgin = [0.0; MAX_FELT_RATE_BRANCHES];
    let zero = empty.advance(&virgin, 3.0).unwrap();
    assert_eq!(zero.strain, virgin);
    assert_eq!(zero.work_density_j_m3, 0.0);
    assert_eq!(zero.dissipated_energy_j_m3, 0.0);
    assert_eq!(empty.compliance_pa_inverse(), 0.0);
    let spectrum = law();
    let p = spectrum.prepare(0.0).unwrap();
    let mut state = virgin;
    state[0] = 0.02;
    state[1] = 0.03;
    let still = p.advance(&state, 1000.0).unwrap();
    assert_eq!(still.strain, state);
    assert_eq!(still.dissipated_energy_j_m3, 0.0);
    assert!(spectrum.prepare(-1.0).is_err());
    assert!(p.advance(&state, -1.0).is_err());
    assert!(p.advance(&state, f64::NAN).is_err());
    let mut invalid = state;
    invalid[7] = 0.1;
    assert!(p.advance(&invalid, 1.0).is_err());
    assert!(
        CreepSpectrum::new(vec![
            CreepBranch {
                modulus_pa: 1.0,
                viscosity_pa_s: 1.0
            };
            9
        ])
        .is_err()
    );
    assert_eq!(p.advance(&state, 1000.0).unwrap(), still);
}

#[test]
fn g1_source_realization_needs_no_nonlinear_envelope() {
    use fs_material::visco::GeneralizedMaxwell;
    let card = GeneralizedMaxwell::new(1.0, vec![(1e20, 1.0)]).unwrap();
    let spectrum = CreepSpectrum::from_prony(&card).unwrap();
    assert_eq!(spectrum.branches().len(), 1);
    close(spectrum.branches()[0].modulus_pa, 1.0, 1e-15);
    close(spectrum.branches()[0].viscosity_pa_s, 1e20, 1e-15);
}
