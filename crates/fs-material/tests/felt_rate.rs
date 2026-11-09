//! G0/G1: original felt parity, exact Kelvin loss and rate-history refinement.
use fs_material::fiber::felt_rate::{CreepBranch, MAX_FELT_RATE_BRANCHES, RateFelt};
use fs_material::fiber::{Uniaxial, WoolFelt};

fn envelope() -> WoolFelt {
    WoolFelt::new(30000.0, 0.2, 2.2, 3.0, 0.15, 0.7).unwrap()
}
fn conservative() -> WoolFelt {
    WoolFelt::new(20000.0, 0.2, 2.0, 2.0, 0.0, 0.7).unwrap()
}
fn arm() -> CreepBranch {
    CreepBranch {
        modulus_pa: 100000.0,
        viscosity_pa_s: 2000.0,
    }
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-100),
        "{a:e} != {b:e}, relative tolerance {relative:e}"
    );
}

#[test]
fn g0_zero_branches_is_bitwise_woolfelt() {
    let wool = envelope();
    let law = RateFelt::new(wool.clone(), vec![])
        .unwrap()
        .prepare(0.003)
        .unwrap();
    let mut state = RateFelt::new(wool.clone(), vec![]).unwrap().initial_state();
    let mut history = wool.initial_state();
    let mut seed = 0x85c431ab1297de61_u64;
    for _ in 0..10000 {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let strain = (seed >> 11) as f64 / (1_u64 << 53) as f64 * 0.7 - 0.05;
        let expected = wool.response(strain, &history).unwrap();
        let trial = law.step_strain_driven(&state, strain).unwrap();
        assert_eq!(trial.stress_pa.to_bits(), expected.stress_pa.to_bits());
        assert_eq!(trial.tangent_pa.to_bits(), expected.tangent_pa.to_bits());
        assert_eq!(trial.state.hysteresis, expected.state);
        assert_eq!(trial.creep_dissipation_j_m3, 0.0);
        assert!(trial.crush_dissipation_j_m3 >= 0.0);
        state = trial.state;
        history = expected.state;
    }
}

#[test]
fn g1_held_stress_creep_and_free_recovery_have_independent_exact_heat() {
    let law = RateFelt::new(conservative(), vec![arm()]).unwrap();
    let prepared = law.prepare(0.03).unwrap();
    let initial = law.initial_state();
    let stress = 2000.0;
    let result = prepared.step_stress_driven(&initial, stress).unwrap();
    let a = (-0.03_f64 / 0.02).exp();
    let z = stress / 100000.0 * (1.0 - a);
    let heat = 0.5 * 100000.0 * (stress / 100000.0_f64).powi(2) * (1.0 - a * a);
    close(result.state.creep_strain[0], z, 5e-15);
    close(result.creep_dissipation_j_m3, heat, 5e-15);
    assert_eq!(result.crush_dissipation_j_m3, 0.0);
    close(
        result.state.elastic_strain,
        0.2 * (stress / 20000.0_f64).sqrt(),
        5e-15,
    );
    close(
        result.work_density_j_m3,
        result.state.elastic_strain * stress / 3.0 + stress * z,
        1e-14,
    );
    let recovery = prepared.step_stress_driven(&result.state, 0.0).unwrap();
    close(recovery.state.creep_strain[0], z * a, 1e-14);
    close(
        recovery.creep_dissipation_j_m3,
        0.5 * 100000.0 * z * z * (1.0 - a * a),
        1e-14,
    );
    assert_eq!(recovery.stress_pa, 0.0);
    assert_eq!(recovery.state.elastic_strain, 0.0);
    assert!(recovery.state.creep_strain[0] > 0.0);
    close(
        recovery.state.work,
        recovery.state.stored + recovery.state.dissipated,
        1e-13,
    );
    assert_eq!(initial, law.initial_state());
}

#[test]
fn g0_work_equals_stored_plus_dissipated_and_tangent_matches_same_trial() {
    let law = RateFelt::new(
        envelope(),
        vec![
            arm(),
            CreepBranch {
                modulus_pa: 70000.0,
                viscosity_pa_s: 14000.0,
            },
        ],
    )
    .unwrap();
    let prepared = law.prepare(0.001).unwrap();
    let mut state = law.initial_state();
    for i in 0..240 {
        let strain = if i < 80 {
            0.3 * (i + 1) as f64 / 80.0
        } else if i < 120 {
            0.3
        } else {
            0.3 - 0.35 * (i - 119) as f64 / 120.0
        };
        let result = prepared.step_strain_driven(&state, strain).unwrap();
        let loss = result.creep_dissipation_j_m3 + result.crush_dissipation_j_m3;
        assert!(loss >= 0.0);
        // Work may nearly cancel storage release and heat. Scale the first-law
        // residual by the physical exchanged terms, not their small net.
        let scale = result
            .work_density_j_m3
            .abs()
            .max(result.stored_energy_change_j_m3.abs())
            .max(loss);
        assert!(
            (result.work_density_j_m3 - result.stored_energy_change_j_m3 - loss).abs()
                <= 1e-12 * scale
        );
        close(
            result.state.work,
            result.state.stored + result.state.dissipated,
            1e-12,
        );
        if i == 50 || i == 100 || i == 160 {
            let delta = 1e-7;
            let plus = prepared.step_strain_driven(&state, strain + delta).unwrap();
            let minus = prepared.step_strain_driven(&state, strain - delta).unwrap();
            let fd = (plus.stress_pa - minus.stress_pa) / (2.0 * delta);
            close(result.tangent_pa, fd, 2e-6);
        }
        state = result.state;
    }
    assert!(state.dissipated > 0.0 && state.creep_strain[0] > 0.0);
}

fn history(t: f64) -> f64 {
    if t <= 0.04 {
        0.25 * t / 0.04
    } else if t <= 0.06 {
        0.25
    } else {
        0.25 * (0.1 - t) / 0.04
    }
}
fn original_derivative(t: f64, z: f64) -> f64 {
    let elastic = (history(t) - z).max(0.0);
    let stress = 20000.0 * (elastic / 0.2).powi(2);
    (stress - 100000.0 * z) / 2000.0
}
#[test]
fn g1_ramp_hold_release_has_first_order_time_refinement() {
    // Independent continuous original scalar constitutive ODE, with all
    // prescribed-history corners aligned to the reference integration grid.
    let intervals = 128000;
    let h = 0.1 / intervals as f64;
    let mut reference = 0.0;
    for i in 0..intervals {
        let t = i as f64 * h;
        let a = original_derivative(t, reference);
        let b = original_derivative(t + h / 2.0, reference + h * a / 2.0);
        let c = original_derivative(t + h / 2.0, reference + h * b / 2.0);
        let d = original_derivative(t + h, reference + h * c);
        reference += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let law = RateFelt::new(conservative(), vec![arm()]).unwrap();
    let mut errors = vec![];
    for steps in [100, 200, 400] {
        let h = 0.1 / steps as f64;
        let prepared = law.prepare(h).unwrap();
        let mut state = law.initial_state();
        for i in 0..steps {
            state = prepared
                .step_strain_driven(&state, history((i + 1) as f64 * h))
                .unwrap()
                .state;
        }
        errors.push((state.creep_strain[0] - reference).abs());
    }
    for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
        let order = ratio.log2();
        assert!(
            (order - 1.0).abs() < 0.2,
            "errors={errors:?}, order={order}"
        );
    }
}

#[test]
fn g1_faster_loading_retains_greater_apparent_stiffness() {
    let law = RateFelt::new(conservative(), vec![arm()]).unwrap();
    let mut stresses = vec![];
    for duration in [0.0001, 0.001, 0.01, 0.1] {
        let prepared = law.prepare(duration / 256.0).unwrap();
        let mut state = law.initial_state();
        let mut stress = 0.0;
        for i in 0..256 {
            let result = prepared
                .step_strain_driven(&state, 0.2 * (i + 1) as f64 / 256.0)
                .unwrap();
            stress = result.stress_pa;
            state = result.state;
        }
        stresses.push(stress);
    }
    assert!(stresses.windows(2).all(|s| s[0] > s[1]), "{stresses:?}");
}

#[test]
fn g0_refusals_preserve_original_history_and_native_limits() {
    assert!(RateFelt::new(envelope(), vec![arm(); MAX_FELT_RATE_BRANCHES]).is_ok());
    assert!(RateFelt::new(envelope(), vec![arm(); MAX_FELT_RATE_BRANCHES + 1]).is_err());
    assert!(
        RateFelt::new(
            envelope(),
            vec![CreepBranch {
                modulus_pa: -1.0,
                ..arm()
            }]
        )
        .is_err()
    );
    let law = RateFelt::new(envelope(), vec![arm()]).unwrap();
    assert!(law.prepare(-0.1).is_err());
    let state = law.initial_state();
    let prepared = law.prepare(0.001).unwrap();
    let expected = prepared.step_strain_driven(&state, 0.2).unwrap();
    assert!(
        prepared
            .step_strain_driven(&state, 0.7 * (1.0 + 1e-12))
            .is_err()
    );
    assert!(prepared.step_stress_driven(&state, -1.0).is_err());
    let mut bad = state.clone();
    bad.creep_strain[7] = 0.1;
    assert!(prepared.step_strain_driven(&bad, 0.2).is_err());
    assert_eq!(prepared.step_strain_driven(&state, 0.2).unwrap(), expected);
    let empty = RateFelt::new(envelope(), vec![]).unwrap();
    assert!(
        empty
            .prepare(0.0)
            .unwrap()
            .step_strain_driven(&empty.initial_state(), 0.7)
            .is_ok()
    );
}

#[test]
fn g0_repeated_peak_stress_preserves_original_compaction_bits() {
    let wool = envelope();
    let law = RateFelt::new(wool.clone(), vec![]).unwrap();
    let prepared = law.prepare(0.001).unwrap();
    for strain in [0.1, 0.23, 0.37, 0.51] {
        let loaded = prepared
            .step_strain_driven(&law.initial_state(), strain)
            .unwrap();
        let repeated = prepared
            .step_stress_driven(&loaded.state, loaded.stress_pa)
            .unwrap();
        assert_eq!(repeated.state.hysteresis, loaded.state.hysteresis);
        assert_eq!(
            repeated.state.elastic_strain.to_bits(),
            loaded.state.elastic_strain.to_bits()
        );
        assert_eq!(repeated.crush_dissipation_j_m3, 0.0);
        assert_eq!(repeated.work_density_j_m3, 0.0);
    }
}
