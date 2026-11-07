//! G0/G1: timed dwell and simultaneous recovery/plasticity on original storage.

use super::*;
use fs_material::plastic::chaboche::rate::ChabocheOverstressJ2;
use fs_material::plastic::chaboche::recovery::ChabocheStaticRecoveryJ2;

#[test]
fn g1_subyield_dwell_releases_branch_energy_without_inventing_plastic_flow() {
    let original = law();
    let loaded = original
        .response(&deviatoric(0.012), &original.initial_state())
        .unwrap();
    let strain = loaded.state.plastic.plastic_strain;
    let unloaded = original.response(&strain, &loaded.state).unwrap();
    assert_eq!(unloaded.state, loaded.state);
    let old = unloaded.state;
    let rates = [2.0, 0.4];
    let duration = 0.7;
    let timed = ChabocheStaticRecoveryJ2::new(original.clone(), rates, duration).unwrap();
    let retained = old.clone();
    let next = timed
        .advance_with_energy(&strain, &old, &strain, TOL)
        .unwrap();
    assert_eq!(
        old, retained,
        "an interval proposal cannot mutate its accepted prefix"
    );
    assert_eq!(
        next.state.plastic, old.plastic,
        "static recovery is not a fabricated creep increment"
    );
    assert_eq!(next.stress_pa, [0.0; 6]);
    let mut released = 0.0;
    let mut physical = 0.0;
    let mut numerical = 0.0;
    for (i, c) in [800.0, 250.0].into_iter().enumerate() {
        let divisor = 1.0 + rates[i] * duration;
        let expected = old.backstress_pa[i].map(|x| x / divisor);
        for (a, b) in next.state.backstress_pa[i].iter().zip(expected) {
            near(*a, b, 1e-14);
        }
        let storage0 = 0.75 * contract(&old.backstress_pa[i], &old.backstress_pa[i]) / c;
        let storage1 = 0.75 * contract(&expected, &expected) / c;
        let difference = core::array::from_fn(|j| expected[j] - old.backstress_pa[i][j]);
        released += storage0 - storage1;
        physical += 2.0 * rates[i] * duration * storage1;
        numerical += 0.75 * contract(&difference, &difference) / c;
    }
    assert!(physical > 0.0 && numerical > 0.0);
    assert_eq!(next.energy.mechanical_work_j_m3, 0.0);
    near(-next.energy.stored_change_j_m3(), released, 1e-14);
    near(next.recovery_dissipation_j_m3, physical, 1e-14);
    near(next.energy.physical_dissipation_j_m3, physical, 1e-14);
    near(next.energy.numerical_dissipation_j_m3, numerical, 1e-14);
    near(released, physical + numerical, 1e-14);
    assert_eq!(
        next,
        timed
            .advance_with_energy(&strain, &old, &strain, TOL)
            .unwrap()
    );

    // Independent continuous dwell limit; numerical BE loss decreases as dt
    // resolves the actual interval. No new plastic history is consumed.
    let exact = old.backstress_pa[0][0] * (-rates[0] * duration).exp();
    let mut last_error = f64::INFINITY;
    let mut last_numerical = f64::INFINITY;
    for subdivisions in [1, 2, 4, 8, 16] {
        let interval = ChabocheStaticRecoveryJ2::new(
            original.clone(),
            rates,
            duration / f64::from(subdivisions),
        )
        .unwrap();
        let mut state = old.clone();
        let mut loss = 0.0;
        for _ in 0..subdivisions {
            let step = interval
                .advance_with_energy(&strain, &state, &strain, TOL)
                .unwrap();
            loss += step.energy.numerical_dissipation_j_m3;
            state = step.state;
        }
        let error = (state.backstress_pa[0][0] - exact).abs();
        assert!(error < last_error);
        assert!(loss < last_numerical);
        assert_eq!(state.plastic, old.plastic);
        last_error = error;
        last_numerical = loss;
    }
}

#[test]
fn g1_joint_recovery_and_perzyna_reuse_full_nonproportional_tangent_and_heat() {
    let original = law();
    let old_strain = deviatoric(0.012);
    let old = original
        .response(&old_strain, &original.initial_state())
        .unwrap()
        .state;
    let strain = [-0.010, 0.005, 0.005, 0.006, -0.004, 0.002];
    let duration = 0.2;
    let rates = [3.0, 0.5];
    for kinetics in [
        None,
        Some(ChabocheOverstressJ2::new_linear(original.clone(), 60.0, duration).unwrap()),
        Some(
            ChabocheOverstressJ2::new_power_law(original.clone(), 60.0, 10.0, 2.5, duration)
                .unwrap(),
        ),
    ] {
        let timed = match &kinetics {
            Some(rate) => ChabocheStaticRecoveryJ2::with_overstress(rate.clone(), rates),
            None => ChabocheStaticRecoveryJ2::new(original.clone(), rates, duration),
        }
        .unwrap();
        let response = timed.response(&strain, &old).unwrap();
        let checked = timed
            .advance_with_energy(&old_strain, &old, &strain, TOL)
            .unwrap();
        assert_eq!(checked.state, response.state);
        let delta = checked.state.plastic.alpha - old.plastic.alpha;
        assert!(delta > 0.0);
        let dp: Voigt = core::array::from_fn(|j| {
            checked.state.plastic.plastic_strain[j] - old.plastic.plastic_strain[j]
        });
        let mut recovery = 0.0;
        for (i, (c, gamma)) in [(800.0, 150.0), (250.0, 20.0)].into_iter().enumerate() {
            for j in 0..6 {
                near(
                    checked.state.backstress_pa[i][j],
                    (old.backstress_pa[i][j] + (2.0 / 3.0) * c * dp[j])
                        / (1.0 + gamma * delta + rates[i] * duration),
                    1e-13,
                );
            }
            let a = checked.state.backstress_pa[i];
            recovery += 1.5 * (gamma * delta + rates[i] * duration) * contract(&a, &a) / c;
        }
        near(checked.recovery_dissipation_j_m3, recovery, 1e-13);
        let viscous = kinetics
            .as_ref()
            .map_or(0.0, |rate| rate.endpoint_overstress_pa(delta).unwrap());
        near(
            checked.energy.physical_dissipation_j_m3,
            3.0 * delta + recovery + viscous * delta,
            1e-13,
        );
        near(
            checked.energy.final_state.total_j_m3(),
            original
                .stored_energy(&strain, &checked.state)
                .unwrap()
                .total_j_m3(),
            1e-13,
        );
        for j in 0..6 {
            let mut plus = strain;
            let mut minus = strain;
            plus[j] += 1e-7;
            minus[j] -= 1e-7;
            let a = timed.response(&plus, &old).unwrap();
            let b = timed.response(&minus, &old).unwrap();
            for i in 0..6 {
                near(
                    response.tangent_pa[i][j],
                    (a.stress_pa[i] - b.stress_pa[i]) / 2e-7,
                    5e-6,
                );
            }
        }

        let zero = match &kinetics {
            Some(rate) => ChabocheStaticRecoveryJ2::with_overstress(rate.clone(), [0.0; 2]),
            None => ChabocheStaticRecoveryJ2::new(original.clone(), [0.0; 2], duration),
        }
        .unwrap();
        let expected = match &kinetics {
            Some(rate) => rate.response(&strain, &old),
            None => original.response(&strain, &old),
        }
        .unwrap();
        assert_eq!(zero.response(&strain, &old).unwrap(), expected);
        let expected_energy = match &kinetics {
            Some(rate) => rate.advance_with_energy(&old_strain, &old, &strain, TOL),
            None => original.advance_with_energy(&old_strain, &old, &strain, TOL),
        }
        .unwrap();
        assert_eq!(
            zero.advance_with_energy(&old_strain, &old, &strain, TOL)
                .unwrap(),
            expected_energy
        );
    }
}

#[test]
fn g0_recovery_preserves_absent_branches_and_refuses_invalid_original_history() {
    let absent =
        ChabocheJ2::new(backbone(60.0), [BackstressBranch::new(0.0, 100.0).unwrap()]).unwrap();
    let timed = ChabocheStaticRecoveryJ2::new(absent.clone(), [20.0], 0.3).unwrap();
    let initial = absent.initial_state();
    let strain = deviatoric(0.012);
    let next = timed
        .advance_with_energy(&[0.0; 6], &initial, &strain, TOL)
        .unwrap();
    let native = absent
        .advance_with_energy(&[0.0; 6], &initial, &strain, TOL)
        .unwrap();
    assert_eq!(next.state, native.state);
    assert_eq!(next.state.backstress_pa, [[0.0; 6]]);
    near(
        next.energy.final_state.total_j_m3(),
        native.energy.final_state.total_j_m3(),
        1e-14,
    );
    assert_eq!(next.recovery_dissipation_j_m3, 0.0);

    let original = law();
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(ChabocheStaticRecoveryJ2::new(original.clone(), [1.0; 2], duration).is_err());
    }
    for rate in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(ChabocheStaticRecoveryJ2::new(original.clone(), [rate, 1.0], 0.1).is_err());
    }
    assert!(ChabocheStaticRecoveryJ2::new(original.clone(), [1e300; 2], 1e300).is_err());
    let timed = ChabocheStaticRecoveryJ2::new(original.clone(), [1e6; 2], 1.0).unwrap();
    let mut invalid = original.initial_state();
    invalid.backstress_pa[0] = deviatoric(0.001);
    let retained = invalid.clone();
    assert!(timed.response(&[0.0; 6], &invalid).is_err());
    assert_eq!(
        invalid, retained,
        "preparation cannot repair an invalid incoming checkpoint"
    );
}
