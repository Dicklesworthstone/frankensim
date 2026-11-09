//! Affine loading against the original tensor differential/algebraic equations.
use super::semidefinite::dae_rate;
use super::*;

const PLANES: [MaterialPlane; 4] = [
    MaterialPlane::Axes12,
    MaterialPlane::Axes21,
    MaterialPlane::Axes13,
    MaterialPlane::Axes23,
];

fn engineering(history: &TensorPronyHistory, plane: MaterialPlane) -> [f64; 3] {
    let inside = axes(plane).0;
    let e = history.strain();
    [e[inside[0]], e[inside[1]], 2.0 * e[inside[2]]]
}

// Independent physical strain/Maxwell-strain state, with the existing original
// viscosity/algebraic traction solve; no production eigenmodes are consulted.
fn integrate_ramp(
    ode: &OriginalOde,
    mut state: [f64; 20],
    plane: MaterialPlane,
    target: [f64; 3],
    duration: f64,
) -> [f64; 20] {
    state[18] = 0.0;
    state[19] = 0.0;
    let end = mandel(prescribed(plane, target));
    let mut imposed = [0.0; 6];
    for i in axes(plane).0 {
        imposed[i] = (end[i] - state[i]) / duration;
    }
    let derivative = |state: &[f64; 20]| {
        let rate = dae_rate(ode, state, plane, imposed);
        let mut derivative = [0.0; 20];
        derivative[..6].copy_from_slice(&rate);
        for arm in 0..2 {
            for i in 0..6 {
                derivative[6 + 6 * arm + i] = rate[i] - state[6 + 6 * arm + i] / TAUS[arm];
            }
        }
        derivative[18] = dot(ode.stress(state, rate), rate);
        derivative[19] = ode.heat_rate(state, rate);
        derivative
    };
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let k1 = derivative(&state);
        let k2 = derivative(&core::array::from_fn(|i| state[i] + 0.5 * dt * k1[i]));
        let k3 = derivative(&core::array::from_fn(|i| state[i] + 0.5 * dt * k2[i]));
        let k4 = derivative(&core::array::from_fn(|i| state[i] + dt * k3[i]));
        for i in 0..20 {
            state[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    for i in axes(plane).0 {
        state[i] = end[i];
    }
    state
}

#[test]
fn affine_ramps_holds_and_reversals_match_original_ode_and_dae_in_every_plane() {
    for viscosity in [viscosity(), diagonal([0.0, 0.0, 0.0, 0.6, 0.9, 1.3])] {
        let (law, ode) = fixture(viscosity);
        for plane in PLANES {
            let mut history = law.history();
            let control = Control::new(
                &history,
                plane,
                TensorKelvinStressControlBudget::default(),
                || false,
            )
            .unwrap();
            let mut reference = [0.0; 20];
            for (target, duration) in [
                ([0.0012, -0.0007, 0.0015], 0.23),
                ([0.0012, -0.0007, 0.0015], 0.31),
                ([-0.0004, 0.0008, -0.0005], 0.17),
            ] {
                let old = reference;
                let old_strain = engineering(&history, plane);
                let step = control
                    .advance(&mut history, target, duration, || false)
                    .unwrap();
                reference = integrate_ramp(&ode, reference, plane, target, duration);
                let prescribed_rate =
                    core::array::from_fn(|i| (target[i] - old_strain[i]) / duration);
                let rate = dae_rate(
                    &ode,
                    &reference,
                    plane,
                    mandel(prescribed(plane, prescribed_rate)),
                );
                ode.compare(&history, &reference, rate);
                near(step.work, reference[18], 3e-13);
                near(step.dissipated, reference[19], 3e-13);
                near(
                    step.stored_energy_change,
                    ode.stored(&reference) - ode.stored(&old),
                    3e-13,
                );
                near(step.balance_residual, 0.0, 3e-13);
                assert!(step.dissipated > 0.0);
                for i in 0..6 {
                    near(step.strain_rate[i], rate[i] / W[i], 3e-11);
                    assert!(step.strain[i] >= step.strain_bounds[0][i] - 1e-15);
                    assert!(step.strain[i] <= step.strain_bounds[1][i] + 1e-15);
                }
                for i in axes(plane).1 {
                    assert_eq!(step.stress[i], 0.0);
                    assert_eq!(step.mean_stress[i], 0.0);
                }
                let mean_work = (0..3)
                    .map(|j| {
                        let i = axes(plane).0[j];
                        step.mean_stress[i] * (target[j] - old_strain[j])
                    })
                    .sum::<f64>();
                near(step.work, mean_work, 3e-17);
            }
            // Original full-strain loading can resume with the same memory.
            let old = mandel(history.strain());
            let target = old.map(|v| 0.8 * v);
            let duration = 0.13;
            history
                .advance(core::array::from_fn(|i| target[i] / W[i]), duration, || {
                    false
                })
                .unwrap();
            for (arm, tau) in TAUS.into_iter().enumerate() {
                let r = duration / tau;
                let decay = (-r).exp();
                let phi = -(-r).exp_m1() / r;
                for i in 0..6 {
                    reference[6 + 6 * arm + i] =
                        decay * reference[6 + 6 * arm + i] + phi * (target[i] - old[i]);
                }
            }
            reference[..6].copy_from_slice(&target);
            ode.compare(
                &history,
                &reference,
                core::array::from_fn(|i| (target[i] - old[i]) / duration),
            );
        }
    }
}

#[test]
fn consistent_endpoint_and_mean_tangents_include_normal_rate_and_all_shear_conventions() {
    let viscosities = [
        viscosity(),
        diagonal([0.0, 0.0, 0.0, 0.6, 0.9, 1.3]),
        diagonal([0.7, 0.4, 1.1, 0.6, 0.0, 0.0]),
        diagonal([0.7, 0.4, 0.0, 0.6, 0.0, 0.0]),
    ];
    for viscosity in viscosities {
        let (law, _) = fixture(viscosity);
        for plane in PLANES {
            let mut history = law.history();
            let control = Control::new(
                &history,
                plane,
                TensorKelvinStressControlBudget::default(),
                || false,
            )
            .unwrap();
            control
                .advance(&mut history, [0.001, 0.0004, -0.0007], 0.21, || false)
                .unwrap();
            let original = history.clone();
            let duration = 0.17;
            let tangent = control.ramp_tangent(duration, || false).unwrap();
            for col in 0..3 {
                let h = 1e-6;
                let mut plus = [0.0004, -0.0003, 0.0002];
                let mut minus = plus;
                plus[col] += h;
                minus[col] -= h;
                let positive = control
                    .advance(&mut history.clone(), plus, duration, || false)
                    .unwrap();
                let negative = control
                    .advance(&mut history.clone(), minus, duration, || false)
                    .unwrap();
                for row in 0..6 {
                    for (analytic, actual) in [
                        (
                            tangent.strain[row][col],
                            (positive.strain[row] - negative.strain[row]) / (2.0 * h),
                        ),
                        (
                            tangent.stress[row][col],
                            (positive.stress[row] - negative.stress[row]) / (2.0 * h),
                        ),
                        (
                            tangent.mean_stress[row][col],
                            (positive.mean_stress[row] - negative.mean_stress[row]) / (2.0 * h),
                        ),
                        (
                            tangent.strain_rate[row][col],
                            (positive.strain_rate[row] - negative.strain_rate[row]) / (2.0 * h),
                        ),
                    ] {
                        near(actual, analytic, 2e-8 * analytic.abs().max(1.0));
                    }
                }
            }
            assert_eq!(history, original);
        }
    }
}

#[test]
fn exact_subdivision_tiny_holds_and_algebraic_kelvin_limits_preserve_physical_energy() {
    let plane = MaterialPlane::Axes12;
    let (law, ode) = fixture(viscosity());
    let (history, state) = seed(&law, plane, BIAS);
    let control = Control::new(
        &history,
        plane,
        TensorKelvinStressControlBudget::default(),
        || false,
    )
    .unwrap();
    let target = [-0.0004, 0.0003, 0.0011];
    let mut whole = history.clone();
    let one = control
        .advance(&mut whole, target, 0.375, || false)
        .unwrap();
    let mut split = history.clone();
    let middle = core::array::from_fn(|j| BIAS[j] + (target[j] - BIAS[j]) / 3.0);
    let left = control
        .advance(&mut split, middle, 0.125, || false)
        .unwrap();
    let right = control.advance(&mut split, target, 0.25, || false).unwrap();
    for i in 0..6 {
        near(whole.strain()[i], split.strain()[i], 2e-14);
        near(whole.stress()[i], split.stress()[i], 2e-12);
    }
    near(one.work, left.work + right.work, 2e-14);
    near(one.dissipated, left.dissipated + right.dissipated, 2e-14);
    let mut tiny = history.clone();
    let duration = 1e-20;
    let step = control
        .advance(&mut tiny, BIAS, duration, || false)
        .unwrap();
    let expected = duration * ode.heat_rate(&state, ode.rate(&state, plane, [0.0; 6]));
    assert!(step.dissipated > 0.0);
    near(step.dissipated, expected, 5e-11 * expected);
    near(step.stored_energy_change, -expected, 5e-11 * expected);
    assert_eq!(step.work, 0.0);

    // No dynamic modes: exact elastic plane stress plus in-plane Kelvin heat.
    let mut eq = diagonal([2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
    eq[0][2] = 1.0;
    eq[2][0] = 1.0;
    let eta = 0.7;
    let law = TensorProny::new_mandel(eq, vec![], 0.1)
        .unwrap()
        .with_kelvin_viscosity_mandel(diagonal([eta, 0.0, 0.0, 0.0, 0.0, 0.0]))
        .unwrap();
    let mut history = law.history();
    let control = Control::new(
        &history,
        plane,
        TensorKelvinStressControlBudget::default(),
        || false,
    )
    .unwrap();
    assert_eq!(control.transverse_viscosity_rank(), 0);
    let (strain, duration) = (0.002, 0.31);
    let step = control
        .advance(&mut history, [strain, 0.0, 0.0], duration, || false)
        .unwrap();
    near(step.strain[2], -strain / 4.0, 2e-18);
    near(step.strain_rate[2], -strain / (4.0 * duration), 2e-17);
    near(step.dissipated, eta * strain * strain / duration, 2e-18);
    near(step.stored_energy, 0.5 * 1.75 * strain * strain, 2e-18);
    near(step.work, step.dissipated + step.stored_energy, 2e-18);
    let tangent = control.ramp_tangent(duration, || false).unwrap();
    near(tangent.stress[0][0], 1.75 + eta / duration, 2e-14);
    near(tangent.mean_stress[0][0], 0.875 + eta / duration, 2e-14);
}

#[test]
fn ramp_refusals_and_cancellation_never_change_accepted_material_memory() {
    let plane = MaterialPlane::Axes12;
    let (law, _) = fixture(viscosity());
    let mut history = law.history();
    let control = Control::new(
        &history,
        plane,
        TensorKelvinStressControlBudget::default(),
        || false,
    )
    .unwrap();
    let original = history.clone();
    for duration in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            control
                .advance(&mut history, [0.001; 3], duration, || false)
                .is_err()
        );
        assert!(control.ramp_tangent(duration, || false).is_err());
        assert_eq!(history, original);
    }
    assert!(
        control
            .advance(&mut law.history(), [0.001; 3], 0.2, || false)
            .is_err()
    );
    assert!(
        control
            .advance(&mut history, [1.0; 3], 0.2, || false)
            .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        control
            .advance(&mut history, [0.001; 3], 0.2, || {
                polls += 1;
                polls > 3
            })
            .unwrap_err(),
        TensorPronyError::Cancelled
    );
    assert_eq!(
        control.ramp_tangent(0.2, || true).unwrap_err(),
        TensorPronyError::Cancelled
    );
    assert_eq!(history, original);
    let (law, _) = fixture(diagonal([0.0, 0.0, 0.0, 0.6, 0.9, 1.3]));
    let (mut incompatible, _) = seed(&law, plane, BIAS);
    let original = incompatible.clone();
    let control = Control::new(
        &incompatible,
        plane,
        TensorKelvinStressControlBudget::default(),
        || false,
    )
    .unwrap();
    assert!(
        control
            .advance(&mut incompatible, BIAS, 0.2, || false)
            .is_err()
    );
    assert_eq!(incompatible, original);
}
