use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_material::visco_tensor::{
    TemperatureRamp, TemperatureRampControls, TensorProny, TensorPronyBranch, TensorPronyError,
};
use fs_material::{Tangent6, Voigt};

fn shift() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 2.0,
            c2_k: 100.0,
        },
        300.0,
        [250.0, 350.0],
    )
    .unwrap()
}
fn controls(limit: f64) -> TemperatureRampControls {
    TemperatureRampControls {
        maximum_log_shift_change: limit,
        maximum_substeps: 10_000,
    }
}
fn ramp(start_k: f64, end_k: f64, duration_s: f64) -> TemperatureRamp {
    TemperatureRamp {
        start_k,
        end_k,
        duration_s,
    }
}
fn diagonal(values: Voigt) -> Tangent6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn scalar_model(tau: Option<f64>, viscosity: f64) -> TensorProny {
    let branches = tau
        .into_iter()
        .map(|tau| {
            TensorPronyBranch::from_gram_rows(diagonal([2.0, 0.0, 0.0, 0.0, 0.0, 0.0]), tau)
                .unwrap()
        })
        .collect();
    let model = TensorProny::new_mandel(diagonal([2.0; 6]), branches, 0.2).unwrap();
    if viscosity == 0.0 {
        model
    } else {
        model
            .with_kelvin_viscosity_mandel(diagonal([viscosity; 6]))
            .unwrap()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}
fn voigt_work(stress: Voigt, strain: Voigt) -> f64 {
    (0..6)
        .map(|i| {
            if i < 3 {
                stress[i] * strain[i]
            } else {
                2.0 * stress[i] * strain[i]
            }
        })
        .sum()
}

#[test]
fn fixed_temperature_matches_original_shifted_law_and_reference_bits() {
    let calibration = shift();
    let model = coupled_model();
    for temperature in [270.0, 300.0, 340.0] {
        let fixed = calibration.shift_tensor(&model, temperature).unwrap();
        let mut expected = fixed.law().history();
        let mut actual = model.history();
        for strain in [
            [0.02, -0.007, 0.009, 0.004, -0.003, 0.006],
            [-0.01, 0.003, 0.001, -0.002, 0.004, 0.0],
            [0.0; 6],
        ] {
            let reference = expected.advance(strain, 0.23, || false).unwrap();
            let result = actual
                .advance_temperature_ramp(
                    strain,
                    ramp(temperature, temperature, 0.23),
                    &calibration,
                    controls(0.01),
                    || false,
                )
                .unwrap();
            assert_eq!(result.response, reference);
            assert_eq!(result.substeps, 1);
            assert_eq!(result.maximum_log_shift_change, 0.0);
            assert_eq!(actual.stress(), expected.stress());
            assert_eq!(actual.work(), expected.work());
            assert_eq!(actual.dissipated(), expected.dissipated());
            if temperature == 300.0 {
                assert_eq!(actual, expected);
            }
        }
    }
}

#[test]
fn kelvin_endpoint_uses_final_temperature_and_mean_uses_physical_time() {
    let model = scalar_model(None, 0.5);
    let mut history = model.history();
    let calibration = shift();
    let strain = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
    let result = history
        .advance_temperature_ramp(
            strain,
            ramp(300.0, 340.0, 0.2),
            &calibration,
            controls(10.0),
            || false,
        )
        .unwrap();
    assert_eq!(result.substeps, 1);
    let rate = strain[0] / 0.2;
    let midpoint = calibration.evaluate(320.0).unwrap().factor();
    let endpoint = calibration.evaluate(340.0).unwrap().factor();
    let stress = 2.0 * strain[0] + 0.5 * endpoint * rate;
    let mean = strain[0] + 0.5 * midpoint * rate;
    close(result.response.stress[0], stress, 2e-16);
    close(history.stress()[0], stress, 2e-16);
    close(result.response.mean_stress[0], mean, 2e-16);
    assert!((stress - (2.0 * strain[0] + 0.5 * midpoint * rate)).abs() > 1e-3);
    close(result.response.work, mean * strain[0], 2e-17);
    close(
        result.response.dissipated,
        0.5 * midpoint * rate * rate * 0.2,
        2e-17,
    );
    close(result.response.balance_residual, 0.0, 2e-17);
}

#[test]
fn thermal_cycles_retain_loaded_memory_and_never_compound_shifts() {
    let model = scalar_model(Some(0.7), 0.0);
    let calibration = shift();
    let mut history = model.history();
    let strain = [0.02, 0.0, 0.0, 0.0, 0.0, 0.0];
    history.advance(strain, 0.0, || false).unwrap();
    let initial_energy = history.branch_stored_energy();
    let mut reduced_time = 0.0;
    let mut total_heat = 0.0;
    for (start, end) in [(300.0, 340.0), (340.0, 280.0), (280.0, 300.0)] {
        let result = history
            .advance_temperature_ramp(
                strain,
                ramp(start, end, 0.3),
                &calibration,
                controls(0.05),
                || false,
            )
            .unwrap();
        let h = 0.3 / result.substeps as f64;
        for i in 0..result.substeps {
            let temperature = start + (end - start) * ((i as f64 + 0.5) / result.substeps as f64);
            // Independent WLF expression, integrating the declared piecewise
            // midpoint kinetics. Loaded memory is exp(-integral dt/(a_T*tau)).
            let factor =
                10.0_f64.powf(-2.0 * (temperature - 300.0) / (100.0 + temperature - 300.0));
            reduced_time += h / factor;
        }
        let remaining = (-reduced_time / 0.7).exp();
        close(
            history.stress()[0],
            2.0 * strain[0] + 4.0 * strain[0] * remaining,
            2e-14,
        );
        close(
            history.branch_stored_energy(),
            initial_energy * remaining * remaining,
            2e-16,
        );
        close(result.response.work, 0.0, 0.0);
        assert!(result.response.dissipated > 0.0);
        total_heat += result.response.dissipated;
        close(
            total_heat + history.branch_stored_energy(),
            initial_energy,
            3e-16,
        );
    }
    // Resume at the ORIGINAL reference temperature: no history replacement or
    // accidental repeated temperature scaling is involved in this continuation.
    let result = history
        .advance_temperature_ramp(
            strain,
            ramp(300.0, 300.0, 0.2),
            &calibration,
            controls(0.05),
            || false,
        )
        .unwrap();
    reduced_time += 0.2;
    close(
        result.response.stress[0],
        2.0 * strain[0] + 4.0 * strain[0] * (-reduced_time / 0.7).exp(),
        2e-14,
    );
    assert!(history.stress()[0] > 2.0 * strain[0] + 1e-4);
}

fn directions() -> [Tangent6; 2] {
    let mut a = [[0.0; 6]; 6];
    a[0] = [3.0, -1.0, 0.5, 0.6, -0.2, 0.1];
    let mut b = [[0.0; 6]; 6];
    b[0] = [0.0, 2.0, 1.0, -0.5, 0.3, 0.0];
    b[1] = [1.0, 0.0, 2.0, 0.0, 0.4, -0.2];
    [a, b]
}
fn equilibrium() -> Tangent6 {
    diagonal([20.0, 18.0, 16.0, 12.0, 10.0, 8.0])
}
fn viscosity() -> Tangent6 {
    let mut result = diagonal([0.5, 0.4, 0.3, 0.2, 0.15, 0.1]);
    result[0][1] = 0.1;
    result[1][0] = 0.1;
    result
}
fn coupled_model() -> TensorProny {
    let branches = directions()
        .into_iter()
        .zip([0.15, 0.8])
        .map(|(factor, tau)| TensorPronyBranch::from_gram_rows(factor, tau).unwrap())
        .collect();
    TensorProny::new_mandel(equilibrium(), branches, 0.2)
        .unwrap()
        .with_kelvin_viscosity_mandel(viscosity())
        .unwrap()
}
fn multiply(matrix: &Tangent6, vector: &Voigt) -> Voigt {
    core::array::from_fn(|i| (0..6).map(|j| matrix[i][j] * vector[j]).sum())
}
fn mandel_strain(strain: Voigt) -> Voigt {
    core::array::from_fn(|i| {
        if i < 3 {
            strain[i]
        } else {
            core::f64::consts::SQRT_2 * strain[i]
        }
    })
}
fn physical_stress(stress: Voigt) -> Voigt {
    core::array::from_fn(|i| {
        if i < 3 {
            stress[i]
        } else {
            stress[i] / core::f64::consts::SQRT_2
        }
    })
}

struct Reference {
    stress: Voigt,
    stored: f64,
    work: f64,
    heat: f64,
}

// Independent physical equations; no PronyRamp, tensor-history update or
// production temperature evaluation is used. Original q_j coordinates satisfy
// q_j' = B_j epsilon' - q_j/(tau_j*a_T). Work and heat have separate ODE states.
fn reference(strain_end: Voigt, duration: f64, steps: usize) -> Reference {
    let endpoint = mandel_strain(strain_end);
    let rate = endpoint.map(|value| value / duration);
    let factors = directions();
    let tau = [0.15, 0.8];
    let eq = equilibrium();
    let kv = viscosity();
    let observe = |time: f64, state: &[f64; 14]| {
        let strain = rate.map(|value| value * time);
        let temperature = 300.0 + 40.0 * time / duration;
        let factor = 10.0_f64.powf(-2.0 * (temperature - 300.0) / (100.0 + temperature - 300.0));
        let mut stress = multiply(&eq, &strain);
        let viscous = multiply(&kv, &rate);
        let mut heat = (0..6).map(|i| factor * rate[i] * viscous[i]).sum::<f64>();
        let mut stored = 0.5 * (0..6).map(|i| strain[i] * stress[i]).sum::<f64>();
        for branch in 0..2 {
            for i in 0..6 {
                let q = state[6 * branch + i];
                stored += 0.5 * q * q;
                heat += q * q / (tau[branch] * factor);
                for j in 0..6 {
                    stress[j] += factors[branch][i][j] * q;
                }
            }
        }
        for i in 0..6 {
            stress[i] += factor * viscous[i];
        }
        (stress, stored, heat, factor)
    };
    let rhs = |time: f64, state: &[f64; 14]| {
        let (stress, _, heat, factor) = observe(time, state);
        let mut derivative = [0.0; 14];
        for branch in 0..2 {
            let drive = multiply(&factors[branch], &rate);
            for i in 0..6 {
                derivative[6 * branch + i] =
                    drive[i] - state[6 * branch + i] / (tau[branch] * factor);
            }
        }
        derivative[12] = (0..6).map(|i| stress[i] * rate[i]).sum();
        derivative[13] = heat;
        derivative
    };
    let h = duration / steps as f64;
    let mut state = [0.0; 14];
    for step in 0..steps {
        let t = step as f64 * h;
        let k1 = rhs(t, &state);
        let k2 = rhs(
            t + h / 2.0,
            &core::array::from_fn(|i| state[i] + h * k1[i] / 2.0),
        );
        let k3 = rhs(
            t + h / 2.0,
            &core::array::from_fn(|i| state[i] + h * k2[i] / 2.0),
        );
        let k4 = rhs(t + h, &core::array::from_fn(|i| state[i] + h * k3[i]));
        for i in 0..14 {
            state[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    let (stress, stored, _, _) = observe(duration, &state);
    Reference {
        stress: physical_stress(stress),
        stored,
        work: state[12],
        heat: state[13],
    }
}

#[test]
fn noncommuting_tensor_memory_work_and_heat_converge_against_independent_equations() {
    let model = coupled_model();
    let calibration = shift();
    let strain = [0.02, -0.007, 0.009, 0.004, -0.003, 0.006];
    let coarse_reference = reference(strain, 0.9, 4096);
    let truth = reference(strain, 0.9, 8192);
    for i in 0..6 {
        close(coarse_reference.stress[i], truth.stress[i], 3e-13);
    }
    close(coarse_reference.stored, truth.stored, 3e-14);
    close(coarse_reference.work, truth.work, 3e-14);
    close(coarse_reference.heat, truth.heat, 3e-14);
    close(truth.work - truth.stored - truth.heat, 0.0, 3e-14);
    let variation = 40.0
        * calibration
            .evaluate(300.0)
            .unwrap()
            .log_derivative_per_k()
            .abs();
    let mut errors = Vec::new();
    for steps in [16_usize, 32, 64] {
        let controls = TemperatureRampControls {
            maximum_log_shift_change: variation / (steps as f64 - 0.25),
            maximum_substeps: steps,
        };
        let mut history = model.history();
        let result = history
            .advance_temperature_ramp(
                strain,
                ramp(300.0, 340.0, 0.9),
                &calibration,
                controls,
                || false,
            )
            .unwrap();
        assert_eq!(result.substeps, steps);
        assert!(result.maximum_log_shift_change <= controls.maximum_log_shift_change);
        assert_eq!(result.response.stress, history.stress());
        for i in 0..6 {
            close(history.strain()[i], strain[i], 2e-17);
        }
        close(
            result.response.work,
            voigt_work(result.response.mean_stress, strain),
            3e-14,
        );
        close(result.response.balance_residual, 0.0, 3e-14);
        assert!(result.response.dissipated > 0.0);
        let stress_error = (0..6)
            .map(|i| (result.response.stress[i] - truth.stress[i]).abs())
            .fold(0.0_f64, f64::max);
        errors.push([
            stress_error,
            (result.response.stored_energy - truth.stored).abs(),
            (result.response.work - truth.work).abs(),
            (result.response.dissipated - truth.heat).abs(),
        ]);
        let mut replay = model.history();
        assert_eq!(
            result,
            replay
                .advance_temperature_ramp(
                    strain,
                    ramp(300.0, 340.0, 0.9),
                    &calibration,
                    controls,
                    || false,
                )
                .unwrap()
        );
        assert_eq!(history, replay);
    }
    for quantity in 0..4 {
        assert!(
            errors[1][quantity] < errors[0][quantity] / 3.2,
            "{errors:?}"
        );
        assert!(
            errors[2][quantity] < errors[1][quantity] / 3.2,
            "{errors:?}"
        );
    }
    assert!(errors[2][0] < 2e-5);
}

#[test]
fn local_thermal_work_and_heat_survive_large_cumulative_ledgers() {
    let model = scalar_model(None, 1.0);
    let mut history = model.history();
    history
        .advance([0.01, 0.0, 0.0, 0.0, 0.0, 0.0], 1e-100, || false)
        .unwrap();
    history.advance([0.0; 6], 1e-100, || false).unwrap();
    let previous_work = history.work();
    let previous_heat = history.dissipated();
    let result = history
        .advance_temperature_ramp(
            [1e-9, 0.0, 0.0, 0.0, 0.0, 0.0],
            ramp(300.0, 340.0, 1.0),
            &shift(),
            controls(0.2),
            || false,
        )
        .unwrap();
    // The cumulative f64 ledgers cannot resolve this increment, but the local
    // physical quadrature still returns it rather than subtracting totals.
    assert_eq!(history.work(), previous_work);
    assert_eq!(history.dissipated(), previous_heat);
    assert!(result.response.work > 1e-18);
    assert!(result.response.dissipated > 1e-19);
    close(
        result.response.work - result.response.stored_energy - result.response.dissipated,
        0.0,
        2e-32,
    );
}

#[test]
fn thermal_budget_calibration_and_final_cancellation_refuse_atomically() {
    let model = coupled_model();
    let calibration = shift();
    let mut history = model.history();
    history
        .advance([0.01, 0.0, 0.0, 0.0, 0.0, 0.0], 0.1, || false)
        .unwrap();
    let before = history.clone();
    let strain = [0.02, -0.007, 0.009, 0.004, -0.003, 0.006];
    for thermal in [
        ramp(300.0, 360.0, 0.9),
        ramp(300.0, 340.0, 0.0),
        ramp(300.0, 340.0, f64::from_bits(1)),
        ramp(300.0, 340.0, f64::NAN),
    ] {
        assert!(matches!(
            history
                .advance_temperature_ramp(strain, thermal, &calibration, controls(0.1), || false,),
            Err(TensorPronyError::Invalid(_))
        ));
        assert_eq!(history, before);
    }
    let limited = TemperatureRampControls {
        maximum_log_shift_change: 0.01,
        maximum_substeps: 1,
    };
    assert!(matches!(
        history.advance_temperature_ramp(
            strain,
            ramp(300.0, 340.0, 0.9),
            &calibration,
            limited,
            || false,
        ),
        Err(TensorPronyError::Invalid(
            "temperature ramp substep budget exhausted"
        ))
    ));
    assert_eq!(history, before);
    assert!(
        history
            .advance_temperature_ramp(
                [0.5; 6],
                ramp(300.0, 340.0, 0.9),
                &calibration,
                controls(0.1),
                || false,
            )
            .is_err()
    );
    assert_eq!(history, before);
    let mut successful = before.clone();
    let mut final_poll = 0;
    successful
        .advance_temperature_ramp(
            strain,
            ramp(300.0, 340.0, 0.9),
            &calibration,
            controls(0.1),
            || {
                final_poll += 1;
                false
            },
        )
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        history.advance_temperature_ramp(
            strain,
            ramp(300.0, 340.0, 0.9),
            &calibration,
            controls(0.1),
            || {
                polls += 1;
                polls == final_poll
            },
        ),
        Err(TensorPronyError::Cancelled)
    ));
    assert_eq!(history, before);

    let mut slow = scalar_model(Some(f64::MAX / 4.0), 0.0).history();
    let before = slow.clone();
    assert!(matches!(
        slow.advance_temperature_ramp(
            [0.0; 6],
            ramp(300.0, 250.0, 1.0),
            &calibration,
            controls(0.1),
            || false,
        ),
        Err(TensorPronyError::Invalid(_))
    ));
    assert_eq!(slow, before);
}
