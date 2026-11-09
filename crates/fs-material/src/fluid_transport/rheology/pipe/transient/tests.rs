use super::super::CarreauPipeControls;
use super::*;
use crate::fluid_transport::rheology::{CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode};
use crate::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}
fn law(mu0: f64, thinning: bool, mode: ViscosityShiftMode) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            mu0,
            if thinning { 0.01 } else { mu0 },
            0.8,
            0.5,
            1.0,
            [0.0, 1e7],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 9000.0,
            },
            300.0,
            [280.0, 350.0],
        )
        .unwrap(),
        mode,
    )
}
fn model(cells: usize, thinning: bool) -> RadialCarreauPipe {
    RadialCarreauPipe::new(
        LaminarCarreauPipe::new(
            law(0.2, thinning, ViscosityShiftMode::Vertical),
            300.0,
            0.02,
            0.7,
            1e5,
        )
        .unwrap(),
        1000.0,
        cells,
    )
    .unwrap()
}
fn controls() -> RadialCarreauControls {
    RadialCarreauControls {
        absolute_velocity_tolerance_m_s: 1e-13,
        relative_velocity_tolerance: 1e-12,
        ..RadialCarreauControls::default()
    }
}
fn balance(step: &RadialCarreauStep) {
    assert!(step.viscous_heat_j >= 0.0);
    let scale =
        step.pressure_work_j.abs() + step.kinetic_energy_change_j.abs() + step.viscous_heat_j;
    assert!(
        step.work_residual_j.abs() <= 3e-9 * scale + 1e-18,
        "{step:?}"
    );
    near(
        step.work_residual_j,
        step.momentum_residual_work_j,
        1e-4,
        5e-15 * scale + 1e-20,
    );
    assert!(step.midpoint_pressure_conductance_m3_s_pa > 0.0);
}

#[test]
fn g1_prescribed_temperature_continues_the_exact_discrete_diffusion_mode() {
    let original = model(2, false);
    let mut state = original
        .state_from_velocities(vec![1.0, 1.0 / 3.0], 0.4, || false)
        .unwrap();
    let mut work = original.workspace().unwrap();
    let mut amplitude = 1.0;
    let mut previous_energy = None;
    for temperature in [300.0, 330.0, 290.0, 315.0] {
        let before = state.velocities_m_s().to_vec();
        let time = state.accepted_time_s();
        let pipe = original.at_temperature(temperature).unwrap();
        assert_eq!(state.velocities_m_s(), before);
        assert_eq!(state.accepted_time_s(), time);
        assert_eq!(pipe.pipe().temperature_k(), temperature);
        assert_eq!(original.pipe().temperature_k(), 300.0);
        assert_eq!(pipe.density_kg_m3(), original.density_kg_m3());
        assert_eq!(pipe.pipe().radius_m(), original.pipe().radius_m());
        assert_eq!(pipe.pipe().length_m(), original.pipe().length_m());
        // Independent Arrhenius coefficient and exact eigenvector [1, 1/3]
        // of the two-ring finite-volume operator, including each held change.
        let shift = (9000.0 / 8.314_462_618_153_24 * (1.0 / temperature - 1.0 / 300.0)).exp();
        near(
            pipe.radial_diffusion_time_s(),
            original.radial_diffusion_time_s() / shift,
            3e-14,
            0.0,
        );
        let lambda = 4.0 * shift / (3.0 * original.radial_diffusion_time_s());
        let dt = 0.02;
        amplitude *= (1.0 - 0.5 * lambda * dt) / (1.0 + 0.5 * lambda * dt);
        let step = pipe
            .step(&mut state, &mut work, 0.0, dt, controls(), || false)
            .unwrap();
        near(state.velocities_m_s()[0], amplitude, 3e-14, 1e-15);
        near(state.velocities_m_s()[1], amplitude / 3.0, 3e-14, 1e-15);
        if let Some(energy) = previous_energy {
            assert_eq!(step.initial_kinetic_energy_j, energy);
        }
        previous_energy = Some(step.final_kinetic_energy_j);
        near(
            step.viscous_heat_j,
            -step.kinetic_energy_change_j,
            3e-14,
            1e-16,
        );
        balance(&step);
    }
}

#[test]
fn g3_thermal_continuation_reuses_the_full_nonlinear_profile_and_workspace() {
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        let source = law(0.2, true, mode);
        let original = RadialCarreauPipe::new(
            LaminarCarreauPipe::new(source, 300.0, 0.02, 0.7, 1e5).unwrap(),
            1000.0,
            24,
        )
        .unwrap();
        let mut state = original.rest_state().unwrap();
        let mut work = original.workspace().unwrap();
        original
            .step(&mut state, &mut work, 70.0, 0.02, controls(), || false)
            .unwrap();
        let mut cold = state.clone();
        let changed = original.at_temperature(330.0).unwrap();
        let fresh = RadialCarreauPipe::new(
            LaminarCarreauPipe::new(source, 330.0, 0.02, 0.7, 1e5).unwrap(),
            1000.0,
            24,
        )
        .unwrap();
        let mut reference = fresh
            .state_from_velocities(
                state.velocities_m_s().to_vec(),
                state.accepted_time_s(),
                || false,
            )
            .unwrap();
        let mut fresh_work = fresh.workspace().unwrap();
        let cold_step = original
            .step(&mut cold, &mut work, 0.0, 0.0005, controls(), || false)
            .unwrap();
        let step = changed
            .step(&mut state, &mut work, 0.0, 0.0005, controls(), || false)
            .unwrap();
        let expected = fresh
            .step(
                &mut reference,
                &mut fresh_work,
                0.0,
                0.0005,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(step, expected);
        assert_eq!(state.velocities_m_s(), reference.velocities_m_s());
        assert_eq!(state.accepted_time_s(), reference.accepted_time_s());
        assert_eq!(
            step.initial_kinetic_energy_j,
            cold_step.initial_kinetic_energy_j
        );
        assert!(step.viscous_heat_j < cold_step.viscous_heat_j);
        assert!(step.final_kinetic_energy_j > cold_step.final_kinetic_energy_j);
        balance(&step);
    }
}

#[test]
fn g0_thermal_continuation_preserves_the_original_source_and_physical_rate_ceiling() {
    let source = ThermalCarreauYasuda::new(
        CarreauYasuda::new(0.2, 0.01, 0.8, 0.5, 1.0, [0.0, 100.0]).unwrap(),
        *law(0.2, true, ViscosityShiftMode::VerticalAndRate).temperature_shift(),
        ViscosityShiftMode::VerticalAndRate,
    );
    let original = RadialCarreauPipe::new(
        LaminarCarreauPipe::new(source, 300.0, 0.02, 0.7, 100.0).unwrap(),
        1000.0,
        8,
    )
    .unwrap();
    let state = original.rest_state().unwrap();
    for temperature in [0.0, 279.0, 351.0, f64::NAN, f64::INFINITY, 290.0] {
        // 290 K is within the thermal range, but its a_T > 1 moves the
        // unchanged physical ceiling beyond the reference-rate support.
        assert!(matches!(
            original.at_temperature(temperature),
            Err(RadialCarreauError::Material(_))
        ));
    }
    let hot = original.at_temperature(330.0).unwrap();
    assert_eq!(hot.pipe().maximum_laminar_wall_rate_per_s(), 100.0);
    assert_eq!(hot.pipe().law(), original.pipe().law());
    assert_eq!(state.velocities_m_s(), &[0.0; 8]);
    assert_eq!(state.accepted_time_s(), 0.0);
    assert_eq!(original.pipe().temperature_k(), 300.0);
}

#[test]
fn immutable_pressure_trials_share_the_original_solve_and_publish_only_complete_outputs() {
    let pipe = model(24, true);
    let mut initial = pipe.rest_state().unwrap();
    let mut work = pipe.workspace().unwrap();
    pipe.step(&mut initial, &mut work, 40.0, 0.02, controls(), || false)
        .unwrap();
    let initial_bits: Vec<_> = initial
        .velocities_m_s()
        .iter()
        .map(|x| x.to_bits())
        .collect();
    let initial_time = initial.accepted_time_s();
    let mut output = pipe.rest_state().unwrap();
    let mut polls = 0;
    let accepted = pipe
        .step_from(
            &initial,
            &mut output,
            &mut work,
            -12.0,
            0.03,
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut ordinary = initial.clone();
    let ordinary_report = pipe
        .step(&mut ordinary, &mut work, -12.0, 0.03, controls(), || false)
        .unwrap();
    assert_eq!(accepted, ordinary_report);
    assert_eq!(output.velocities_m_s(), ordinary.velocities_m_s());
    assert_eq!(output.accepted_time_s(), ordinary.accepted_time_s());
    assert_eq!(initial.accepted_time_s(), initial_time);
    assert_eq!(
        initial
            .velocities_m_s()
            .iter()
            .map(|x| x.to_bits())
            .collect::<Vec<_>>(),
        initial_bits
    );
    for stop in [1, polls / 2, polls] {
        let mut rejected = output.clone();
        let before = rejected.velocities_m_s().to_vec();
        let before_time = rejected.accepted_time_s();
        let mut calls = 0;
        assert_eq!(
            pipe.step_from(
                &initial,
                &mut rejected,
                &mut work,
                -12.0,
                0.03,
                controls(),
                || {
                    calls += 1;
                    calls == stop
                }
            ),
            Err(RadialCarreauError::Interrupted)
        );
        assert_eq!(rejected.velocities_m_s(), before);
        assert_eq!(rejected.accepted_time_s(), before_time);
        assert_eq!(
            pipe.step_from(
                &initial,
                &mut rejected,
                &mut work,
                -12.0,
                0.03,
                controls(),
                || false
            )
            .unwrap(),
            accepted
        );
        assert_eq!(rejected.velocities_m_s(), output.velocities_m_s());
    }
    let limited = RadialCarreauControls {
        maximum_material_evaluations: 1,
        ..controls()
    };
    let before = output.velocities_m_s().to_vec();
    let time = output.accepted_time_s();
    assert_eq!(
        pipe.step_from(
            &initial,
            &mut output,
            &mut work,
            20.0,
            0.03,
            limited,
            || false
        ),
        Err(RadialCarreauError::MaterialBudget)
    );
    assert_eq!(output.velocities_m_s(), before);
    assert_eq!(output.accepted_time_s(), time);
    let other = model(24, true);
    let mut foreign = other.rest_state().unwrap();
    assert_eq!(
        pipe.step_from(
            &initial,
            &mut foreign,
            &mut work,
            20.0,
            0.03,
            controls(),
            || false
        ),
        Err(RadialCarreauError::StateMismatch)
    );
}

// Independent convergent Bessel series; the production radial discretization
// uses neither Bessel functions nor pre-fitted exponential memory arms.
fn j0(x: f64) -> f64 {
    let mut term = 1.0;
    let mut sum = 1.0;
    for k in 1..=64 {
        term *= -0.25 * x * x / (k * k) as f64;
        sum += term;
        if term.abs() < 1e-19 {
            break;
        }
    }
    sum
}

#[test]
fn newtonian_startup_and_bessel_decay_converge_with_radial_resolution() {
    // Exact pressure-step mean flow: Q/Qss = 1-32 sum(exp(-j0k²*tau)/j0k⁴).
    // At tau=.1 the omitted modes beyond these eight are below1e-25.
    let roots: [f64; 8] = [
        2.404825557695773,
        5.520078110286311,
        8.653727912911013,
        11.791534439014281,
        14.930917708487787,
        18.071063967910922,
        21.211636629879258,
        24.352471530749302,
    ];
    let pressure = 70.0;
    let radius: f64 = 0.02;
    let mu = 0.2;
    let duration = 0.1 * 1000.0 * radius * radius / mu;
    let exact_flow = PI * radius.powi(4) * pressure / (8.0 * mu * 0.7)
        * (1.0
            - 32.0
                * roots
                    .iter()
                    .map(|a| (-0.1 * a * a).exp() / a.powi(4))
                    .sum::<f64>());
    let mut flow_errors = Vec::new();
    let mut field_errors = Vec::new();
    for cells in [16, 32, 64] {
        let pipe = model(cells, false);
        let mut work = pipe.workspace().unwrap();
        let mut state = pipe.rest_state().unwrap();
        let mut last = None;
        for _ in 0..512 {
            let report = pipe
                .step(
                    &mut state,
                    &mut work,
                    pressure,
                    duration / 512.0,
                    controls(),
                    || false,
                )
                .unwrap();
            balance(&report);
            last = Some(report);
        }
        flow_errors.push((last.unwrap().final_volume_flow_m3_s / exact_flow - 1.0).abs());
        let start: Vec<f64> = (0..cells)
            .map(|i| j0(roots[0] * pipe.cell_radius_m(i).unwrap() / radius))
            .collect();
        let mut state = pipe.state_from_velocities(start, 0.0, || false).unwrap();
        for _ in 0..512 {
            balance(
                &pipe
                    .step(
                        &mut state,
                        &mut work,
                        0.0,
                        duration / 512.0,
                        controls(),
                        || false,
                    )
                    .unwrap(),
            );
        }
        let mut error = 0.0_f64;
        for (i, &velocity) in state.velocities_m_s().iter().enumerate() {
            let exact = j0(roots[0] * pipe.cell_radius_m(i).unwrap() / radius)
                * (-0.1 * roots[0].powi(2)).exp();
            error = error.max((velocity - exact).abs());
        }
        field_errors.push(error);
    }
    for errors in [flow_errors, field_errors] {
        assert!(
            errors[0] / errors[1] > 3.5 && errors[1] / errors[2] > 3.5,
            "{errors:?}"
        );
        assert!(errors[2] < 6e-4, "{errors:?}");
    }
}

#[test]
fn midpoint_has_second_order_time_convergence_and_exact_discrete_mode_energy() {
    // For two annuli, M^-1 K = (nu/dr²)*[[2,-2],[-2/3,10/3]].
    // [1,1/3] is its slow eigenvector, lambda=4nu/(3dr²).
    let pipe = model(2, false);
    let lambda = 4.0 / (3.0 * pipe.radial_diffusion_time_s());
    let duration = 1.0 / lambda;
    let exact = (-1.0_f64).exp();
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let mut state = pipe
            .state_from_velocities(vec![1.0, 1.0 / 3.0], 0.0, || false)
            .unwrap();
        let mut work = pipe.workspace().unwrap();
        let dt = duration / steps as f64;
        let factor = (1.0 - 0.5 * lambda * dt) / (1.0 + 0.5 * lambda * dt);
        for step in 1..=steps {
            let report = pipe
                .step(&mut state, &mut work, 0.0, dt, controls(), || false)
                .unwrap();
            near(state.velocities_m_s()[0], factor.powi(step), 3e-14, 1e-15);
            near(
                state.velocities_m_s()[1],
                factor.powi(step) / 3.0,
                3e-14,
                1e-15,
            );
            near(
                report.viscous_heat_j,
                -report.kinetic_energy_change_j,
                3e-14,
                1e-16,
            );
            balance(&report);
        }
        errors.push((state.velocities_m_s()[0] - exact).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9,
        "{errors:?}"
    );
}

#[test]
fn nonlinear_steady_profile_limit_and_discrete_stationarity_retain_original_law() {
    let mut errors = Vec::new();
    for cells in [16, 32, 64] {
        let pipe = model(cells, true);
        let original = pipe
            .pipe()
            .profile_at_pressure(
                250.0,
                CarreauPipeControls {
                    relative_quadrature_tolerance: 1e-10,
                    ..CarreauPipeControls::default()
                },
            )
            .unwrap();
        // Integrate the independently available original steady radial gradient
        // across the finite-volume face distances. The actual transient owner
        // must preserve this stationary discrete flux exactly.
        let mut velocity = vec![0.0; cells];
        let wall = original.sample(pipe.pipe().radius_m()).unwrap();
        velocity[cells - 1] = -0.5 * pipe.spacing_m() * wall.radial_velocity_gradient_per_s;
        for i in (0..cells - 1).rev() {
            let sample = original.sample((i + 1) as f64 * pipe.spacing_m()).unwrap();
            velocity[i] =
                velocity[i + 1] - pipe.spacing_m() * sample.radial_velocity_gradient_per_s;
        }
        let before = velocity.clone();
        let mut state = pipe.state_from_velocities(velocity, 0.0, || false).unwrap();
        let mut work = pipe.workspace().unwrap();
        let report = pipe
            .step(&mut state, &mut work, 250.0, 0.01, controls(), || false)
            .unwrap();
        for (&a, &b) in state.velocities_m_s().iter().zip(&before) {
            near(a, b, 2e-9, 1e-13);
        }
        near(report.pressure_work_j, report.viscous_heat_j, 1e-9, 1e-15);
        balance(&report);
        errors.push(
            (report.final_volume_flow_m3_s / original.flow().volumetric_flow_m3_s - 1.0).abs(),
        );
        assert_eq!(*pipe.pipe().law(), *original.pipe().law());
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
    assert!(errors[2] < 9e-4, "{errors:?}");
}

#[test]
fn nonlinear_startup_reversal_and_tangent_use_actual_radial_history() {
    let pipe = model(40, true);
    let mut state = pipe.rest_state().unwrap();
    let mut work = pipe.workspace().unwrap();
    for _ in 0..5 {
        balance(
            &pipe
                .step(&mut state, &mut work, 150.0, 0.01, controls(), || false)
                .unwrap(),
        );
    }
    let initial = state.clone();
    let report = pipe
        .step(&mut state, &mut work, -150.0, 0.01, controls(), || false)
        .unwrap();
    balance(&report);
    assert!(report.initial_volume_flow_m3_s > 0.0);
    assert!(
        report.final_volume_flow_m3_s > 0.0,
        "pressure reversal does not instantly erase radial momentum"
    );
    assert!(report.pressure_work_j < 0.0);
    assert!(report.kinetic_energy_change_j < 0.0);
    for _ in 0..12 {
        balance(
            &pipe
                .step(&mut state, &mut work, -150.0, 0.01, controls(), || false)
                .unwrap(),
        );
    }
    assert!(state.velocities_m_s().iter().all(|u| *u < 0.0));
    let mut plus = initial.clone();
    let mut minus = initial.clone();
    let perturbation = 1e-3;
    let upper = pipe
        .step(
            &mut plus,
            &mut work,
            -150.0 + perturbation,
            0.01,
            controls(),
            || false,
        )
        .unwrap();
    let lower = pipe
        .step(
            &mut minus,
            &mut work,
            -150.0 - perturbation,
            0.01,
            controls(),
            || false,
        )
        .unwrap();
    near(
        report.midpoint_pressure_conductance_m3_s_pa,
        (upper.midpoint_volume_flow_m3_s - lower.midpoint_volume_flow_m3_s) / (2.0 * perturbation),
        2e-6,
        0.0,
    );
    let reflected: Vec<_> = initial.velocities_m_s().iter().map(|u| -*u).collect();
    let mut mirror = pipe
        .state_from_velocities(reflected, initial.accepted_time_s(), || false)
        .unwrap();
    let mirrored = pipe
        .step(&mut mirror, &mut work, 150.0, 0.01, controls(), || false)
        .unwrap();
    near(
        mirrored.midpoint_volume_flow_m3_s,
        -report.midpoint_volume_flow_m3_s,
        2e-13,
        1e-20,
    );
    near(mirrored.viscous_heat_j, report.viscous_heat_j, 2e-13, 1e-20);
}

#[test]
fn equal_mean_flow_does_not_erase_radial_memory() {
    let pipe = model(32, false);
    let shape = |power: i32| -> Vec<f64> {
        (0..pipe.radial_cells())
            .map(|i| {
                0.1 * (1.0 - (pipe.cell_radius_m(i).unwrap() / pipe.pipe().radius_m()).powi(power))
            })
            .collect()
    };
    let first = shape(2);
    let mut second = shape(4);
    let sum = |u: &[f64]| -> f64 {
        u.iter()
            .enumerate()
            .map(|(i, v)| (2 * i + 1) as f64 * v)
            .sum()
    };
    let ratio = sum(&first) / sum(&second);
    for u in &mut second {
        *u *= ratio;
    }
    let mut a = pipe.state_from_velocities(first, 0.0, || false).unwrap();
    let mut b = pipe.state_from_velocities(second, 0.0, || false).unwrap();
    let mut work = pipe.workspace().unwrap();
    let ra = pipe
        .step(&mut a, &mut work, 0.0, 0.01, controls(), || false)
        .unwrap();
    let rb = pipe
        .step(&mut b, &mut work, 0.0, 0.01, controls(), || false)
        .unwrap();
    near(
        ra.initial_volume_flow_m3_s,
        rb.initial_volume_flow_m3_s,
        1e-14,
        0.0,
    );
    assert!(
        (ra.final_volume_flow_m3_s - rb.final_volume_flow_m3_s).abs()
            > 0.005 * ra.initial_volume_flow_m3_s
    );
    assert!((ra.viscous_heat_j - rb.viscous_heat_j).abs() > 0.01 * ra.viscous_heat_j);
    balance(&ra);
    balance(&rb);
}

#[test]
fn thermal_modes_reuse_the_original_shift_and_support() {
    let temperature = 325.0;
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        let source = law(0.2, true, mode);
        let factor = source
            .temperature_shift()
            .evaluate(temperature)
            .unwrap()
            .factor();
        let equivalent = ThermalCarreauYasuda::new(
            CarreauYasuda::new(
                factor * 0.2,
                factor * 0.01,
                if mode == ViscosityShiftMode::Vertical {
                    0.8
                } else {
                    0.8 * factor
                },
                0.5,
                1.0,
                [0.0, 1e7],
            )
            .unwrap(),
            *source.temperature_shift(),
            ViscosityShiftMode::Vertical,
        );
        let make = |law, temp| {
            RadialCarreauPipe::new(
                LaminarCarreauPipe::new(law, temp, 0.02, 0.7, 1e5).unwrap(),
                1000.0,
                24,
            )
            .unwrap()
        };
        let shifted = make(source, temperature);
        let reference = make(equivalent, 300.0);
        let mut a = shifted.rest_state().unwrap();
        let mut b = reference.rest_state().unwrap();
        let mut wa = shifted.workspace().unwrap();
        let mut wb = reference.workspace().unwrap();
        let ra = shifted
            .step(&mut a, &mut wa, 100.0, 0.02, controls(), || false)
            .unwrap();
        let rb = reference
            .step(&mut b, &mut wb, 100.0, 0.02, controls(), || false)
            .unwrap();
        near(
            ra.final_volume_flow_m3_s,
            rb.final_volume_flow_m3_s,
            2e-12,
            1e-18,
        );
        near(ra.viscous_heat_j, rb.viscous_heat_j, 2e-12, 1e-18);
        assert_eq!(shifted.pipe().temperature_k(), temperature);
        assert_eq!(*shifted.pipe().law(), source);
    }
}

#[test]
fn source_work_budgets_and_late_cancellation_preserve_state_and_retry() {
    let pipe = model(24, true);
    let mut state = pipe.rest_state().unwrap();
    let initial = state.clone();
    let mut work = pipe.workspace().unwrap();
    assert_eq!(
        pipe.step(
            &mut state,
            &mut work,
            120.0,
            0.03,
            RadialCarreauControls {
                maximum_material_evaluations: 1,
                ..controls()
            },
            || false
        ),
        Err(RadialCarreauError::MaterialBudget)
    );
    assert_eq!(state.velocities_m_s(), initial.velocities_m_s());
    assert_eq!(
        state.accepted_time_s().to_bits(),
        initial.accepted_time_s().to_bits()
    );
    assert_eq!(
        pipe.step(
            &mut state,
            &mut work,
            120.0,
            0.03,
            RadialCarreauControls {
                maximum_newton_iterations: 1,
                ..controls()
            },
            || false
        ),
        Err(RadialCarreauError::NewtonBudget)
    );
    let mut checks = 0;
    let expected = pipe
        .step(&mut state, &mut work, 120.0, 0.03, controls(), || {
            checks += 1;
            false
        })
        .unwrap();
    let accepted = state.clone();
    let mut retried = initial.clone();
    let mut calls = 0;
    assert_eq!(
        pipe.step(&mut retried, &mut work, 120.0, 0.03, controls(), || {
            calls += 1;
            calls == checks
        }),
        Err(RadialCarreauError::Interrupted)
    );
    assert_eq!(retried.velocities_m_s(), initial.velocities_m_s());
    assert_eq!(
        retried.accepted_time_s().to_bits(),
        initial.accepted_time_s().to_bits()
    );
    let actual = pipe
        .step(&mut retried, &mut work, 120.0, 0.03, controls(), || false)
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(retried.velocities_m_s(), accepted.velocities_m_s());
    assert_eq!(
        retried.accepted_time_s().to_bits(),
        accepted.accepted_time_s().to_bits()
    );
    let other = model(24, true);
    assert_eq!(
        other.step(&mut state, &mut work, 1.0, 0.01, controls(), || false),
        Err(RadialCarreauError::StateMismatch)
    );
    assert!(RadialCarreauPipe::new(*pipe.pipe(), 0.0, 24).is_err());
    assert!(RadialCarreauPipe::new(*pipe.pipe(), 1000.0, 1).is_err());
    let capped = RadialCarreauPipe::new(
        LaminarCarreauPipe::new(*pipe.pipe().law(), 300.0, 0.02, 0.7, 1.0).unwrap(),
        1000.0,
        24,
    )
    .unwrap();
    assert!(
        capped
            .state_from_velocities(vec![1.0; 24], 0.0, || false)
            .is_err()
    );
    let mut capped_state = capped.rest_state().unwrap();
    let mut capped_work = capped.workspace().unwrap();
    assert!(
        capped
            .step(
                &mut capped_state,
                &mut capped_work,
                1e6,
                0.1,
                controls(),
                || false
            )
            .is_err()
    );
    assert!(capped_state.velocities_m_s().iter().all(|u| *u == 0.0));
    assert_eq!(capped_state.accepted_time_s(), 0.0);
}
