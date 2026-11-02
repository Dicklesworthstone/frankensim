//! G1/G3/G4: work-conjugate inertia, hereditary heat, and transactional admission.
use super::*;
use fs_material::thermovisco::HistoryIntegration;

fn exact_sources(cp: f64, density: f64, branch: bool, thermal_shift: bool) -> Sources {
    let mut source = sources(cp, density, branch, thermal_shift);
    source.mechanical = source
        .mechanical
        .with_history_integration(HistoryIntegration::ExactStrainRamp);
    source
}

#[test]
fn g1_elastic_midpoint_is_second_order_and_conserves_energy_without_heat() {
    let source = exact_sources(10.0, 1000.0, false, false);
    let g = geometry(0.0);
    let initial =
        Oscillator::try_new(&source.mechanical, &source.caloric, g, 300.0, 0.0, 0.2).unwrap();
    let stiffness = g.cross_section_m2 * 1e6 / g.length_m;
    let omega = (stiffness / initial.effective_mass_kg()).sqrt();
    let end = 1.3 / omega;
    let force = 0.03;
    let equilibrium = force / stiffness;
    let exact_u = equilibrium * (1.0 - (omega * end).cos())
        + initial.velocity_m_s() / omega * (omega * end).sin();
    let exact_v =
        equilibrium * omega * (omega * end).sin() + initial.velocity_m_s() * (omega * end).cos();
    let mut errors = Vec::new();
    for count in [10, 20, 40] {
        let mut state = initial.clone();
        for _ in 0..count {
            let step = state
                .step_midpoint(
                    force,
                    end / f64::from(count),
                    config(),
                    adiabatic,
                    &CancelGate::new(),
                )
                .unwrap();
            assert_eq!(step.material_heat_j, 0.0);
            assert_eq!(step.numerical_loss_j, 0.0);
            assert_eq!(step.inertial_numerical_loss_j, 0.0);
            assert_eq!(
                step.next.point().caloric_state(),
                initial.point().caloric_state()
            );
            state = step.next;
        }
        errors
            .push((state.extension_m() - exact_u).hypot((state.velocity_m_s() - exact_v) / omega));
        close(
            state.kinetic_energy_j()
                + state.point().volume_m3() * state.point().mechanical().stored_energy_j_m3(),
            initial.kinetic_energy_j() + force * state.extension_m(),
            1e-13,
        );
    }
    for adjacent in errors.windows(2) {
        let ratio = adjacent[0] / adjacent[1];
        assert!(
            (3.7..4.3).contains(&ratio),
            "elastic refinement ratio {ratio}: {errors:?}"
        );
    }
    assert!(errors[2] < 1e-7, "elastic finest error: {}", errors[2]);

    let mut state = initial.clone();
    for _ in 0..1000 {
        let step = state
            .step_midpoint(0.0, 0.03 / omega, config(), adiabatic, &CancelGate::new())
            .unwrap();
        assert_eq!(step.material_heat_j, 0.0);
        assert_eq!(step.numerical_loss_j, 0.0);
        assert_eq!(step.inertial_numerical_loss_j, 0.0);
        assert_eq!(
            step.next.point().caloric_state(),
            initial.point().caloric_state()
        );
        state = step.next;
        close(
            state.kinetic_energy_j()
                + state.point().volume_m3() * state.point().mechanical().stored_energy_j_m3(),
            initial.kinetic_energy_j(),
            1e-13,
        );
    }
}

// Independent continuous fixed-coefficient ODE: displacement, velocity, Maxwell
// spring strain, and integrated dashpot heat. No production constitutive update
// or time integrator enters this reference calculation.
fn maxwell_rk4(g: Geometry, mass: f64, force: f64, initial_v: f64, end: f64) -> [f64; 4] {
    let rate = |x: [f64; 4]| {
        let stress = 1e6 * x[0] / g.length_m + 3e6 * x[2];
        [
            x[1],
            (force - g.cross_section_m2 * stress) / mass,
            x[1] / g.length_m - x[2] / 0.005,
            g.cross_section_m2 * g.length_m * 3e6 * x[2] * x[2] / 0.005,
        ]
    };
    let add = |x: [f64; 4], y: [f64; 4], h: f64| core::array::from_fn(|i| x[i] + h * y[i]);
    let count = 20_000;
    let dt = end / f64::from(count);
    let mut x = [0.0, initial_v, 0.0, 0.0];
    for _ in 0..count {
        let k1 = rate(x);
        let k2 = rate(add(x, k1, 0.5 * dt));
        let k3 = rate(add(x, k2, 0.5 * dt));
        let k4 = rate(add(x, k3, dt));
        for i in 0..4 {
            x[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    x
}

#[test]
fn g1_maxwell_motion_memory_and_physical_heat_converge_to_independent_rk4() {
    // The mechanical source is independent of T. Adiabatic caloric storage can
    // therefore record the actual heat without changing the fixed-T ODE oracle.
    let source = exact_sources(1000.0, 1000.0, true, false);
    let g = geometry(0.0);
    let initial =
        Oscillator::try_new(&source.mechanical, &source.caloric, g, 300.0, 0.0, 0.2).unwrap();
    let (end, force) = (0.0073, 0.03);
    let exact = maxwell_rk4(
        g,
        initial.effective_mass_kg(),
        force,
        initial.velocity_m_s(),
        end,
    );
    let omega = (g.cross_section_m2 * 4e6 / (g.length_m * initial.effective_mass_kg())).sqrt();
    let mut errors = Vec::new();
    let mut heat_errors = Vec::new();
    for count in [32, 64, 128] {
        let mut state = initial.clone();
        let mut heat = 0.0;
        for _ in 0..count {
            let step = state
                .step_midpoint(
                    force,
                    end / f64::from(count),
                    config(),
                    adiabatic,
                    &CancelGate::new(),
                )
                .unwrap();
            assert!(step.material_heat_j >= 0.0);
            assert_eq!(step.numerical_loss_j, 0.0);
            assert_eq!(step.inertial_numerical_loss_j, 0.0);
            heat += step.material_heat_j;
            state = step.next;
        }
        let z = state.point().mechanical().branch_strains()[0];
        let error = ((state.extension_m() - exact[0]) / g.length_m)
            .hypot((state.velocity_m_s() - exact[1]) / (g.length_m * omega))
            .hypot(z - exact[2]);
        errors.push(error);
        heat_errors.push((heat - exact[3]).abs());
        let caloric_heat = state.point().mass_kg()
            * (state.point().caloric_state().specific_enthalpy_j_kg()
                - initial.point().caloric_state().specific_enthalpy_j_kg());
        close(caloric_heat, heat, 1e-10);
        close(
            state.kinetic_energy_j()
                + state.point().volume_m3() * state.point().mechanical().stored_energy_j_m3()
                + heat,
            initial.kinetic_energy_j() + force * state.extension_m(),
            1e-11,
        );
        if count == 128 {
            close(state.extension_m(), exact[0], 8e-7);
            close(state.velocity_m_s(), exact[1], 8e-4);
            close(z, exact[2], 8e-6);
            close(heat, exact[3], 3e-8);
        }
    }
    for adjacent in errors.windows(2) {
        let ratio = adjacent[0] / adjacent[1];
        assert!(
            (3.3..4.8).contains(&ratio),
            "Maxwell refinement ratio {ratio}: {errors:?}"
        );
    }
    for adjacent in heat_errors.windows(2) {
        let ratio = adjacent[0] / adjacent[1];
        assert!(
            (3.3..4.8).contains(&ratio),
            "heat refinement ratio {ratio}: {heat_errors:?}"
        );
    }
}

#[test]
fn g3_only_heat_capacity_changes_self_heated_relaxation_and_actual_motion() {
    let low_cp = exact_sources(0.02, 1000.0, true, true);
    let high_cp = exact_sources(2.0, 1000.0, true, true);
    assert_eq!(low_cp.mechanical.law(), high_cp.mechanical.law());
    let mut low = Oscillator::try_new(
        &low_cp.mechanical,
        &low_cp.caloric,
        geometry(0.0),
        300.0,
        0.0,
        0.5,
    )
    .unwrap();
    let mut high = Oscillator::try_new(
        &high_cp.mechanical,
        &high_cp.caloric,
        geometry(0.0),
        300.0,
        0.0,
        0.5,
    )
    .unwrap();
    let initial_energy = low.kinetic_energy_j();
    let initial_low_h = low.point().caloric_state().specific_enthalpy_j_kg();
    let initial_high_h = high.point().caloric_state().specific_enthalpy_j_kg();
    let mut largest_velocity_difference: f64 = 0.0;
    for _ in 0..600 {
        let a = low
            .step_midpoint(0.0, 5e-5, config(), adiabatic, &CancelGate::new())
            .unwrap();
        let b = high
            .step_midpoint(0.0, 5e-5, config(), adiabatic, &CancelGate::new())
            .unwrap();
        assert_eq!(a.numerical_loss_j, 0.0);
        assert_eq!(b.numerical_loss_j, 0.0);
        assert!(a.material_heat_j >= 0.0 && b.material_heat_j >= 0.0);
        low = a.next;
        high = b.next;
        largest_velocity_difference =
            largest_velocity_difference.max((low.velocity_m_s() - high.velocity_m_s()).abs());
    }
    assert!(
        largest_velocity_difference > 1e-6,
        "heat capacity must change the mechanical trajectory"
    );
    assert!(
        low.point().mechanical().temperature_k() > high.point().mechanical().temperature_k() + 0.02
    );
    for (state, initial_h) in [(&low, initial_low_h), (&high, initial_high_h)] {
        close(
            state.point().mechanical().temperature_k(),
            state.point().caloric_state().temperature_k(),
            1e-9,
        );
        close(
            state.kinetic_energy_j()
                + state.point().volume_m3() * state.point().mechanical().stored_energy_j_m3()
                + state.point().mass_kg()
                    * (state.point().caloric_state().specific_enthalpy_j_kg() - initial_h),
            initial_energy,
            1e-8,
        );
    }
}

#[test]
fn g4_failed_midpoint_trials_and_cancellation_preserve_warm_history_and_retry() {
    let source = exact_sources(1000.0, 1000.0, true, false);
    let cold = Oscillator::try_new(
        &source.mechanical,
        &source.caloric,
        geometry(0.0),
        300.0,
        0.0,
        0.2,
    )
    .unwrap();
    let warm = cold
        .step_midpoint(0.03, 1e-4, config(), adiabatic, &CancelGate::new())
        .unwrap()
        .next;
    assert!(warm.point().mechanical().branch_strains()[0] != 0.0);
    let snapshot = warm.clone();
    assert!(
        warm.step_midpoint(1e9, 1e-4, config(), adiabatic, &CancelGate::new())
            .is_err()
    );

    let cancelled = CancelGate::new();
    cancelled.request();
    let calls = core::cell::Cell::new(0);
    assert!(matches!(
        warm.step_midpoint(
            0.03,
            1e-4,
            config(),
            |input| {
                calls.set(calls.get() + 1);
                adiabatic(input)
            },
            &cancelled
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(calls.get(), 0);

    let during_transport = CancelGate::new();
    assert!(matches!(
        warm.step_midpoint(
            0.03,
            1e-4,
            config(),
            |input| {
                let output = adiabatic(input);
                during_transport.request();
                output
            },
            &during_transport
        ),
        Err(Error::Feedback(FeedbackError::Cancelled))
    ));
    assert_eq!(warm.point(), snapshot.point());
    assert_eq!(warm.velocity_m_s(), snapshot.velocity_m_s());
    let retry = warm
        .step_midpoint(-0.01, 1e-4, config(), adiabatic, &CancelGate::new())
        .unwrap();
    let control = snapshot
        .step_midpoint(-0.01, 1e-4, config(), adiabatic, &CancelGate::new())
        .unwrap();
    assert_eq!(retry.next.point(), control.next.point());
    assert_eq!(retry.next.velocity_m_s(), control.next.velocity_m_s());
    assert_eq!(retry.material_heat_j, control.material_heat_j);
    close(retry.next.point().time_s(), 2e-4, 1e-16);
}

#[test]
fn g4_midpoint_checks_endpoint_rate_even_when_mean_strain_rate_is_admitted() {
    let source = exact_sources(1000.0, 1000.0, false, false);
    let g = geometry(0.0);
    let old = Oscillator::try_new(&source.mechanical, &source.caloric, g, 300.0, 0.0, 0.0).unwrap();
    let snapshot = old.clone();
    let dt = 1e-4;
    let rate_limit = 10.0;
    let force_for_mean_rate = |fraction: f64| {
        let strain_change = fraction * rate_limit * dt;
        2.0 * old.effective_mass_kg() * g.length_m * strain_change / (dt * dt)
            + g.cross_section_m2 * 1e6 * strain_change / 2.0
    };
    // Mean strain rate is 7.5/s (admitted), but endpoint velocity is 1.5 m/s,
    // hence its instantaneous strain rate is 15/s, outside the source limit.
    let error = old.step_midpoint(
        force_for_mean_rate(0.75),
        dt,
        config(),
        adiabatic,
        &CancelGate::new(),
    );
    assert!(matches!(
        error,
        Err(Error::Invalid(
            "endpoint velocity exceeds the material strain-rate envelope"
        ))
    ));
    assert_eq!(old.point(), snapshot.point());
    assert_eq!(old.velocity_m_s(), snapshot.velocity_m_s());
    assert!(
        Oscillator::try_new(
            &source.mechanical,
            &source.caloric,
            g,
            300.0,
            0.0,
            1.5 * rate_limit * g.length_m,
        )
        .is_err()
    );
    let admitted = old
        .step_midpoint(
            force_for_mean_rate(0.45),
            dt,
            config(),
            adiabatic,
            &CancelGate::new(),
        )
        .unwrap();
    close(
        admitted.next.velocity_m_s(),
        0.9 * rate_limit * g.length_m,
        1e-12,
    );
}

#[test]
fn g4_backward_euler_source_refuses_midpoint_before_transport() {
    let source = sources(1000.0, 1000.0, true, false);
    let old = Oscillator::try_new(
        &source.mechanical,
        &source.caloric,
        geometry(0.0),
        300.0,
        0.0,
        0.2,
    )
    .unwrap();
    let snapshot = old.clone();
    let calls = core::cell::Cell::new(0);
    let result = old.step_midpoint(
        0.03,
        1e-4,
        config(),
        |input| {
            calls.set(calls.get() + 1);
            adiabatic(input)
        },
        &CancelGate::new(),
    );
    assert!(matches!(result, Err(Error::Invalid(_))));
    assert_eq!(calls.get(), 0);
    assert_eq!(old.point(), snapshot.point());
    assert_eq!(old.velocity_m_s(), snapshot.velocity_m_s());
}
