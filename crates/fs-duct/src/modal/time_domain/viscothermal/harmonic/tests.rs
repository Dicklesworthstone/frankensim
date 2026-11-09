//! G1/G2/G3: exact wall-state elimination and independent pressure/work checks.
use super::*;
use fs_material::gas::GasSpec;
use fs_phs::harmonic::HarmonicModel;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 124,
        max_dense_bytes: 32 * 1024 * 1024,
        max_dense_work: 200_000_000,
        relative_residual_tolerance: 2e-10,
    }
}
fn resolution(shear: f64, poles: usize) -> ModalViscothermalResolution {
    let gas = gas();
    let omega = shear * shear * gas.dynamic_viscosity / (gas.density * 0.002_f64.powi(2));
    ModalViscothermalResolution {
        poles,
        band_rad_s: [0.99 * omega, 1.01 * omega],
        maximum_diffusive_relative_error: 0.5,
        maximum_boundary_layer_ratio: 0.1,
    }
}
fn source(elements: usize, modes: usize, length: f64, end: ModalEnd) -> ModalGalerkin {
    assemble_lossless(
        &RadiusProfile::exponential(0.002, 0.0, length).unwrap(),
        &gas(),
        ModeCount(modes),
        elements,
        end,
        124,
        || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}
fn vector_close(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected
        .iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max)
        .max(1e-30);
    for (i, (a, b)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (*a - *b).abs() <= tolerance * scale,
            "coordinate {i}: {a:?} != {b:?}; scale {scale:e}"
        );
    }
}

#[test]
fn condensed_response_recovers_every_original_wall_state_and_observable() {
    // At shear 500 both retained radial modes propagate; the lower-frequency
    // case also verifies evanescent coordinates and their physical wall work.
    for (shear, poles) in [(100.0, 16), (500.0, 38)] {
        for end in [ModalEnd::Rigid, ModalEnd::PressureRelease, ModalEnd::Ports] {
            let resolution = resolution(shear, poles);
            let omega = resolution.band_rad_s[0] / 0.99;
            let full = source(1, 2, 0.006, end)
                .with_viscothermal(
                    resolution,
                    ModalViscothermalBudget {
                        maximum_states: 124,
                        maximum_payload_bytes: 32 * 1024 * 1024,
                        maximum_preparation_work: 200_000_000,
                    },
                    || false,
                )
                .unwrap();
            let condensed = source(1, 2, 0.006, end)
                .into_viscothermal_harmonic(resolution, budget(), || false)
                .unwrap();
            let reference = HarmonicModel::from_ode(full.system(), budget(), || false).unwrap();
            let forcing: Vec<_> = (0..full.source().port_count())
                .map(|i| C64::new((i + 1) as f64 * 1e-7, (i % 2) as f64 * -2e-7))
                .collect();
            let exact = reference
                .factor(omega, || false)
                .unwrap()
                .response(&forcing, || false)
                .unwrap();
            let actual = condensed
                .factor(omega, || false)
                .unwrap()
                .response(&forcing, || false)
                .unwrap();
            let n = full.source().coordinates();
            for range in [0..n, n..2 * n, 2 * n..exact.state.len()] {
                vector_close(&actual.state[range.clone()], &exact.state[range], 2e-8);
            }
            vector_close(&actual.port_pressure_pa, &exact.output, 2e-8);
            close(actual.stored_energy_j, exact.mean_stored_energy_j, 2e-8);
            close(actual.supplied_power_w, exact.supplied_power_w, 2e-8);
            close(
                actual.heat.viscous_w + actual.heat.thermal_w,
                exact.dissipated_power_w,
                2e-8,
            );
            assert!(actual.relative_residual < budget().relative_residual_tolerance);
        }
    }
}

#[test]
fn condensed_wall_load_is_reciprocal_with_live_signed_cross_mode_coupling() {
    let resolution = resolution(200.0, 32);
    let omega = resolution.band_rad_s[0] / 0.99;
    let model = source(6, 3, 0.002, ModalEnd::Ports)
        .into_viscothermal_harmonic(resolution, budget(), || false)
        .unwrap();
    let factor = model.factor(omega, || false).unwrap();
    let ports = model.source().port_count();
    let mut impedance = vec![C64::ZERO; ports * ports];
    for column in 0..ports {
        let mut flow = vec![C64::ZERO; ports];
        flow[column] = C64::from_re(1e-7);
        let response = factor.response(&flow, || false).unwrap();
        assert!(response.heat.viscous_w >= 0.0 && response.heat.thermal_w >= 0.0);
        close(
            response.supplied_power_w,
            response.heat.viscous_w + response.heat.thermal_w,
            2e-8,
        );
        for row in 0..ports {
            impedance[row * ports + column] = response.port_pressure_pa[row].scale(1e7);
        }
    }
    let scale = model.source().inlet_reference_impedance();
    for row in 0..ports {
        for column in 0..ports {
            assert!(
                (impedance[row * ports + column] - impedance[column * ports + row]).abs()
                    < 2e-8 * scale
            );
        }
    }
    assert!(impedance[1].abs() > 1e-5 * scale);
    for sign in [1.0, -1.0] {
        let mut flow = vec![C64::ZERO; ports];
        flow[0] = C64::from_re(1e-7);
        flow[1] = C64::from_re(sign * 1e-7);
        let response = factor.response(&flow, || false).unwrap();
        close(
            response.supplied_power_w,
            response.heat.viscous_w + response.heat.thermal_w,
            2e-8,
        );
    }
}

#[test]
fn all_wall_histories_survive_beyond_the_complete_harmonic_dimension_cap() {
    let resolution = resolution(200.0, 64);
    let omega = resolution.band_rad_s[0] / 0.99;
    let start = std::time::Instant::now();
    let model = source(8, 3, 0.002, ModalEnd::Ports)
        .into_viscothermal_harmonic(resolution, budget(), || false)
        .unwrap();
    assert_eq!(model.solve_dimension(), 54);
    assert!(model.source().state_dimension() > MAX_HARMONIC_DIMENSION);
    assert_eq!(
        model
            .source()
            .viscothermal_model()
            .unwrap()
            .series_terms()
            .len(),
        64
    );
    assert_eq!(
        model
            .source()
            .viscothermal_model()
            .unwrap()
            .shunt_terms()
            .len(),
        64
    );
    let response = model
        .factor(omega, || false)
        .unwrap()
        .response(
            &[
                C64::from_re(1e-7),
                C64::from_re(-2e-7),
                C64::ZERO,
                C64::ZERO,
                C64::new(0.0, 1e-7),
                C64::ZERO,
            ],
            || false,
        )
        .unwrap();
    assert_eq!(response.state.len(), 1142);
    assert!(response.state[54..].iter().any(|v| v.abs() > 0.0));
    close(
        response.supplied_power_w,
        response.heat.viscous_w + response.heat.thermal_w,
        2e-8,
    );
    println!(
        "modal wall: original_states={}, factor_dimension={}, residual={:e}, elapsed_ms={}",
        response.state.len(),
        model.solve_dimension(),
        response.relative_residual,
        start.elapsed().as_millis()
    );
}

#[test]
fn condensed_plane_response_spatially_converges_to_independent_bessel_duct() {
    let gas = gas();
    let resolution = resolution(100.0, 32);
    let omega = resolution.band_rad_s[0] / 0.99;
    let length = 0.7 * gas.sound_speed / omega;
    let duct = crate::Duct {
        segments: vec![crate::Segment::Cylinder {
            radius: 0.002,
            length,
        }],
    };
    let reference =
        C64::from_re(gas.characteristic_impedance / (core::f64::consts::PI * 0.002_f64.powi(2)));
    for (end, termination) in [
        (ModalEnd::Rigid, crate::Termination::Closed),
        (ModalEnd::PressureRelease, crate::Termination::IdealOpen),
    ] {
        let exact = crate::input_impedance(&duct, &gas, omega, LossModel::Bessel, termination)
            .unwrap()
            .impedance;
        let exact_reflection = (exact - reference) / (exact + reference);
        let mut errors = Vec::new();
        for elements in [2, 4, 8] {
            let model = source(elements, 1, length, end)
                .into_viscothermal_harmonic(resolution, budget(), || false)
                .unwrap();
            let response = model
                .factor(omega, || false)
                .unwrap()
                .response(&[C64::ONE], || false)
                .unwrap();
            let z = response.port_pressure_pa[0];
            errors.push(((z - reference) / (z + reference) - exact_reflection).abs());
        }
        assert!(errors[0] > 1e-4, "{errors:?}");
        assert!(
            errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 0.001, "{errors:?}");
    }
}

#[test]
fn pole_ladder_keeps_the_same_acoustic_modes_and_signed_wall_channels() {
    let solve = |poles| {
        let resolution = resolution(100.0, poles);
        let omega = resolution.band_rad_s[0] / 0.99;
        let model = source(6, 2, 0.006, ModalEnd::Ports)
            .into_viscothermal_harmonic(resolution, budget(), || false)
            .unwrap();
        assert_eq!(model.source().modes(), ModeCount(2));
        model
            .factor(omega, || false)
            .unwrap()
            .response(
                &[
                    C64::from_re(1e-7),
                    C64::from_re(-1e-7),
                    C64::ZERO,
                    C64::ZERO,
                ],
                || false,
            )
            .unwrap()
            .port_pressure_pa
    };
    let exact = solve(32);
    let scale = exact.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let errors: Vec<_> = [8, 12, 16]
        .into_iter()
        .map(|poles| {
            solve(poles)
                .iter()
                .zip(&exact)
                .map(|(a, b)| (*a - *b).abs() / scale)
                .fold(0.0_f64, f64::max)
        })
        .collect();
    assert!(errors[0] > 1e-5, "{errors:?}");
    assert!(
        errors[1] < 0.25 * errors[0] && errors[2] < 0.1 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 1e-4, "{errors:?}");
}

#[test]
fn condensed_solver_preserves_band_capacity_and_cancellation_refusals() {
    let resolution = resolution(100.0, 16);
    let make = || source(2, 2, 0.006, ModalEnd::Ports);
    let mut small = budget();
    small.max_dimension = 11;
    assert!(
        make()
            .into_viscothermal_harmonic(resolution, small, || false)
            .is_err()
    );
    small = budget();
    small.max_dense_bytes = 1;
    assert!(
        make()
            .into_viscothermal_harmonic(resolution, small, || false)
            .is_err()
    );
    small = budget();
    small.max_dense_work = 1;
    assert!(
        make()
            .into_viscothermal_harmonic(resolution, small, || false)
            .is_err()
    );
    assert!(
        make()
            .into_viscothermal_harmonic(resolution, budget(), || true)
            .is_err()
    );
    let model = make()
        .into_viscothermal_harmonic(resolution, budget(), || false)
        .unwrap();
    for omega in [
        0.0,
        f64::NAN,
        resolution.band_rad_s[0] * 0.999,
        resolution.band_rad_s[1] * 1.001,
    ] {
        assert!(model.factor(omega, || false).is_err());
    }
    let omega = resolution.band_rad_s[0];
    assert!(matches!(
        model.factor(omega, || true),
        Err(HarmonicError::Cancelled)
    ));
    let mut calls = 0;
    model
        .factor(omega, || {
            calls += 1;
            false
        })
        .unwrap();
    let mut current = 0;
    assert!(matches!(
        model.factor(omega, || {
            current += 1;
            current == calls
        }),
        Err(HarmonicError::Cancelled)
    ));
    let factor = model.factor(omega, || false).unwrap();
    assert!(factor.response(&[], || false).is_err());
    assert!(
        factor
            .response(&[C64::new(f64::NAN, 0.0); 4], || false)
            .is_err()
    );
    assert!(matches!(
        factor.response(&[C64::ONE; 4], || true),
        Err(HarmonicError::Cancelled)
    ));
    let zero = factor.response(&[C64::ZERO; 4], || false).unwrap();
    assert!(zero.state.iter().all(|v| *v == C64::ZERO));
    assert_eq!(zero.heat.viscous_w, 0.0);
    assert_eq!(zero.heat.thermal_w, 0.0);
    assert_eq!(zero.supplied_power_w, 0.0);
}
