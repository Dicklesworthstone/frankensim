use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::plate::impact::radiation::{Model, Pole};

fn mode(omega: f64, zeta: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: zeta,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}

fn network(air_loss: f64) -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let config = ModalCouplingConfig {
        max_modes: 16,
        max_connections: 8,
        max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-9,
    };
    CoupledModalSystem::new_with_columns_preparation(
        vec![mode(500.0, 0.02), mode(900.0, 0.03)],
        vec![ModalColumnConnection {
            column: vec![1.0, -0.4],
            stiffness_n_m: 50_000.0,
            damping_n_s_m: 12.0,
            rest_extension_m: 1e-5,
        }],
        &[ConstantLossConnection {
            connection: 0,
            reference_omega_rad_s: 700.0,
        }],
        config,
        ModalPreparation::HarmonicOnly(PreparedModalBudget {
            max_modes: 2,
            max_dense_bytes: 2_000_000,
            max_multiply_terms: 100_000_000,
        }),
        &gate,
    )
    .unwrap()
    .with_material_relaxation(
        vec![
            fs_phs::RelaxationBranch {
                projection: vec![1.0, 0.0, -0.3, 0.0],
                stiffness: 30_000.0,
                relaxation_time_s: 0.002,
            },
            fs_phs::RelaxationBranch {
                projection: vec![0.4, 0.0, 0.8, 0.0],
                stiffness: 8_000.0,
                relaxation_time_s: 0.0007,
            },
        ],
        Some(vec![0.004, -0.007]),
        &gate,
    )
    .unwrap()
    .with_radiation_load(
        &Model {
            ports: 2,
            poles: vec![Pole {
                omega: 700.0,
                zeta: air_loss,
                coupling: vec![100.0, -30.0],
            }],
        },
        &[0, 1],
        1,
        &gate,
    )
    .unwrap()
}

fn close_scalar(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-7 * actual.abs().max(expected.abs()) + 2e-15,
        "{actual} != {expected}"
    );
}

fn compare_vectors(derivative: &[C64], minus: &[C64], plus: &[C64], step: f64) {
    assert_eq!(derivative.len(), minus.len());
    assert_eq!(derivative.len(), plus.len());
    for ((&actual, &lo), &hi) in derivative.iter().zip(minus).zip(plus) {
        let expected = (hi - lo).scale(0.5 / step);
        close_scalar(actual.re, expected.re);
        close_scalar(actual.im, expected.im);
    }
}

#[test]
fn coupled_frequency_derivative_matches_the_complete_frequency_dependent_physics() {
    // Independent central differences check both the implicit solve and every
    // explicit-frequency observation. Anchored damping has a different direction
    // from its spring; two signed hereditary arms and reciprocal air all remain.
    let gate = CancelGate::new_clock_free();
    let system = network(0.03);
    let memories = system.material_memory().to_vec();
    let air = system.radiation_state().to_vec();
    let force = [C64::new(0.3, 0.1), C64::new(-0.2, 0.05)];
    let slope = [C64::new(0.0003, -0.0001), C64::new(-0.0002, 0.0004)];
    let step = 0.0005;
    for omega in [550.0, 650.0, 700.0, 980.0] {
        let result = system
            .harmonic_response_with_frequency_derivative(omega, &force, Some(&slope), 10_000, &gate)
            .unwrap();
        assert_eq!(
            result.response,
            system
                .harmonic_response(omega, &force, 10_000, &gate)
                .unwrap()
        );
        let minus_force: Vec<_> = force
            .iter()
            .zip(&slope)
            .map(|(&f, &df)| f - df.scale(step))
            .collect();
        let plus_force: Vec<_> = force
            .iter()
            .zip(&slope)
            .map(|(&f, &df)| f + df.scale(step))
            .collect();
        let lo = system
            .harmonic_response(omega - step, &minus_force, 10_000, &gate)
            .unwrap();
        let hi = system
            .harmonic_response(omega + step, &plus_force, 10_000, &gate)
            .unwrap();
        for (actual, low, high) in [
            (
                &result.displacement_derivative,
                &lo.displacement_m_sqrt_kg,
                &hi.displacement_m_sqrt_kg,
            ),
            (
                &result.velocity_derivative,
                &lo.velocity_m_sqrt_kg_per_s,
                &hi.velocity_m_sqrt_kg_per_s,
            ),
            (
                &result.material_memory_derivative,
                &lo.material_memory,
                &hi.material_memory,
            ),
            (
                &result.radiation_state_derivative,
                &lo.radiation_state,
                &hi.radiation_state,
            ),
            (
                &result.connection_forces_derivative,
                &lo.connection_forces_n,
                &hi.connection_forces_n,
            ),
            (
                &result.connection_spring_forces_derivative,
                &lo.connection_spring_forces_n,
                &hi.connection_spring_forces_n,
            ),
            (
                &result.connection_damping_forces_derivative,
                &lo.connection_damping_forces_n,
                &hi.connection_damping_forces_n,
            ),
            (
                &result.generalized_connection_forces_derivative,
                &lo.generalized_connection_forces_n_per_sqrt_kg,
                &hi.generalized_connection_forces_n_per_sqrt_kg,
            ),
        ] {
            compare_vectors(actual, low, high, step);
        }
        for (actual, low, high) in [
            (
                result.mean_stored_energy_derivative,
                lo.mean_stored_energy_j,
                hi.mean_stored_energy_j,
            ),
            (
                result.external_power_derivative,
                lo.external_power_w,
                hi.external_power_w,
            ),
            (
                result.component_dissipation_derivative,
                lo.component_dissipation_w,
                hi.component_dissipation_w,
            ),
            (
                result.connection_dissipation_derivative,
                lo.connection_dissipation_w,
                hi.connection_dissipation_w,
            ),
            (
                result.material_dissipation_derivative,
                lo.material_dissipation_w,
                hi.material_dissipation_w,
            ),
            (
                result.radiation_dissipation_derivative,
                lo.radiation_dissipation_w,
                hi.radiation_dissipation_w,
            ),
        ] {
            close_scalar(actual, (high - low) * 0.5 / step);
        }
        for (actual, low, high) in [
            (
                &result.material_branch_dissipation_derivative,
                &lo.material_branch_dissipation_w,
                &hi.material_branch_dissipation_w,
            ),
            (
                &result.material_branch_mean_storage_derivative,
                &lo.material_branch_mean_storage_j,
                &hi.material_branch_mean_storage_j,
            ),
        ] {
            for ((&value, &low), &high) in actual.iter().zip(low).zip(high) {
                close_scalar(value, (high - low) * 0.5 / step);
            }
        }
        assert!(result.relative_derivative_residual <= 1e-11);
        assert!(
            result.power_residual_derivative.abs()
                <= 5e-11 * result.external_power_derivative.abs() + 1e-18
        );
        let observation = [0.4, -0.7];
        let expected = projection(&observation, &result.displacement_derivative)
            / projection(&observation, &result.response.displacement_m_sqrt_kg);
        assert_eq!(
            result
                .displacement_logarithmic_derivative(&observation)
                .unwrap(),
            Some(expected)
        );
    }
    assert_eq!(system.material_memory(), memories);
    assert_eq!(system.radiation_state(), air);
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn zero_force_with_nonzero_drive_slope_has_a_linear_tangent_and_zero_power_slope() {
    let gate = CancelGate::new_clock_free();
    let system = network(0.0);
    let slope = [C64::new(0.3, 0.1), C64::new(-0.2, 0.05)];
    let result = system
        .harmonic_response_with_frequency_derivative(
            700.0,
            &[C64::ZERO; 2],
            Some(&slope),
            10_000,
            &gate,
        )
        .unwrap();
    let expected = system
        .harmonic_response(700.0, &slope, 10_000, &gate)
        .unwrap();
    assert_eq!(
        result.displacement_derivative,
        expected.displacement_m_sqrt_kg
    );
    assert_eq!(
        result.velocity_derivative,
        expected.velocity_m_sqrt_kg_per_s
    );
    assert_eq!(result.external_power_derivative, 0.0);
    assert_eq!(result.total_dissipation_derivative(), 0.0);
    assert_eq!(result.mean_stored_energy_derivative, 0.0);
    assert_eq!(
        result
            .displacement_logarithmic_derivative(&[1.0, -0.3])
            .unwrap(),
        None
    );
    assert!(result.displacement_logarithmic_derivative(&[1.0]).is_err());
    assert!(
        result
            .displacement_logarithmic_derivative(&[f64::NAN, 0.0])
            .is_err()
    );
}

#[test]
fn tangent_refusals_preserve_all_accepted_material_and_air_histories() {
    let gate = CancelGate::new_clock_free();
    let system = network(0.03);
    let histories = system.material_memory().to_vec();
    let air = system.radiation_state().to_vec();
    let states: Vec<_> = system
        .components()
        .iter()
        .map(|model| model.states().to_vec())
        .collect();
    let force = [C64::new(0.3, 0.1), C64::ZERO];
    // d=3, G=1, C=1, M=2: 3^2*(3+1+2*(1+2+6))=198.
    assert!(
        system
            .harmonic_response_with_frequency_derivative(650.0, &force, None, 197, &gate)
            .is_err()
    );
    system
        .harmonic_response_with_frequency_derivative(650.0, &force, None, 198, &gate)
        .unwrap();
    assert!(
        system
            .harmonic_response_with_frequency_derivative(
                650.0,
                &force,
                Some(&force[..1]),
                10_000,
                &gate
            )
            .is_err()
    );
    assert!(
        system
            .harmonic_response_with_frequency_derivative(
                650.0,
                &force,
                Some(&[C64::ZERO, C64::new(f64::INFINITY, 0.0)]),
                10_000,
                &gate
            )
            .is_err()
    );
    assert!(
        system
            .harmonic_response_with_frequency_derivative(0.0, &force, None, 10_000, &gate)
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.harmonic_response_with_frequency_derivative(650.0, &force, None, 10_000, &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    assert_eq!(system.material_memory(), histories);
    assert_eq!(system.radiation_state(), air);
    assert_eq!(
        states,
        system
            .components()
            .iter()
            .map(|model| model.states().to_vec())
            .collect::<Vec<_>>()
    );
    assert_eq!(system.samples_rendered(), 0);
}
