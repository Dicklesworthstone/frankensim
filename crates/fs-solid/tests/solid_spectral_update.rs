//! Parameter changes retain physical vibration, clock and independent accounts.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::spectral::{
    CompleteSpectralSolidControls, CompleteSpectralSolidError, DampedSpectralSolid,
};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::{Coo, Csr};

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 64,
                kernel_id: 35,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn controls() -> CompleteSpectralSolidControls {
    CompleteSpectralSolidControls {
        maximum_free_dofs: 128,
        maximum_preparation_bytes: 8 * 1024 * 1024,
        maximum_preparation_multiply_terms: 100_000_000,
        maximum_step_multiply_terms: 1_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_absolute_eigen_residual: 1e-13,
        maximum_relative_state_error: 1e-12,
        maximum_work_error_j: 1e-10,
        maximum_steps: 1000,
    }
}
fn matrix(values: &[f64], n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            coo.push(i, j, values[i * n + j]);
        }
    }
    coo.assemble()
}
fn pencil(k: &[f64], m: &[f64], n: usize) -> TetElasticAssembly {
    TetElasticAssembly {
        stiffness: matrix(k, n),
        mass: matrix(m, n),
        free_dofs: (0..n).collect(),
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0 / 6.0],
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} +/- {tolerance:e}"
    );
}
fn state_bits(s: &DampedSpectralSolid) -> Vec<u64> {
    s.displacement_m()
        .iter()
        .chain(s.velocity_mps())
        .copied()
        .chain([
            s.time_s(),
            s.kinetic_energy_j(),
            s.cumulative_work_j(),
            s.cumulative_dissipated_j(),
        ])
        .map(f64::to_bits)
        .chain([s.steps()])
        .collect()
}
fn all_bits(s: &DampedSpectralSolid) -> Vec<u64> {
    let mut result = state_bits(s);
    result.extend(
        s.damping_rates_per_s()
            .iter()
            .copied()
            .chain([
                s.stored_energy_j(),
                s.cumulative_parameter_work_j(),
                s.relative_damping_projection_error(),
            ])
            .map(f64::to_bits),
    );
    result
}
// Phi=[[1,1/2],[0,1]]. Original matrices are Phi^-T diag(d) Phi^-1;
// its off-diagonal consistent mass is not interchangeable with a lumped mass.
fn operator(a: f64, b: f64) -> [f64; 4] {
    [a, -0.5 * a, -0.5 * a, 0.25 * a + b]
}
fn ring_down(lambda: f64, rate: f64, q: f64, v: f64, t: f64) -> (f64, f64) {
    if lambda == 0.0 {
        let decay = (-rate * t).exp();
        return (q - v * (-rate * t).exp_m1() / rate, v * decay);
    }
    let half = 0.5 * rate;
    let frequency = (lambda - half * half).sqrt();
    let (sine, cosine) = (frequency * t).sin_cos();
    let decay = (-half * t).exp();
    (
        decay * (q * cosine + (v + half * q) * sine / frequency),
        decay * (v * cosine - (half * v + lambda * q) * sine / frequency),
    )
}

#[test]
fn degenerate_rotations_crossings_and_signed_modes_keep_the_physical_state() {
    with_cx(false, |cx| {
        let (mass, stiffness, viscosity) =
            (operator(1.0, 1.0), operator(4.0, 4.0), operator(0.5, 0.5));
        let mut solid = DampedSpectralSolid::prepare_with_basis(
            &pencil(&stiffness, &mass, 2),
            &matrix(&viscosity, 2),
            &[4.0, 4.0],
            &[1.0, 0.5, 0.0, 1.0],
            0.07,
            &[0.3, -0.2],
            &[0.1, 0.25],
            controls(),
            cx,
        )
        .unwrap();
        solid
            .advance_affine(&[0.7, -0.2], &[-0.1, 0.6], cx)
            .unwrap();
        let before = state_bits(&solid);
        let stored = solid.stored_energy_j();
        let previous_heat = solid.cumulative_dissipated_j();
        // A genuine rotation within the repeated eigenspace, not sign matching.
        let rotated = [1.0, -0.5, 0.8, 0.6];
        let report = solid
            .replace_parameters(
                &matrix(&stiffness, 2),
                &matrix(&viscosity, 2),
                &[4.0, 4.0],
                &rotated,
                0.0,
                cx,
            )
            .unwrap();
        assert_eq!(state_bits(&solid), before);
        assert_eq!(solid.stored_energy_j().to_bits(), stored.to_bits());
        assert_eq!(report.parameter_work_j, 0.0);
        assert_eq!(solid.cumulative_parameter_work_j(), 0.0);
        let q = [
            solid.displacement_m()[0] - 0.5 * solid.displacement_m()[1],
            solid.displacement_m()[1],
        ];
        let v = [
            solid.velocity_mps()[0] - 0.5 * solid.velocity_mps()[1],
            solid.velocity_mps()[1],
        ];
        // The first physical mode stiffens, the second becomes truly rigid.
        // New runtime columns are [-phi_1,phi_0], so sorted index copying fails.
        let changed = operator(9.0, 0.0);
        let report = solid
            .replace_parameters(
                &matrix(&changed, 2),
                &matrix(&viscosity, 2),
                &[0.0, 9.0],
                &[-0.5, 1.0, -1.0, 0.0],
                10.0,
                cx,
            )
            .unwrap();
        assert_eq!(state_bits(&solid), before);
        let expected_work = 2.5 * q[0] * q[0] - 2.0 * q[1] * q[1];
        near(report.parameter_work_j, expected_work, 3e-16);
        near(solid.cumulative_parameter_work_j(), expected_work, 3e-16);
        near(report.balance_error_j, 0.0, 3e-16);
        let frame = solid.advance(&[0.0, 0.0], cx).unwrap();
        let first = ring_down(9.0, 0.5, q[0], v[0], 0.07);
        let second = ring_down(0.0, 0.5, q[1], v[1], 0.07);
        near(solid.displacement_m()[0], first.0 + 0.5 * second.0, 3e-15);
        near(solid.displacement_m()[1], second.0, 3e-15);
        near(solid.velocity_mps()[0], first.1 + 0.5 * second.1, 4e-15);
        near(solid.velocity_mps()[1], second.1, 4e-15);
        near(frame.energy.cumulative_balance_error_j, 0.0, 2e-14);
        assert!(solid.cumulative_dissipated_j() > previous_heat);
    });
}

#[test]
fn zero_amplitude_and_damping_events_preserve_old_heat_through_free_motion_and_impulse() {
    with_cx(false, |cx| {
        let mut solid = DampedSpectralSolid::prepare_with_basis(
            &pencil(&[16.0], &[4.0], 1),
            &matrix(&[2.0], 1),
            &[4.0],
            &[0.5],
            0.1,
            &[0.0],
            &[0.3],
            controls(),
            cx,
        )
        .unwrap();
        let before = state_bits(&solid);
        let zero = solid
            .replace_parameters(
                &matrix(&[36.0], 1),
                &matrix(&[4.0], 1),
                &[9.0],
                &[-0.5],
                0.0,
                cx,
            )
            .unwrap();
        assert_eq!(zero.parameter_work_j, 0.0);
        assert_eq!(state_bits(&solid), before);
        solid.advance(&[2.0], cx).unwrap();
        let before = state_bits(&solid);
        let (heat, work) = (solid.cumulative_dissipated_j(), solid.cumulative_work_j());
        assert!(heat > 0.0 && work > 0.0);
        let no_drag = solid
            .replace_parameters(
                &matrix(&[36.0], 1),
                &matrix(&[0.0], 1),
                &[9.0],
                &[0.5],
                0.0,
                cx,
            )
            .unwrap();
        assert_eq!(no_drag.parameter_work_j, 0.0);
        assert_eq!(state_bits(&solid), before);
        let (u, v) = (solid.displacement_m()[0], solid.velocity_mps()[0]);
        let release = solid
            .replace_parameters(
                &matrix(&[0.0], 1),
                &matrix(&[0.0], 1),
                &[0.0],
                &[0.5],
                10.0,
                cx,
            )
            .unwrap();
        near(release.parameter_work_j, -18.0 * u * u, 2e-17);
        assert_eq!(state_bits(&solid), before);
        let free = solid.advance_affine(&[0.8], &[-0.2], cx).unwrap();
        let expected_u = u + 0.1 * v + 0.8 * 0.01 / 8.0 - 0.01 / 24.0;
        let expected_v = v + 0.6 * 0.1 / 8.0;
        near(solid.displacement_m()[0], expected_u, 2e-16);
        near(solid.velocity_mps()[0], expected_v, 2e-16);
        assert_eq!(free.energy.dissipated_j, 0.0);
        assert_eq!(
            free.energy.cumulative_dissipated_j.to_bits(),
            heat.to_bits()
        );
        near(free.energy.cumulative_balance_error_j, 0.0, 2e-15);
        let parameter_total = solid.cumulative_parameter_work_j();
        let impulse = solid.apply_momentum_impulse(&[0.8], 1e-13, cx).unwrap();
        near(solid.velocity_mps()[0], expected_v + 0.2, 2e-16);
        near(
            impulse.energy.supplied_work_j,
            0.8 * (expected_v + 0.1),
            2e-16,
        );
        assert_eq!(
            impulse.energy.cumulative_dissipated_j.to_bits(),
            heat.to_bits()
        );
        assert_eq!(
            solid.cumulative_parameter_work_j().to_bits(),
            parameter_total.to_bits()
        );
        near(impulse.energy.cumulative_balance_error_j, 0.0, 2e-15);
        let (u, v) = (solid.displacement_m()[0], solid.velocity_mps()[0]);
        let reclamp = solid
            .replace_parameters(
                &matrix(&[36.0], 1),
                &matrix(&[0.0], 1),
                &[9.0],
                &[-0.5],
                10.0,
                cx,
            )
            .unwrap();
        near(reclamp.parameter_work_j, 18.0 * u * u, 2e-17);
        let frame = solid.advance(&[0.0], cx).unwrap();
        near(
            solid.displacement_m()[0],
            u * 0.3_f64.cos() + v * 0.3_f64.sin() / 3.0,
            2e-16,
        );
        near(
            solid.velocity_mps()[0],
            v * 0.3_f64.cos() - 3.0 * u * 0.3_f64.sin(),
            3e-16,
        );
        assert_eq!(
            frame.energy.cumulative_dissipated_j.to_bits(),
            heat.to_bits()
        );
        near(frame.energy.cumulative_balance_error_j, 0.0, 2e-15);
        assert_eq!(solid.steps(), 3);
        near(solid.time_s(), 0.3, 6e-17);
    });
}

#[test]
fn rejected_parameter_candidates_never_replace_state_operators_or_history() {
    with_cx(false, |cx| {
        let (mass, stiffness, viscosity) =
            (operator(1.0, 1.0), operator(4.0, 9.0), operator(0.5, 0.5));
        let mut solid = DampedSpectralSolid::prepare_with_basis(
            &pencil(&stiffness, &mass, 2),
            &matrix(&viscosity, 2),
            &[4.0, 9.0],
            &[1.0, 0.5, 0.0, 1.0],
            0.1,
            &[0.4, 0.1],
            &[-0.2, 0.3],
            controls(),
            cx,
        )
        .unwrap();
        solid.advance(&[0.0, 0.0], cx).unwrap();
        let before = all_bits(&solid);
        let changed = matrix(&operator(8.0, 18.0), 2);
        // Late physical-work refusal after complete candidate preparation.
        assert!(matches!(
            solid.replace_parameters(
                &changed,
                &matrix(&viscosity, 2),
                &[8.0, 18.0],
                &[1.0, 0.5, 0.0, 1.0],
                0.0,
                cx
            ),
            Err(CompleteSpectralSolidError::Energy { .. })
        ));
        assert_eq!(all_bits(&solid), before);
        // Euclidean normalization cannot replace the original mass norm.
        assert!(
            solid
                .replace_parameters(
                    &matrix(&stiffness, 2),
                    &matrix(&viscosity, 2),
                    &[4.0, 9.0],
                    &[1.0, 0.0, 0.0, 1.0],
                    10.0,
                    cx
                )
                .is_err()
        );
        assert!(
            solid
                .replace_parameters(
                    &changed,
                    &matrix(&[0.5, 0.0, 0.0, 0.5], 2),
                    &[8.0, 18.0],
                    &[1.0, 0.5, 0.0, 1.0],
                    10.0,
                    cx
                )
                .is_err()
        );
        assert!(
            solid
                .replace_parameters(
                    &changed,
                    &matrix(&viscosity, 2),
                    &[8.0],
                    &[1.0, 0.5, 0.0, 1.0],
                    10.0,
                    cx
                )
                .is_err()
        );
        assert!(
            solid
                .replace_parameters(
                    &changed,
                    &matrix(&viscosity, 2),
                    &[8.0, 18.0],
                    &[1.0, 0.5, 0.0, 1.0],
                    f64::NAN,
                    cx
                )
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert!(matches!(
                solid.replace_parameters(
                    &changed,
                    &matrix(&viscosity, 2),
                    &[8.0, 18.0],
                    &[1.0, 0.5, 0.0, 1.0],
                    10.0,
                    cancelled
                ),
                Err(CompleteSpectralSolidError::Cancelled)
            ))
        });
        assert_eq!(all_bits(&solid), before);
        let q = [
            solid.displacement_m()[0] - 0.5 * solid.displacement_m()[1],
            solid.displacement_m()[1],
        ];
        let v = [
            solid.velocity_mps()[0] - 0.5 * solid.velocity_mps()[1],
            solid.velocity_mps()[1],
        ];
        solid.advance(&[0.0, 0.0], cx).unwrap();
        let first = ring_down(4.0, 0.5, q[0], v[0], 0.1);
        let second = ring_down(9.0, 0.5, q[1], v[1], 0.1);
        near(solid.displacement_m()[0], first.0 + 0.5 * second.0, 3e-15);
        near(solid.velocity_mps()[1], second.1, 3e-15);
        let mut limit = controls();
        limit.maximum_preparation_bytes = 6000;
        let mut small = DampedSpectralSolid::prepare_with_basis(
            &pencil(&stiffness, &mass, 2),
            &matrix(&viscosity, 2),
            &[4.0, 9.0],
            &[1.0, 0.5, 0.0, 1.0],
            0.1,
            &[0.4, 0.1],
            &[-0.2, 0.3],
            limit,
            cx,
        )
        .unwrap();
        let before = all_bits(&small);
        assert!(matches!(
            small.replace_parameters(
                &changed,
                &matrix(&viscosity, 2),
                &[8.0, 18.0],
                &[1.0, 0.5, 0.0, 1.0],
                10.0,
                cx
            ),
            Err(CompleteSpectralSolidError::Budget { .. })
        ));
        assert_eq!(all_bits(&small), before);
    });
}

#[test]
fn parameter_work_does_not_subtract_large_unchanged_stored_energy() {
    with_cx(false, |cx| {
        let mut limit = controls();
        limit.maximum_work_error_j = 2.0; // Explicit endpoint storage-rounding allowance.
        let mut solid = DampedSpectralSolid::prepare_with_basis(
            &pencil(&[1e16, 0.0, 0.0, 1.0], &[1.0, 0.0, 0.0, 1.0], 2),
            &matrix(&[0.0; 4], 2),
            &[1e16, 1.0],
            &[1.0, 0.0, 0.0, 1.0],
            1e-9,
            &[1.0, 1.0],
            &[0.0, 0.0],
            limit,
            cx,
        )
        .unwrap();
        let before = all_bits(&solid);
        let stiffness = matrix(&[1e16, 0.0, 0.0, 3.0], 2);
        assert!(
            solid
                .replace_parameters(
                    &stiffness,
                    &matrix(&[0.0; 4], 2),
                    &[1e16, 3.0],
                    &[1.0, 0.0, 0.0, 1.0],
                    0.75,
                    cx
                )
                .is_err()
        );
        assert_eq!(all_bits(&solid), before);
        let report = solid
            .replace_parameters(
                &stiffness,
                &matrix(&[0.0; 4], 2),
                &[1e16, 3.0],
                &[1.0, 0.0, 0.0, 1.0],
                1.0,
                cx,
            )
            .unwrap();
        assert_eq!(report.parameter_work_j, 1.0);
        assert_eq!(solid.cumulative_parameter_work_j(), 1.0);
        assert_ne!(
            report.stored_after_j - report.stored_before_j,
            report.parameter_work_j
        );
        assert_eq!(solid.cumulative_work_j(), 0.0);
        assert_eq!(solid.cumulative_dissipated_j(), 0.0);
    });
}
