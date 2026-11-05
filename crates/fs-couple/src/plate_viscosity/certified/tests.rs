use super::*;
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
};
use crate::render::schedule::force::coupled::{
    ModalCouplingError, PreparedModalBudget,
};
use fs_math::c64::C64;

fn pole_budget() -> ModalDampingBudget {
    ModalDampingBudget {
        max_rel_decay_error: 1e-3,
        max_freq_error_cents: 0.01,
        max_cluster_size: 4,
    }
}

fn config(modes: usize, connections: usize) -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: modes,
        max_connections: connections,
        max_setup_terms: 2_000_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 1.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 1e5,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-13,
        energy_relative_tolerance: 1e-9,
    }
}

fn prepared(modes: usize) -> PreparedModalBudget {
    PreparedModalBudget {
        max_modes: modes,
        max_dense_bytes: 64 * 1024 * 1024,
        max_multiply_terms: 1_000_000_000,
    }
}

fn model(omega: &[f64]) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        48_000,
        omega.iter().enumerate().map(|(i, &angular_frequency_rad_s)| ModalAcousticMode {
            angular_frequency_rad_s,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::new(0.25 * (i + 1) as f64, -0.1),
        }).collect(),
        ModalAcousticTimeBudget::audible_reference(),
    ).unwrap()
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "{actual:e} versus {expected:e}",
    );
}

#[test]
fn g1_full_cluster_matches_rotated_scalar_motion_and_independent_heat() {
    let gate = CancelGate::new_clock_free();
    let omega = [100.0, 100.0];
    let damping = ModalPlateViscosity::from_psd_modal_matrix(
        vec![4.0, 3.0, 3.0, 4.0], 2, 1.0, &gate,
    ).unwrap();
    let image = damping.certify_for_time(&omega, &pole_budget(), &gate).unwrap();
    assert!(matches!(
        &image.certificate().blocks[..],
        [ModalDampingBlock::Full { modes, .. }] if modes == &[0, 1]
    ));
    assert_eq!(image.diagonal_rates_per_s(), &[0.0, 0.0]);
    assert_eq!(image.coupled_channel_count(), 2);

    let initial = [
        ModalAcousticState {
            displacement_m_sqrt_kg: 0.001,
            velocity_m_sqrt_kg_per_s: 0.2,
        },
        ModalAcousticState {
            displacement_m_sqrt_kg: -0.0003,
            velocity_m_sqrt_kg_per_s: -0.1,
        },
    ];
    let root2 = core::f64::consts::SQRT_2;
    for preparation in [
        ModalPreparation::Ordinary,
        ModalPreparation::Prepared(prepared(2)),
    ] {
        let mut original = model(&omega);
        original.restore_states(&initial).unwrap();
        let mut actual = image.prepare_coupled(
            original, config(2, 2), preparation, &gate,
        ).unwrap();
        assert_eq!(actual.components()[0].states(), &initial);
        let initial_energy = actual.total_energy_j().unwrap();

        // Independent analytic diagonalization: C has orthogonal [1,1] and
        // [1,-1] directions with rates 7 and 1, while K=100^2 I is unchanged.
        // The reference uses scalar exact formulas, not the matrix exponential.
        let mut reference = model(&omega).with_damping_rates(&[7.0, 1.0]).unwrap();
        reference.restore_states(&[
            ModalAcousticState {
                displacement_m_sqrt_kg: (initial[0].displacement_m_sqrt_kg
                    + initial[1].displacement_m_sqrt_kg) / root2,
                velocity_m_sqrt_kg_per_s: (initial[0].velocity_m_sqrt_kg_per_s
                    + initial[1].velocity_m_sqrt_kg_per_s) / root2,
            },
            ModalAcousticState {
                displacement_m_sqrt_kg: (initial[0].displacement_m_sqrt_kg
                    - initial[1].displacement_m_sqrt_kg) / root2,
                velocity_m_sqrt_kg_per_s: (initial[0].velocity_m_sqrt_kg_per_s
                    - initial[1].velocity_m_sqrt_kg_per_s) / root2,
            },
        ]).unwrap();

        let (mut work, mut heat) = (0.0, 0.0);
        for sample in 0..64 {
            let force = if sample < 32 { [0.4, -0.2] } else { [0.0; 2] };
            let frame = actual.step_under_gate(&force, &gate).unwrap().clone();
            let truth = reference.step(&[
                (force[0] + force[1]) / root2,
                (force[0] - force[1]) / root2,
            ]).unwrap();
            for (i, sign) in [1.0, -1.0].into_iter().enumerate() {
                let a = actual.components()[0].states()[i];
                let p = reference.states()[0];
                let m = reference.states()[1];
                close(a.displacement_m_sqrt_kg,
                    (p.displacement_m_sqrt_kg + sign * m.displacement_m_sqrt_kg) / root2,
                    1e-9, 1e-13);
                close(a.velocity_m_sqrt_kg_per_s,
                    (p.velocity_m_sqrt_kg_per_s + sign * m.velocity_m_sqrt_kg_per_s) / root2,
                    1e-9, 1e-13);
            }
            assert_eq!(frame.component_dissipation_j, 0.0);
            assert_eq!(frame.material_dissipation_j, 0.0);
            assert_eq!(frame.numerical_dissipation_j, 0.0);
            assert!(frame.connection_dissipation_j >= 0.0);
            close(frame.connection_dissipation_j, truth.viscous_dissipation_j, 1e-9, 2e-15);
            close(frame.external_work_j, truth.input_work_j, 1e-9, 2e-15);
            work += frame.external_work_j;
            heat += frame.connection_dissipation_j;
        }
        assert!(heat > 0.0);
        close(actual.total_energy_j().unwrap() + heat - work, initial_energy, 1e-10, 1e-13);
    }
}

#[test]
fn g1_many_mode_dense_and_windowed_images_realize_only_required_clusters() {
    let gate = CancelGate::new_clock_free();
    let n = 96;
    // The first three coordinates are the existing certifier's clustered
    // fixture. Tiny, physically nonzero neighbouring dashpot entries connect
    // the remaining modes, so the original complete matrix has one large island.
    let omega: Vec<_> = (0..n).map(|i| match i {
        0 => 1000.0,
        1 => 1003.0,
        _ => 1700.0 + 100.0 * (i - 2) as f64,
    }).collect();
    let mut matrix = vec![0.0; n * n];
    for i in 0..n { matrix[n * i + i] = 3.0; }
    matrix[0] = 4.0;
    matrix[n + 1] = 0.1;
    for (i, j, value) in [(0, 1, 0.6), (0, 2, 0.05), (1, 2, 0.02)] {
        matrix[n * i + j] = value;
        matrix[n * j + i] = value;
    }
    for i in 2..n - 1 {
        matrix[n * i + i + 1] = 1e-6;
        matrix[n * (i + 1) + i] = 1e-6;
    }
    let damping = ModalPlateViscosity::from_psd_modal_matrix(
        matrix.clone(), n, 2.0, &gate,
    ).unwrap();
    let dense = damping.certify_for_time(&omega, &pole_budget(), &gate).unwrap();

    let window = 1000.0;
    let near: Vec<Vec<(usize, f64)>> = (0..n).map(|i| {
        (0..n).filter(|&j| i == j || (omega[i] - omega[j]).abs() < window)
            .map(|j| (j, matrix[n * i + j])).collect()
    }).collect();
    let far: Vec<_> = (0..n).map(|i| {
        (0..n).filter(|&j| i != j && (omega[i] - omega[j]).abs() >= window)
            .map(|j| matrix[n * i + j].powi(2)).sum::<f64>().sqrt()
    }).collect();
    assert!(far.iter().all(|value| *value == 0.0));
    let source = WindowedPlateViscosity {
        omega: omega.clone(),
        window,
        near,
        far,
        column: (0..n).map(|i| {
            (0..n).map(|j| matrix[n * i + j].powi(2)).sum::<f64>().sqrt()
        }).collect(),
        defect: 0.0,
    };
    let windowed = source.certify_for_time(&pole_budget(), &gate).unwrap();
    for image in [&dense, &windowed] {
        assert_eq!(image.omega_rad_s(), omega);
        assert_eq!(image.clusters.len(), 1);
        assert_eq!(image.clusters[0].modes, [0, 1]);
        assert_eq!(image.clusters[0].damping.matrix_per_s(), &[4.0, 0.6, 0.6, 0.1]);
        assert_eq!(image.coupled_channel_count(), 2);
        assert_eq!(image.diagonal_rates_per_s()[..2], [0.0; 2]);
        assert!(image.diagonal_rates_per_s()[2..].iter().all(|rate| *rate == 3.0));
    }
    assert_eq!(windowed.certificate().far_coupling.window_rad_s, Some(window));
    let original = model(&omega);
    let full = CoupledModalSystem::new(
        vec![original.clone()], damping.connections(0), config(n, n), &gate,
    );
    assert!(matches!(
        full,
        Err(ModalCouplingError::Propagation(fs_phs::PhsError::LinearSplitBudget { .. }))
    ), "the nonzero far links must keep the full source beyond ordinary dense capacity");

    let mut a = dense.prepare_coupled(
        original.clone(), config(n, 2), ModalPreparation::Ordinary, &gate,
    ).unwrap();
    let mut b = windowed.prepare_coupled(
        original.clone(), config(n, 2), ModalPreparation::Ordinary, &gate,
    ).unwrap();
    for system in [&a, &b] {
        assert_eq!(system.mode_count(), n);
        assert_eq!(system.sample_period_s(), original.sample_period_s());
        assert_eq!(system.components()[0].states(), original.states());
        for (left, right) in system.components()[0].modes().iter().zip(original.modes()) {
            assert_eq!(left.angular_frequency_rad_s, right.angular_frequency_rad_s);
            assert_eq!(left.pressure_per_modal_velocity, right.pressure_per_modal_velocity);
        }
        for channel in 0..2 {
            assert!(system.connection_column(channel).unwrap()[2..].iter().all(|v| *v == 0.0));
        }
    }
    let force: Vec<_> = (0..n).map(|i| 0.1 / (i + 1) as f64).collect();
    for _ in 0..4 {
        let fa = a.step(&force).unwrap().clone();
        let fb = b.step(&force).unwrap().clone();
        close(fa.total_dissipation_j(), fb.total_dissipation_j(), 1e-9, 1e-18);
        assert!(fa.component_dissipation_j > 0.0 && fa.connection_dissipation_j > 0.0);
        for (left, right) in a.components()[0].states().iter().zip(b.components()[0].states()) {
            close(left.displacement_m_sqrt_kg, right.displacement_m_sqrt_kg, 1e-9, 1e-17);
            close(left.velocity_m_sqrt_kg_per_s, right.velocity_m_sqrt_kg_per_s, 1e-9, 1e-14);
        }
    }
}

#[test]
fn g0_scalar_owner_and_exact_lossless_modes_keep_native_fast_path() {
    let gate = CancelGate::new_clock_free();
    let omega = [100.0, 250.0, 700.0];
    let zeta = [0.0, 0.001, 0.001 + 1e-17];
    let damping = ModalPlateViscosity::from_psd_modal_matrix(
        vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.4], 3, 1.0, &gate,
    ).unwrap().with_diagonal_zeta(zeta.to_vec());
    let image = damping.certify_for_time(&omega, &pole_budget(), &gate).unwrap();
    assert_eq!(image.coupled_channel_count(), 0);
    for i in 0..3 {
        assert_eq!(image.diagonal_rates_per_s()[i].to_bits(), (2.0 * (zeta[i] * omega[i])).to_bits());
    }
    let original = model(&omega);
    let mut reference = original.with_damping_rates(image.diagonal_rates_per_s()).unwrap();
    let mut system = image.prepare_coupled(
        original, config(3, 0), ModalPreparation::Ordinary, &gate,
    ).unwrap();
    for _ in 0..16 {
        let truth = reference.step(&[0.1, 0.2, -0.3]).unwrap();
        let actual = system.step(&[0.1, 0.2, -0.3]).unwrap().clone();
        assert_eq!(system.components()[0].states(), reference.states());
        assert_eq!(actual.connection_dissipation_j, 0.0);
        assert_eq!(actual.numerical_dissipation_j, 0.0);
        close(actual.component_dissipation_j, truth.viscous_dissipation_j, 1e-12, 1e-20);
        assert_eq!(truth.modal_viscous_losses[0].energy_j, 0.0);
    }
    let repeated = ModalPlateViscosity::from_psd_modal_matrix(vec![0.0; 4], 2, 1.0, &gate)
        .unwrap().certify_for_time(&[100.0; 2], &pole_budget(), &gate).unwrap();
    assert_eq!(repeated.coupled_channel_count(), 0);
    assert_eq!(repeated.diagonal_rates_per_s(), &[0.0; 2]);
}

#[test]
fn g0_binding_refuses_duplicate_loss_wrong_basis_and_original_capacity_limits() {
    let gate = CancelGate::new_clock_free();
    let omega = [1000.0, 1003.0];
    let damping = ModalPlateViscosity::from_psd_modal_matrix(
        vec![4.0, 0.6, 0.6, 0.1], 2, 1.0, &gate,
    ).unwrap();
    assert!(damping.certify_for_time(
        &omega, &ModalDampingBudget { max_cluster_size: 1, ..pole_budget() }, &gate,
    ).is_err());
    let image = damping.certify_for_time(&omega, &pole_budget(), &gate).unwrap();
    let original = model(&omega);
    for invalid_model in [
        original.with_damping_rates(&[0.1, 0.0]).unwrap(),
        model(&[1003.0, 1000.0]),
    ] {
        assert!(image.prepare_coupled(
            invalid_model, config(2, 2), ModalPreparation::Ordinary, &gate,
        ).is_err());
    }
    assert!(image.prepare_coupled(
        original.clone(), config(2, 0), ModalPreparation::Ordinary, &gate,
    ).is_err());
    assert!(image.prepare_coupled(
        original.clone(), config(2, 2), ModalPreparation::Prepared(prepared(1)), &gate,
    ).is_err());
    assert!(image.prepare_coupled(
        original.clone(), config(2, 2), ModalPreparation::HarmonicOnly(prepared(2)), &gate,
    ).is_err());

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        damping.certify_for_time(&omega, &pole_budget(), &cancelled),
        Err(AcousticRealizeError::Cancelled)
    ));
    assert!(matches!(
        image.prepare_coupled(
            original.clone(), config(2, 2), ModalPreparation::Ordinary, &cancelled,
        ),
        Err(AcousticRealizeError::Cancelled)
    ));
    let retry = image.prepare_coupled(
        original.clone(), config(2, 2), ModalPreparation::Ordinary, &gate,
    ).unwrap();
    assert_eq!(retry.samples_rendered(), 0);
    assert_eq!(retry.sample_period_s(), original.sample_period_s());
    assert_eq!(retry.components()[0].states(), original.states());
}
