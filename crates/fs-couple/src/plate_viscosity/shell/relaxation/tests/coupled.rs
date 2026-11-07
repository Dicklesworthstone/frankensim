//! Complete public attachment coverage; uses the original coupled runtime.
use super::*;
use crate::render::schedule::force::coupled::ModalCouplingConfig;
use fs_math::c64::C64;
use fs_modal::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};

fn coupled(omega: &[f64], intrinsic_damping: f64, gate: &CancelGate) -> CoupledModalSystem {
    let mut shell = ModalAcousticTimeModel::try_new(
        192_000,
        omega
            .iter()
            .map(|&frequency| ModalAcousticMode {
                angular_frequency_rad_s: frequency,
                damping_ratio: intrinsic_damping,
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect(),
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    shell
        .restore_states(
            &(0..omega.len())
                .map(|i| ModalAcousticState {
                    displacement_m_sqrt_kg: (i as f64 + 1.0) * 1e-6,
                    velocity_m_sqrt_kg_per_s: 0.0,
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
    let mut preceding = ModalAcousticTimeModel::try_new_with_free_coordinates(
        192_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 0.0,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    preceding
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.002,
            velocity_m_sqrt_kg_per_s: 0.0,
        }])
        .unwrap();
    let mut projection = vec![0.0; 2 * (omega.len() + 1)];
    projection[0] = 2.0;
    CoupledModalSystem::new(
        vec![preceding, shell],
        vec![],
        ModalCouplingConfig {
            max_modes: 12,
            max_connections: 8,
            max_setup_terms: 20_000_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e6,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-11,
            energy_relative_tolerance: 1e-8,
        },
        gate,
    )
    .unwrap()
    .with_material_relaxation(
        vec![RelaxationBranch {
            projection,
            stiffness: 1.0,
            relaxation_time_s: 0.1,
        }],
        Some(vec![0.001]),
        gate,
    )
    .unwrap()
}

#[test]
fn public_shell_binding_preserves_prior_histories_and_both_initial_material_states() {
    let (mesh, model, pairs, viscosity) = fixture(CstDktMass::Consistent);
    let phi: Vec<_> = pairs.iter().map(|mode| mode.phi.clone()).collect();
    let omega: Vec<_> = pairs.iter().map(|mode| mode.lambda.sqrt()).collect();
    let gate = CancelGate::new_clock_free();
    let image = ModalShellRelaxation::from_cst_dkt_shell(
        &mesh,
        &spectra(&viscosity),
        &model,
        &phi,
        limits(),
        &gate,
    )
    .unwrap();
    for initial in [InitialPlateMemory::Relaxed, InitialPlateMemory::Unrelaxed] {
        let base = coupled(&omega, 0.0, &gate);
        let earlier_energy = base.material_energy_j().unwrap();
        let q: Vec<_> = base.components()[1]
            .states()
            .iter()
            .map(|s| s.displacement_m_sqrt_kg)
            .collect();
        let mut system = image.attach_to(base, 1, initial, &gate).unwrap();
        assert_eq!(system.material_memory()[0], 0.001);
        assert_eq!(
            system.material_memory().len(),
            1 + image.coupling_rows_per_s().len()
        );
        let mut added_energy = 0.0;
        for ((row, &actual), arm) in image
            .coupling_rows_per_s()
            .iter()
            .zip(&system.material_memory()[1..])
            .zip(&system.material_relaxation_branches()[1..])
        {
            let extension: f64 = row.iter().zip(&q).map(|(c, q)| c * q).sum();
            let expected = if initial == InitialPlateMemory::Relaxed {
                extension
            } else {
                0.0
            };
            assert_eq!(actual, expected);
            added_energy += 0.5 * (extension - expected).powi(2);
            assert_eq!(arm.projection[0], 0.0);
            assert_eq!(arm.projection[1], 0.0);
        }
        near(
            system.material_energy_j().unwrap(),
            earlier_energy + added_energy,
            1e-12,
            earlier_energy + added_energy,
        );
        let initial_energy = system.total_energy_j().unwrap();
        let mut heat = 0.0;
        for _ in 0..16 {
            let frame = system.step(&vec![0.0; image.mode_count() + 1]).unwrap();
            assert!(frame.material_dissipation_j >= 0.0);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert_eq!(frame.component_dissipation_j, 0.0);
            assert_eq!(frame.connection_dissipation_j, 0.0);
            heat += frame.material_dissipation_j;
        }
        near(
            system.total_energy_j().unwrap() + heat,
            initial_energy,
            1e-8,
            initial_energy,
        );
        assert!(image.attach_to(system, 1, initial, &gate).is_err());
    }
    assert!(
        image
            .attach_to(
                coupled(&omega, 0.01, &gate),
                1,
                InitialPlateMemory::Relaxed,
                &gate
            )
            .is_err()
    );
}
