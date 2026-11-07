use super::*;
use crate::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingConfig};
use fs_math::c64::C64;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
    ModePair,
};
use fs_plate::PlateSection;
use fs_plate::shell::{CstDktMass, ShellSupport, assemble_shell_sections_with_mass};

pub(super) fn fixture(mass: CstDktMass) -> (ShellMesh, ShellModel, Vec<ModePair>, Vec<ShellViscosity>) {
    let mesh = ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [0.2, 0.0, 0.01],
            [0.0, 0.1, 0.02],
            [0.18, 0.09, 0.035],
        ],
        vec![[0, 1, 3], [0, 3, 2]],
    )
    .unwrap();
    let sections = [
        PlateSection::isotropic(1e6, 0.25, 0.005, 1000.0).unwrap(),
        PlateSection::isotropic(2e6, 0.3, 0.003, 1200.0).unwrap(),
    ];
    let model = assemble_shell_sections_with_mass(
        &mesh,
        &sections,
        &[0, 1, 2],
        ShellSupport::Clamped,
        mass,
    )
    .unwrap();
    let n = model.free;
    let mut stiffness = vec![0.0; n * n];
    let mut inertia = stiffness.clone();
    for i in 0..n {
        for j in 0..n {
            stiffness[n * i + j] = model.k.get(i, j);
            inertia[n * i + j] = model.m.get(i, j);
        }
    }
    let modes = fs_modal::eigh_gen_dense(&stiffness, &inertia, n).unwrap();
    let viscosity = vec![
        ShellViscosity::new_coupled(
            [0.8, 0.1, 0.02, 0.1, 1.2, 0.04, 0.02, 0.04, 0.5],
            [0.0001, 0.0, 0.0, 0.0, -0.0002, 0.0, 0.0, 0.0, 0.00005],
            [2e-6, 0.2e-6, 0.0, 0.2e-6, 4e-6, 0.3e-6, 0.0, 0.3e-6, 1e-6],
        )
        .unwrap(),
        ShellViscosity::new(
            [4.0, 0.2, 0.1, 0.2, 0.7, 0.0, 0.1, 0.0, 0.4],
            [8e-6, 0.1e-6, 0.0, 0.1e-6, 1e-6, 0.0, 0.0, 0.0, 3e-6],
        )
        .unwrap(),
    ];
    (mesh, model, modes, viscosity)
}

fn limits() -> PlateViscosityLimits {
    PlateViscosityLimits {
        max_modes: 12,
        max_projection_terms: 2_000_000,
        reference_resistance_n_s_m: 2.5,
    }
}

fn basis(modes: &[ModePair]) -> Vec<Vec<f64>> {
    modes.iter().map(|mode| mode.phi.clone()).collect()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * scale.max(1e-15),
        "{actual:e} != {expected:e}, scale={scale:e}"
    );
}

#[test]
fn full_shell_channels_match_original_sparse_material_power_with_both_mass_models() {
    for mass in [CstDktMass::Lumped, CstDktMass::Consistent] {
        let (mesh, model, modes, viscosity) = fixture(mass);
        let phi = basis(&modes);
        let gate = CancelGate::new_clock_free();
        let law = ModalPlateViscosity::from_cst_dkt_shell(
            &mesh,
            &viscosity,
            &model,
            &phi,
            limits(),
            &gate,
        )
        .unwrap();
        let n = modes.len();
        let matrix = law.matrix_per_s();
        assert!(
            matrix
                .iter()
                .enumerate()
                .any(|(i, value)| i / n != i % n && value.abs() > 1e-5)
        );
        let spatial = fs_plate::shell::viscosity::assemble_constitutive_operator(
            &mesh,
            &viscosity,
            &model,
            ShellViscosityBudget {
                max_elements: 2,
                max_full_dofs: 24,
                max_dofs: 24,
                max_modes: 12,
                max_work: 2_000_000,
            },
            || false,
        )
        .unwrap();
        for velocity in [vec![0.2, -0.1, 0.07, -0.03, 0.04, -0.02], vec![1.0; n]] {
            let mut physical = vec![0.0; model.free];
            for (shape, &rate) in phi.iter().zip(&velocity) {
                for (entry, &value) in physical.iter_mut().zip(shape) {
                    *entry += rate * value;
                }
            }
            let mut force = vec![0.0; model.free];
            spatial.spmv(&physical, &mut force);
            let expected: f64 = physical.iter().zip(force).map(|(v, f)| v * f).sum();
            let modal: f64 = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| velocity[i] * matrix[n * i + j] * velocity[j])
                        .sum::<f64>()
                })
                .sum();
            let channels: f64 = law
                .connections(0)
                .iter()
                .map(|link| {
                    let slip: f64 = link
                        .left
                        .shapes
                        .iter()
                        .zip(&velocity)
                        .map(|(s, v)| s * v)
                        .sum();
                    link.damping_n_s_m * slip * slip
                })
                .sum();
            assert!(expected > 0.0);
            near(modal, expected, expected);
            near(channels, expected, expected);
        }
        let reordered: Vec<_> = phi
            .iter()
            .rev()
            .map(|mode| mode.iter().map(|x| -x).collect())
            .collect();
        let changed = ModalPlateViscosity::from_cst_dkt_shell(
            &mesh,
            &viscosity,
            &model,
            &reordered,
            limits(),
            &gate,
        )
        .unwrap();
        let scale = matrix.iter().map(|v| v.abs()).fold(0.0, f64::max);
        for i in 0..n {
            for j in 0..n {
                near(
                    changed.matrix_per_s()[n * i + j],
                    matrix[n * (n - 1 - i) + n - 1 - j],
                    scale,
                );
            }
        }
    }
}

fn rhs(state: &[f64], omega2: &[f64], damping: &[f64]) -> Vec<f64> {
    let n = omega2.len();
    let mut out = vec![0.0; 2 * n];
    for i in 0..n {
        out[i] = state[n + i];
        out[n + i] = -omega2[i] * state[i];
        for j in 0..n {
            out[n + i] -= damping[n * i + j] * state[n + j];
        }
    }
    out
}

fn rk4(state: &mut [f64], omega2: &[f64], damping: &[f64], dt: f64) {
    let stage = |rate: &[f64], scale| {
        state
            .iter()
            .zip(rate)
            .map(|(x, k)| x + scale * k)
            .collect::<Vec<_>>()
    };
    let k1 = rhs(state, omega2, damping);
    let k2 = rhs(&stage(&k1, 0.5 * dt), omega2, damping);
    let k3 = rhs(&stage(&k2, 0.5 * dt), omega2, damping);
    let k4 = rhs(&stage(&k3, dt), omega2, damping);
    for i in 0..state.len() {
        state[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

#[test]
fn shell_material_loss_drives_existing_coupled_motion_and_closes_physical_heat() {
    let (mesh, model, modes, viscosity) = fixture(CstDktMass::Consistent);
    let gate = CancelGate::new_clock_free();
    let law = ModalPlateViscosity::from_cst_dkt_shell(
        &mesh,
        &viscosity,
        &model,
        &basis(&modes),
        limits(),
        &gate,
    )
    .unwrap();
    let n = modes.len();
    let omega2: Vec<_> = modes.iter().map(|mode| mode.lambda).collect();
    let mut component = ModalAcousticTimeModel::try_new(
        48_000,
        modes
            .iter()
            .map(|mode| ModalAcousticMode {
                angular_frequency_rad_s: mode.lambda.sqrt(),
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect(),
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let initial_velocity = [0.2, -0.1, 0.07, -0.03, 0.04, -0.02];
    component
        .restore_states(
            &initial_velocity
                .iter()
                .map(|&velocity| ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: velocity,
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
    let mut system = CoupledModalSystem::new(
        vec![component],
        law.connections(0),
        ModalCouplingConfig {
            max_modes: 12,
            max_connections: 12,
            max_setup_terms: 2_000_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1.0,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-11,
            energy_relative_tolerance: 1e-8,
        },
        &gate,
    )
    .unwrap();
    let initial = system.total_energy_j().unwrap();
    let mut reference = vec![0.0; 2 * n];
    reference[n..].copy_from_slice(&initial_velocity);
    let mut diagonal_reference = reference.clone();
    let mut diagonal = law.matrix_per_s().to_vec();
    for i in 0..n {
        for j in 0..n {
            if i != j {
                diagonal[n * i + j] = 0.0;
            }
        }
    }
    let mut heat = 0.0;
    for _ in 0..128 {
        let frame = system.step_under_gate(&vec![0.0; n], &gate).unwrap();
        assert!(frame.connection_dissipation_j >= 0.0);
        assert_eq!(frame.component_dissipation_j, 0.0);
        assert_eq!(frame.numerical_dissipation_j, 0.0);
        heat += frame.connection_dissipation_j;
        for _ in 0..16 {
            rk4(
                &mut reference,
                &omega2,
                law.matrix_per_s(),
                1.0 / (48_000.0 * 16.0),
            );
            rk4(
                &mut diagonal_reference,
                &omega2,
                &diagonal,
                1.0 / (48_000.0 * 16.0),
            );
        }
    }
    assert!(heat > 1e-7);
    near(system.total_energy_j().unwrap() + heat, initial, initial);
    for (i, state) in system.components()[0].states().iter().enumerate() {
        near(state.displacement_m_sqrt_kg, reference[i], 0.01);
        near(state.velocity_m_sqrt_kg_per_s, reference[n + i], 1.0);
    }
    assert!(
        reference
            .iter()
            .zip(diagonal_reference)
            .any(|(full, diagonal)| (full - diagonal).abs() > 1e-6)
    );
}

#[test]
fn zero_material_loss_and_invalid_basis_budget_cancellation_are_explicit() {
    let (mesh, model, modes, viscosity) = fixture(CstDktMass::Consistent);
    let gate = CancelGate::new_clock_free();
    let phi = basis(&modes);
    let zeros = vec![ShellViscosity::new([0.0; 9], [0.0; 9]).unwrap(); 2];
    let law = ModalPlateViscosity::from_cst_dkt_shell(&mesh, &zeros, &model, &phi, limits(), &gate)
        .unwrap();
    assert_eq!(law.channel_count(), 0);
    assert!(law.matrix_per_s().iter().all(|&value| value == 0.0));
    let mut unnormalized = phi.clone();
    for value in &mut unnormalized[0] {
        *value *= 2.0;
    }
    assert!(
        ModalPlateViscosity::from_cst_dkt_shell(
            &mesh,
            &viscosity,
            &model,
            &unnormalized,
            limits(),
            &gate
        )
        .is_err()
    );
    let mut work = limits();
    work.max_projection_terms = 1;
    assert!(
        ModalPlateViscosity::from_cst_dkt_shell(&mesh, &viscosity, &model, &phi, work, &gate)
            .is_err()
    );
    let mut count = limits();
    count.max_modes = modes.len() - 1;
    assert!(
        ModalPlateViscosity::from_cst_dkt_shell(&mesh, &viscosity, &model, &phi, count, &gate)
            .is_err()
    );
    gate.request();
    assert!(matches!(
        ModalPlateViscosity::from_cst_dkt_shell(&mesh, &viscosity, &model, &phi, limits(), &gate),
        Err(AcousticRealizeError::Cancelled)
    ));
}
