use super::*;
use crate::plate_viscosity::shell::tests::fixture;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::{PortHamiltonian, QuadraticStorage};
use fs_plate::shell::CstDktMass;
use fs_plate::shell::viscosity::{ShellViscosity, assemble_constitutive_operator};

mod coupled;

fn limits() -> PlateRelaxationLimits {
    PlateRelaxationLimits {
        max_modes: 12,
        max_projection_terms: 10_000_000,
        max_memory_states: 64,
    }
}

fn spectra(viscosity: &[ShellViscosity]) -> Vec<Vec<ShellRelaxationBranch>> {
    viscosity
        .iter()
        .enumerate()
        .map(|(element, law)| {
            let (a, b, d) = law.matrices();
            [
                (100.0, 0.0002),
                (40.0, if element == 0 { 0.004 } else { 0.02 }),
            ]
            .into_iter()
            .map(|(scale, time)| {
                ShellRelaxationBranch::new_coupled(
                    a.map(|x| scale * x),
                    b.map(|x| scale * x),
                    d.map(|x| scale * x),
                    time,
                )
                .unwrap()
            })
            .collect()
        })
        .collect()
}

fn near(actual: f64, expected: f64, relative: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= relative * scale.max(1e-20),
        "{actual:e} != {expected:e}, scale={scale:e}"
    );
}

#[test]
fn every_projected_time_preserves_original_spatial_and_channel_storage() {
    for mass in [CstDktMass::Lumped, CstDktMass::Consistent] {
        let (mesh, model, pairs, viscosity) = fixture(mass);
        let phi: Vec<_> = pairs.iter().map(|mode| mode.phi.clone()).collect();
        let source = spectra(&viscosity);
        let gate = CancelGate::new_clock_free();
        let image =
            ModalShellRelaxation::from_cst_dkt_shell(&mesh, &source, &model, &phi, limits(), &gate)
                .unwrap();
        assert_eq!(image.projected_branches().len(), 3);
        assert_eq!(image.mode_count(), phi.len());
        assert_eq!(
            image.coupling_rows_per_s().len(),
            image.relaxation_times_s().len()
        );
        let q = [0.2, -0.1, 0.07, -0.03, 0.04, -0.02];
        let mut physical = vec![0.0; model.free];
        for (mode, &value) in phi.iter().zip(&q) {
            for (out, &shape) in physical.iter_mut().zip(mode) {
                *out += shape * value;
            }
        }
        for branch in image.projected_branches() {
            let sections: Vec<_> = source
                .iter()
                .map(|section| {
                    let mut a = [0.0; 9];
                    let mut b = a;
                    let mut d = a;
                    for arm in section
                        .iter()
                        .filter(|arm| arm.relaxation_time_s() == branch.relaxation_time_s())
                    {
                        let (aa, bb, dd) = arm.matrices();
                        for i in 0..9 {
                            a[i] += aa[i];
                            b[i] += bb[i];
                            d[i] += dd[i];
                        }
                    }
                    ShellViscosity::new_coupled(a, b, d).unwrap()
                })
                .collect();
            let spatial = assemble_constitutive_operator(
                &mesh,
                &sections,
                &model,
                ShellViscosityBudget {
                    max_elements: 2,
                    max_full_dofs: 24,
                    max_dofs: 24,
                    max_modes: 12,
                    max_work: 10_000_000,
                },
                || false,
            )
            .unwrap();
            let mut force = vec![0.0; model.free];
            spatial.spmv(&physical, &mut force);
            let expected: f64 = 0.5 * physical.iter().zip(force).map(|(u, f)| u * f).sum::<f64>();
            let channels = 0.5
                * image
                    .coupling_rows_per_s()
                    .iter()
                    .zip(image.relaxation_times_s())
                    .filter(|(_, time)| **time == branch.relaxation_time_s())
                    .map(|(row, _)| {
                        row.iter()
                            .zip(q)
                            .map(|(c, value)| c * value)
                            .sum::<f64>()
                            .powi(2)
                    })
                    .sum::<f64>();
            let modal = 0.5
                * (0..q.len())
                    .map(|i| {
                        (0..q.len())
                            .map(|j| q[i] * branch.matrix()[q.len() * i + j] * q[j])
                            .sum::<f64>()
                    })
                    .sum::<f64>();
            near(channels, expected, 3e-9, expected);
            near(modal, expected, 3e-9, expected);
            assert!(
                branch
                    .matrix()
                    .iter()
                    .enumerate()
                    .any(|(i, x)| i / q.len() != i % q.len() && x.abs() > 1e-6)
            );
        }
        let embedded = image.branches(1, phi.len() + 2, &gate).unwrap();
        for (arm, row) in embedded.iter().zip(image.coupling_rows_per_s()) {
            assert_eq!(arm.projection[0], 0.0);
            assert_eq!(arm.projection[2 * (phi.len() + 1)], 0.0);
            assert!(arm.projection.iter().skip(1).step_by(2).all(|x| *x == 0.0));
            for (i, &value) in row.iter().enumerate() {
                assert_eq!(arm.projection[2 * (1 + i)], value);
            }
        }
    }
}

fn phs(omega2: &[f64], image: &ModalShellRelaxation, gate: &CancelGate) -> PortHamiltonian {
    let n = omega2.len();
    let dimension = 2 * n;
    let mut storage = vec![0.0; dimension * dimension];
    let mut skew = storage.clone();
    let mut ports = vec![0.0; dimension * n];
    for i in 0..n {
        storage[dimension * (2 * i) + 2 * i] = omega2[i];
        storage[dimension * (2 * i + 1) + 2 * i + 1] = 1.0;
        skew[dimension * (2 * i) + 2 * i + 1] = 1.0;
        skew[dimension * (2 * i + 1) + 2 * i] = -1.0;
        ports[n * (2 * i + 1) + i] = 1.0;
    }
    PortHamiltonian::new(
        dimension,
        n,
        skew,
        vec![0.0; dimension * dimension],
        ports,
        Box::new(QuadraticStorage::new(storage, dimension).unwrap()),
    )
    .unwrap()
    .with_relaxation_branches(image.branches(0, n, gate).unwrap())
    .unwrap()
}

// Independent continuous equations, including direct irreversible power and
// input work. They never infer heat from the difference of stored energies.
fn derivative(x: &[f64], omega2: &[f64], image: &ModalShellRelaxation, force: &[f64]) -> Vec<f64> {
    let n = omega2.len();
    let states = 2 * n + image.coupling_rows_per_s().len();
    let mut out = vec![0.0; states + 2];
    for i in 0..n {
        out[2 * i] = x[2 * i + 1];
        out[2 * i + 1] = force[i] - omega2[i] * x[2 * i];
        out[states + 1] += force[i] * x[2 * i + 1];
    }
    for (arm, (row, &time)) in image
        .coupling_rows_per_s()
        .iter()
        .zip(image.relaxation_times_s())
        .enumerate()
    {
        let extension = row
            .iter()
            .enumerate()
            .map(|(i, c)| c * x[2 * i])
            .sum::<f64>()
            - x[2 * n + arm];
        for i in 0..n {
            out[2 * i + 1] -= row[i] * extension;
        }
        out[2 * n + arm] = extension / time;
        out[states] += extension * extension / time;
    }
    out
}

fn rk4(x: &mut [f64], omega2: &[f64], image: &ModalShellRelaxation, force: &[f64], dt: f64) {
    let stage = |rate: &[f64], scale: f64| {
        x.iter()
            .zip(rate)
            .map(|(a, b)| a + scale * b)
            .collect::<Vec<_>>()
    };
    let k1 = derivative(x, omega2, image, force);
    let k2 = derivative(&stage(&k1, 0.5 * dt), omega2, image, force);
    let k3 = derivative(&stage(&k2, 0.5 * dt), omega2, image, force);
    let k4 = derivative(&stage(&k3, dt), omega2, image, force);
    for i in 0..x.len() {
        x[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

#[test]
fn hereditary_shell_branches_drive_existing_exact_storage_work_and_heat() {
    let (mesh, model, pairs, viscosity) = fixture(CstDktMass::Consistent);
    let phi: Vec<_> = pairs.iter().map(|mode| mode.phi.clone()).collect();
    let omega2: Vec<_> = pairs.iter().map(|mode| mode.lambda).collect();
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
    let system = phs(&omega2, &image, &gate);
    let n = image.mode_count();
    let states = system.state_dim();
    let dt = 1e-4;
    let prepared = PreparedExactZoh::new(
        &system,
        dt,
        LinearSplitBudget {
            max_modes: states,
            max_omega_dt: 0.25,
        },
    )
    .unwrap();
    let forces = [0.1, -0.03, 0.04, -0.02, 0.01, 0.015];
    for relaxed in [false, true] {
        let mut initial = vec![0.0; states];
        for i in 0..n {
            initial[2 * i] = (i as f64 + 1.0) * 1e-5;
            initial[2 * i + 1] = if i % 2 == 0 { 0.01 } else { -0.02 };
        }
        if relaxed {
            for (arm, row) in image.coupling_rows_per_s().iter().enumerate() {
                initial[2 * n + arm] = row
                    .iter()
                    .enumerate()
                    .map(|(i, c)| c * initial[2 * i])
                    .sum();
            }
        }
        let mechanical: f64 = (0..n)
            .map(|i| 0.5 * (omega2[i] * initial[2 * i].powi(2) + initial[2 * i + 1].powi(2)))
            .sum();
        let start_energy = system.storage().hamiltonian(&initial);
        if relaxed {
            near(start_energy, mechanical, 1e-13, mechanical);
        } else {
            assert!(start_energy > mechanical);
        }
        let mut exact = initial.clone();
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..64 {
            let frame = prepared.step(&exact, &forces).unwrap();
            assert!(frame.dissipated >= 0.0);
            assert!(frame.balance_residual().abs() <= 2e-9 * start_energy.max(1e-12));
            heat += frame.dissipated;
            work += frame.supplied;
            exact = frame.x;
        }
        assert!(heat > 0.0);
        near(
            system.storage().hamiltonian(&exact) + heat - work,
            start_energy,
            3e-9,
            start_energy,
        );
        let mut errors = Vec::new();
        for subdivisions in [16, 32, 64] {
            let mut reference = initial.clone();
            reference.extend([0.0, 0.0]);
            for _ in 0..64 * subdivisions {
                rk4(
                    &mut reference,
                    &omega2,
                    &image,
                    &forces,
                    dt / subdivisions as f64,
                );
            }
            let error = exact
                .iter()
                .zip(&reference)
                .enumerate()
                .map(|(i, (a, b))| {
                    let weight = if i < 2 * n && i % 2 == 0 {
                        omega2[i / 2].sqrt()
                    } else {
                        1.0
                    };
                    (weight * (a - b)).powi(2)
                })
                .sum::<f64>()
                .sqrt();
            errors.push(error);
            if subdivisions == 64 {
                near(reference[states], heat, 2e-7, heat);
                near(reference[states + 1], work, 2e-7, work.abs());
                assert!(error <= 2e-7 * (2.0 * start_energy).sqrt());
            }
        }
        assert!(errors[1] < 0.15 * errors[0], "RK4 errors {errors:?}");
        assert!(errors[2] < 0.15 * errors[1], "RK4 errors {errors:?}");
    }
}

#[test]
fn zero_memory_and_bad_basis_state_work_embedding_and_cancellation_are_explicit() {
    let (mesh, model, pairs, viscosity) = fixture(CstDktMass::Consistent);
    let phi: Vec<_> = pairs.iter().map(|mode| mode.phi.clone()).collect();
    let source = spectra(&viscosity);
    let gate = CancelGate::new_clock_free();
    let mut no_memory = limits();
    no_memory.max_memory_states = 0;
    let empty = ModalShellRelaxation::from_cst_dkt_shell(
        &mesh,
        &[vec![], vec![]],
        &model,
        &phi,
        no_memory,
        &gate,
    )
    .unwrap();
    assert!(empty.coupling_rows_per_s().is_empty());
    assert!(empty.projected_branches().is_empty());
    assert!(
        ModalShellRelaxation::from_cst_dkt_shell(&mesh, &source, &model, &phi, no_memory, &gate)
            .is_err()
    );
    let mut work = limits();
    work.max_projection_terms = 1;
    assert!(
        ModalShellRelaxation::from_cst_dkt_shell(&mesh, &source, &model, &phi, work, &gate)
            .is_err()
    );
    let mut wrong = phi.clone();
    wrong[0].iter_mut().for_each(|value| *value *= 2.0);
    assert!(
        ModalShellRelaxation::from_cst_dkt_shell(&mesh, &source, &model, &wrong, limits(), &gate)
            .is_err()
    );
    let image =
        ModalShellRelaxation::from_cst_dkt_shell(&mesh, &source, &model, &phi, limits(), &gate)
            .unwrap();
    assert!(image.branches(1, image.mode_count(), &gate).is_err());
    assert!(image.branches(0, usize::MAX, &gate).is_err());
    gate.request();
    assert!(matches!(
        image.branches(0, image.mode_count(), &gate),
        Err(AcousticRealizeError::Cancelled)
    ));
    assert!(matches!(
        ModalShellRelaxation::from_cst_dkt_shell(&mesh, &source, &model, &phi, limits(), &gate),
        Err(AcousticRealizeError::Cancelled)
    ));
}
