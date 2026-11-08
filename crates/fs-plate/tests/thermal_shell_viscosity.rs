//! G0/G1/G3: full spatial shell damping, independent director decay and heat.
//! Synthetic resultants check numerical laws, not a measured material spectrum.

use fs_material::sheet::SheetMatrix;
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::relaxation::ShearShellRelaxationBranch;
use fs_plate::shell::mitc::thermal::transient::{EnrichedMitcShell, MitcRelaxationBudget};
use fs_plate::shell::mitc::viscosity::{ShearShellViscosity, project_viscosity};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_sparse::Csr;
use std::sync::Arc;

const WORK: usize = 1_000_000;
const I2: f64 = 0.04;

fn section(scale: f64) -> ShearShellSection {
    ShearShellSection::new_coupled(
        SheetMatrix::orthotropic(2e4 * scale, 1e4 * scale, 0.2 * scale, 3e3 * scale).unwrap(),
        [0.0; 9],
        SheetMatrix::orthotropic(20.0 * scale, 10.0 * scale, 0.2 * scale, 3.0 * scale).unwrap(),
        [120.0 * scale, 0.0, 80.0 * scale],
        1.0,
        2.0,
        0.1,
        I2,
        0.0,
    )
    .unwrap()
}

fn geometry() -> (ShellMesh, NodalDirectors) {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    (mesh, directors)
}

fn model(interior_only: bool, scale: f64) -> Arc<MitcShellModel> {
    let (mesh, directors) = geometry();
    let constraints: Vec<_> = if interior_only {
        (0..3)
            .map(|node| ShellConstraint {
                node,
                // One numerical drilling coordinate remains free; it has no
                // physical inertia or material loss. Only the two directors move.
                mask: [true, true, true, true, true, node != 0],
                frame: ConstraintFrame::Global,
            })
            .collect()
    } else {
        vec![]
    };
    Arc::new(
        assemble_mitc3p(
            &mesh,
            &directors,
            &[section(scale)],
            &constraints,
            ShellMass::Consistent,
            MitcBudget {
                max_nodes: 3,
                max_triangles: 1,
            },
            1e-3,
            || false,
        )
        .unwrap(),
    )
}

fn entry(matrix: &Csr, row: usize, column: usize) -> f64 {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .zip(values)
        .find_map(|(&c, &v)| (c == column).then_some(v))
        .unwrap_or(0.0)
}

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}

fn proportional(beta: f64) -> ShearShellViscosity {
    let a = [2e4, 0.2, 0.0, 0.2, 1e4, 0.0, 0.0, 0.0, 3e3];
    let d = [20.0, 0.2, 0.0, 0.2, 10.0, 0.0, 0.0, 0.0, 3.0];
    ShearShellViscosity::new_coupled(
        a.map(|x| beta * x),
        [0.0; 9],
        d.map(|x| beta * x),
        [beta * 120.0, 0.0, beta * 80.0],
    )
    .unwrap()
}

fn coupled_loss() -> ShearShellViscosity {
    ShearShellViscosity::new_coupled(
        [2.0, 0.2, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0, 0.3],
        [0.05, 0.0, 0.0, 0.0, -0.02, 0.0, 0.0, 0.0, 0.01],
        [0.05, 0.004, 0.0, 0.004, 0.03, 0.0, 0.0, 0.0, 0.02],
        [0.8, 0.1, 0.6],
    )
    .unwrap()
}

fn memory(law: ShearShellViscosity, tau: f64) -> ShearShellRelaxationBranch {
    let (a, b, d, s) = law.coupled_matrices();
    ShearShellRelaxationBranch::new_coupled(a, b, d, [s[0][0], s[0][1], s[1][1]], tau).unwrap()
}

fn memory_budget(branches: usize) -> MitcRelaxationBudget {
    MitcRelaxationBudget {
        max_element_branches: branches,
        max_assembly_terms: branches * WORK,
        max_payload_bytes: branches * 5000,
    }
}

// G1 independent continuous Maxwell/Kelvin oscillator, including both heat
// channels as differential states. No production Prony coefficients are used.
fn memory_reference(omega: f64, mass: f64, beta: f64, tau: [f64; 2], duration: f64) -> [f64; 6] {
    let alpha = [0.4, 0.7];
    let stiffness = mass * omega * omega;
    let derivative = |x: [f64; 6]| {
        [
            x[1],
            -omega * omega * (x[0] + beta * x[1] + alpha[0] * x[2] + alpha[1] * x[3]),
            x[1] - x[2] / tau[0],
            x[1] - x[3] / tau[1],
            stiffness * (alpha[0] * x[2] * x[2] / tau[0] + alpha[1] * x[3] * x[3] / tau[1]),
            beta * stiffness * x[1] * x[1],
        ]
    };
    let dt = duration / 20000.0;
    let mut x = [0.0, 0.2, 0.0, 0.0, 0.0, 0.0];
    for _ in 0..20000 {
        let a = derivative(x);
        let b = derivative(core::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| x[i] + dt * c[i]));
        x = core::array::from_fn(|i| x[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}

#[test]
fn g1_original_director_maxwell_memory_matches_continuous_motion_and_separate_heat() {
    let original = EnrichedMitcShell::try_new(model(true, 1.0), || false).unwrap();
    let n = original.model().pencil().free;
    let mass = 0.5 * I2 * 81.0 / 280.0;
    let a = entry(original.stiffness(), n, n);
    let b = entry(original.stiffness(), n, n + 1);
    let d = entry(original.stiffness(), n + 1, n + 1);
    let lambda = 0.5 * (a + d) + (0.5 * (a - d)).hypot(b);
    let vector = if b.abs() > 1e-14 {
        let length = b.hypot(lambda - a);
        [b / length, (lambda - a) / length]
    } else if a >= d {
        [1.0, 0.0]
    } else {
        [0.0, 1.0]
    };
    let omega = (lambda / mass).sqrt();
    let (beta, tau, duration) = (0.25 / omega, [0.3 / omega, 1.7 / omega], 3.0 / omega);
    let shell = original
        .with_viscosity(&[proportional(beta)], WORK, || false)
        .unwrap()
        .with_relaxation(
            &[vec![
                memory(proportional(0.4), tau[0]),
                memory(proportional(0.7), tau[1]),
            ]],
            memory_budget(2),
            || false,
        )
        .unwrap();
    let reference = memory_reference(omega, mass, beta, tau, duration);
    let mut errors = Vec::new();
    for steps in [40, 80, 160] {
        let solver = shell
            .prepare_midpoint(duration / f64::from(steps), || false)
            .unwrap();
        let initial = shell
            .state(
                0.0,
                &[0.0; 18],
                &[0.0; 18],
                &[[0.0; 2]],
                &[vector.map(|x| 0.2 * x)],
                || false,
            )
            .unwrap();
        let mut state = shell
            .with_relaxation_history(&initial, &[[0.0; 20]; 2], || false)
            .unwrap();
        let energy = state.kinetic_energy_j();
        let mut heat = [0.0; 2];
        let force = vec![0.0; shell.dimension()];
        for _ in 0..steps {
            let next = solver.advance(&state, &force, &force, || false).unwrap();
            assert!(next.report.relaxation_dissipated_energy_j >= 0.0);
            assert!(next.report.viscous_dissipated_energy_j >= 0.0);
            near(next.report.energy_residual_j, 0.0, 0.0, 3e-11 * energy);
            heat[0] += next.report.relaxation_dissipated_energy_j;
            heat[1] += next.report.viscous_dissipated_energy_j;
            state = next.state;
        }
        let amplitude = |q: &[f64; 2]| q[0] * vector[0] + q[1] * vector[1];
        let q = amplitude(&state.bubble_coordinates()[0]);
        let v = amplitude(&state.bubble_velocity()[0]);
        let history = state.relaxation_history().unwrap();
        let actual = [
            q,
            v,
            amplitude(&[history[0][18], history[0][19]]),
            amplitude(&[history[1][18], history[1][19]]),
            heat[0],
            heat[1],
        ];
        let scale = [
            omega / 0.2,
            5.0,
            omega / 0.2,
            omega / 0.2,
            1.0 / energy,
            1.0 / energy,
        ];
        errors.push(
            actual
                .iter()
                .zip(reference)
                .zip(scale)
                .map(|((a, b), scale)| (a - b).abs() * scale)
                .fold(0.0_f64, f64::max),
        );
        near(
            state.kinetic_energy_j()
                + 0.5 * lambda * q * q
                + state.relaxation_energy_j()
                + heat[0]
                + heat[1],
            energy,
            2e-10,
            1e-18,
        );
        assert!(heat[0] > 0.0 && heat[1] > 0.0);
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
    assert!(errors[2] < 3e-4, "{errors:?}");
}

#[test]
fn g3_coupled_spatial_memory_retains_null_motion_and_transactional_histories() {
    let original = EnrichedMitcShell::try_new(model(false, 1.0), || false).unwrap();
    let shell = original
        .with_viscosity(&[coupled_loss()], WORK, || false)
        .unwrap()
        .with_relaxation(
            &[vec![
                memory(coupled_loss(), 1e-4),
                memory(proportional(0.3), 0.002),
            ]],
            memory_budget(2),
            || false,
        )
        .unwrap();
    let solver = shell.prepare_midpoint(1e-4, || false).unwrap();
    let mut velocity = [0.0; 18];
    velocity[6] = 0.4;
    velocity[13] = -0.2;
    let initial = shell
        .state(
            0.0,
            &[0.0; 18],
            &velocity,
            &[[0.0; 2]],
            &[[0.3, -0.2]],
            || false,
        )
        .unwrap();
    let mut history = [[0.0; 20]; 2];
    history[0][6] = 0.005;
    history[1][18] = 0.015;
    history[1][19] = -0.02;
    let initial = shell
        .with_relaxation_history(&initial, &history, || false)
        .unwrap();
    let energy = initial.kinetic_energy_j() + initial.relaxation_energy_j();
    let force = vec![0.0; shell.dimension()];
    let mut polls = 0;
    let next = solver
        .advance(&initial, &force, &force, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(next.report.relaxation_dissipated_energy_j > 0.0);
    assert!(next.report.viscous_dissipated_energy_j > 0.0);
    near(next.report.energy_residual_j, 0.0, 0.0, 1e-9 * energy);
    near(
        next.facet_dissipated_energy_j[0],
        next.facet_viscous_dissipated_energy_j[0] + next.facet_relaxation_dissipated_energy_j[0],
        1e-14,
        0.0,
    );
    let mut calls = 0;
    assert!(matches!(
        solver.advance(&initial, &force, &force, || {
            calls += 1;
            calls >= polls - 1
        }),
        Err(PlateError::Cancelled)
    ));
    assert_eq!(initial.relaxation_history().unwrap(), &history);
    assert_eq!(initial.time_s(), 0.0);
    let mut translation = [0.0; 18];
    for node in 0..3 {
        translation[6 * node] = 0.4;
        translation[6 * node + 1] = -0.3;
    }
    let rigid = shell
        .state(
            0.0,
            &[0.0; 18],
            &translation,
            &[[0.0; 2]],
            &[[0.0; 2]],
            || false,
        )
        .unwrap();
    let rigid = shell
        .with_relaxation_history(&rigid, &[[0.0; 20]; 2], || false)
        .unwrap();
    let moved = solver.advance(&rigid, &force, &force, || false).unwrap();
    assert!(moved.report.dissipated_energy_j < 1e-25);
    assert!(moved.state.relaxation_energy_j() < 1e-25);
    for (&a, &b) in moved.state.nodal_velocity_6n().iter().zip(&translation) {
        near(a, b, 1e-10, 1e-11);
    }
}

#[test]
fn g0_memory_initialization_bounds_and_material_refresh_are_explicit() {
    let original = EnrichedMitcShell::try_new(model(true, 1.0), || false).unwrap();
    let spectrum = [vec![memory(proportional(0.4), 0.01)]];
    for budget in [
        MitcRelaxationBudget {
            max_element_branches: 0,
            ..memory_budget(1)
        },
        MitcRelaxationBudget {
            max_assembly_terms: WORK - 1,
            ..memory_budget(1)
        },
        MitcRelaxationBudget {
            max_payload_bytes: 0,
            ..memory_budget(1)
        },
    ] {
        assert!(
            original
                .with_relaxation(&spectrum, budget, || false)
                .is_err()
        );
    }
    let shell = original
        .with_relaxation(&spectrum, memory_budget(1), || false)
        .unwrap();
    let initial = shell
        .state(
            0.0,
            &[0.0; 18],
            &[0.0; 18],
            &[[0.01, -0.02]],
            &[[0.1, 0.2]],
            || false,
        )
        .unwrap();
    let force = vec![0.0; shell.dimension()];
    assert!(
        shell
            .prepare_midpoint(1e-4, || false)
            .unwrap()
            .advance(&initial, &force, &force, || false)
            .is_err()
    );
    assert!(
        shell
            .with_relaxation_history(&initial, &[], || false)
            .is_err()
    );
    let initial = shell
        .with_relaxation_history(&initial, &[[0.0; 20]], || false)
        .unwrap();
    let refreshed = shell
        .with_material_model(model(true, 1.3), || false)
        .unwrap();
    let solver = refreshed.prepare_midpoint(1e-4, || false).unwrap();
    assert!(
        solver
            .advance_from(&original, &initial, &force, &force, || false)
            .is_err()
    );
    let next = solver
        .advance_from(&shell, &initial, &force, &force, || false)
        .unwrap();
    assert!(next.report.relaxation_dissipated_energy_j > 0.0);
    assert!(next.report.energy_residual_j.abs() < 1e-10);
    assert_eq!(initial.relaxation_history().unwrap(), &[[0.0; 20]]);
}

#[test]
fn g1_complete_spatial_operator_retains_existing_nodal_projection_and_cross_loss() {
    let original = model(false, 1.0);
    let shell = EnrichedMitcShell::try_new(Arc::clone(&original), || false).unwrap();
    let law = coupled_loss();
    let damped = shell.with_viscosity(&[law], WORK, || false).unwrap();
    let n = original.pencil().free;
    let modes: Vec<_> = (0..n)
        .map(|i| {
            let mut coordinate = vec![0.0; n];
            coordinate[i] = 1.0;
            original.expand(&coordinate).unwrap()
        })
        .collect();
    let (mesh, directors) = geometry();
    let projected = project_viscosity(
        &mesh,
        &directors,
        &[section(1.0)],
        &[law],
        &modes,
        MitcBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        2_000_000,
        || false,
    )
    .unwrap();
    let c = damped.material_viscosity().unwrap();
    for i in 0..n {
        for j in 0..n {
            near(entry(c, i, j), projected[n * i + j], 2e-12, 2e-14);
        }
    }
    assert!((0..n).any(|i| entry(c, i, n).abs() > 1e-6));
    assert!(entry(c, n, n) > 0.0 && entry(c, n + 1, n + 1) > 0.0);
    // Physical director-axis rotations are absent from every material strain.
    for row in [5, 11, 17] {
        for column in 0..damped.dimension() {
            assert_eq!(entry(c, row, column), 0.0);
        }
    }
}

#[test]
fn g1_internal_director_decay_converges_to_independent_damped_oscillator() {
    let shell = EnrichedMitcShell::try_new(model(true, 1.0), || false).unwrap();
    let n = shell.model().pencil().free;
    assert_eq!(n, 1);
    // Independent closed 2x2 symmetric eigenpair; physical bubble inertia is
    // area*I2*integral(b^2) = (1/2)*I2*(81/280), with b=27*l0*l1*l2.
    let mass = 0.5 * I2 * 81.0 / 280.0;
    let a = entry(shell.stiffness(), n, n);
    let b = entry(shell.stiffness(), n, n + 1);
    let d = entry(shell.stiffness(), n + 1, n + 1);
    let lambda = 0.5 * (a + d) + (0.5 * (a - d)).hypot(b);
    let vector = if b.abs() > 1e-14 {
        let length = b.hypot(lambda - a);
        [b / length, (lambda - a) / length]
    } else if a >= d {
        [1.0, 0.0]
    } else {
        [0.0, 1.0]
    };
    let omega = (lambda / mass).sqrt();
    let beta = 0.4 / omega;
    let sigma = 0.5 * beta * omega * omega;
    let wd = (omega * omega - sigma * sigma).sqrt();
    let duration = 3.0 / omega;
    let initial_rate = 0.2;
    let envelope = (-sigma * duration).exp();
    let exact_q = initial_rate * envelope * (wd * duration).sin() / wd;
    let exact_v =
        initial_rate * envelope * ((wd * duration).cos() - sigma / wd * (wd * duration).sin());
    let damped = shell
        .with_viscosity(&[proportional(beta)], WORK, || false)
        .unwrap();
    let mut errors = Vec::new();
    for steps in [40, 80, 160] {
        let dt = duration / f64::from(steps);
        let solver = damped.prepare_midpoint(dt, || false).unwrap();
        let mut state = damped
            .state(
                0.0,
                &[0.0; 18],
                &[0.0; 18],
                &[[0.0; 2]],
                &[vector.map(|x| initial_rate * x)],
                || false,
            )
            .unwrap();
        let initial_energy = state.kinetic_energy_j();
        let mut heat = 0.0;
        let zero = vec![0.0; damped.dimension()];
        for _ in 0..steps {
            let next = solver.advance(&state, &zero, &zero, || false).unwrap();
            assert!(next.report.dissipated_energy_j >= 0.0);
            assert_eq!(
                next.facet_dissipated_energy_j,
                vec![next.report.dissipated_energy_j]
            );
            assert!(next.report.energy_residual_j.abs() < 2e-11 * initial_energy);
            heat += next.report.dissipated_energy_j;
            state = next.state;
        }
        let q = state.bubble_coordinates()[0];
        let v = state.bubble_velocity()[0];
        let amplitude = q[0] * vector[0] + q[1] * vector[1];
        let rate = v[0] * vector[0] + v[1] * vector[1];
        let energy = 0.5 * mass * (rate * rate + omega * omega * amplitude * amplitude);
        near(heat + energy, initial_energy, 2e-10, 1e-18);
        errors.push((omega * (amplitude - exact_q)).hypot(rate - exact_v) / initial_rate);
        assert!(heat > 0.0);
    }
    assert!(
        errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9,
        "{errors:?}"
    );
    assert!(errors[2] < 8e-5, "{errors:?}");
}

#[test]
fn g3_nonproportional_heat_closes_work_and_preserves_rigid_translation() {
    let shell = EnrichedMitcShell::try_new(model(false, 1.0), || false).unwrap();
    let shell = shell
        .with_viscosity(&[coupled_loss()], WORK, || false)
        .unwrap();
    let solver = shell.prepare_midpoint(1e-4, || false).unwrap();
    let mut velocity = [0.0; 18];
    velocity[6] = 0.4;
    velocity[13] = -0.2;
    let initial = shell
        .state(
            0.0,
            &[0.0; 18],
            &velocity,
            &[[0.0; 2]],
            &[[0.3, -0.2]],
            || false,
        )
        .unwrap();
    let force = vec![0.0; shell.dimension()];
    let next = solver.advance(&initial, &force, &force, || false).unwrap();
    assert!(next.report.dissipated_energy_j > 0.0);
    assert!(next.report.energy_residual_j.abs() < 1e-10 * initial.kinetic_energy_j());
    let initial_coordinates = initial.coordinates().to_vec();
    let mut translation = [0.0; 18];
    for node in 0..3 {
        translation[6 * node] = 0.4;
        translation[6 * node + 1] = -0.3;
    }
    let rigid = shell
        .state(
            0.0,
            &[0.0; 18],
            &translation,
            &[[0.0; 2]],
            &[[0.0; 2]],
            || false,
        )
        .unwrap();
    let moved = solver.advance(&rigid, &force, &force, || false).unwrap();
    assert!(moved.report.dissipated_energy_j < 1e-25);
    for (&actual, &expected) in moved.state.nodal_velocity_6n().iter().zip(&translation) {
        near(actual, expected, 1e-10, 1e-11);
    }
    assert_eq!(initial.coordinates(), initial_coordinates);
}

#[test]
fn g0_lossless_parity_material_refresh_and_transactional_refusals() {
    let shell = EnrichedMitcShell::try_new(model(true, 1.0), || false).unwrap();
    let zero_law = ShearShellViscosity::new([0.0; 9], [0.0; 9], [0.0; 3]).unwrap();
    let zero = shell.with_viscosity(&[zero_law], WORK, || false).unwrap();
    let initial = shell
        .state(
            0.0,
            &[0.0; 18],
            &[0.0; 18],
            &[[0.01, -0.02]],
            &[[0.1, 0.2]],
            || false,
        )
        .unwrap();
    let force = vec![0.0; shell.dimension()];
    let a = shell
        .prepare_midpoint(1e-4, || false)
        .unwrap()
        .advance(&initial, &force, &force, || false)
        .unwrap();
    let b = zero
        .prepare_midpoint(1e-4, || false)
        .unwrap()
        .advance(&initial, &force, &force, || false)
        .unwrap();
    assert_eq!(a.state.coordinates(), b.state.coordinates());
    assert_eq!(a.state.velocity(), b.state.velocity());
    assert_eq!(a.report, b.report);
    assert_eq!(a.report.dissipated_energy_j, 0.0);
    assert!(a.facet_dissipated_energy_j.is_empty());
    assert!(shell.with_viscosity(&[], WORK, || false).is_err());
    assert!(
        shell
            .with_viscosity(&[coupled_loss()], WORK - 1, || false)
            .is_err()
    );
    assert!(matches!(
        shell.with_viscosity(&[coupled_loss()], WORK, || true),
        Err(PlateError::Cancelled)
    ));
    let damped = shell
        .with_viscosity(&[coupled_loss()], WORK, || false)
        .unwrap();
    let refreshed = damped
        .with_material_model(model(true, 1.3), || false)
        .unwrap();
    assert!(std::ptr::eq(
        damped.material_viscosity().unwrap(),
        refreshed.material_viscosity().unwrap()
    ));
    let stepper = refreshed.prepare_midpoint(1e-4, || false).unwrap();
    assert!(
        stepper
            .advance_from(&shell, &initial, &force, &force, || false)
            .is_err()
    );
    let next = stepper
        .advance_from(&damped, &initial, &force, &force, || false)
        .unwrap();
    assert!(next.report.dissipated_energy_j > 0.0);
    assert!(next.report.energy_residual_j.abs() < 1e-10);
    assert!(matches!(
        stepper.advance_from(&damped, &initial, &force, &force, || true),
        Err(PlateError::Cancelled)
    ));
    assert_eq!(initial.bubble_coordinates(), &[[0.01, -0.02]]);
}
