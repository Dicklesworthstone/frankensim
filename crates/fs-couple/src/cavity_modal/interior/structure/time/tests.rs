use super::*;
use super::super::{NormalDisplacementTrace, StructuralMatrices};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousHelmholtzFe;
use fs_feec::acoustics::helmholtz::heterogeneous::time::ResistiveBoundary;
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_feec::fixtures::single_tet;
use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicResponseControls;
use fs_sparse::{Coo, Csr};

#[path = "source_tests.rs"]
mod source_tests;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 83, kernel_id: 21, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 128,
        maximum_payload_bytes: 16 * 1024 * 1024,
        maximum_relative_residual: 1e-10,
    }
}

fn dense(values: &[f64], n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            if values[i * n + j] != 0.0 { coo.push(i, j, values[i * n + j]); }
        }
    }
    coo.assemble()
}

fn uniform<'f, 'c>(geometry: &'f HelmholtzFe<'c>, medium: Medium)
    -> HeterogeneousHelmholtzFe<'f, 'c>
{
    HeterogeneousHelmholtzFe::new(
        geometry, geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(medium), &mut || false,
    ).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual {actual:.16e}, expected {expected:.16e}");
}

#[test]
fn coupled_compression_trajectory_matches_an_independent_damped_oscillator() {
    let (complex, points) = single_tet();
    let medium = Medium { density_kg_m3: 1.2, sound_speed_m_s: 3.0 };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    // The affine dilation u(x)=q*(x-centroid) has B_i=3*V/4=1/8
    // in this original P1 pressure space. B is proportional to M*1, so
    // constant compression is an EXACT invariant subspace of this FE fixture.
    // This is a semidiscrete analytic check, not a continuum mass-loading claim.
    let traces: Vec<NormalDisplacementTrace> = geometry.boundary_faces().into_iter().map(|face| {
        let (vertices, normal, _) = geometry.boundary_geometry(face).unwrap();
        let shape = vertices.map(|point| (0..3).map(|axis| {
            (point[axis] - 0.25) * normal[axis]
        }).sum());
        (face, shape)
    }).collect();
    let inertia = 2.0;
    let spring = 3.0;
    let drag = 0.4;
    let force = 0.8;
    let duration = 0.23;
    let mass = dense(&[inertia], 1);
    let stiffness = dense(&[spring], 1);
    let damping = dense(&[drag], 1);
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior, StructuralMatrices { mass: &mass, stiffness: &stiffness, damping: Some(&damping) },
            &[traces], &[], controls(), cx,
        ).unwrap();
        for &b in model.coupling() { close(b, 0.125, 1e-14); }
        let prepared = model.prepare_exact_time(duration, StructuralInteriorTimeLimits::default(), cx).unwrap();
        let initial = prepared.state(0.0, &[0.0], &[0.0], &[0.0; 4], &[0.0; 4], cx).unwrap();
        let step = prepared.step(&initial, &[force], &[0.0; 4], cx).unwrap();

        let bulk = medium.density_kg_m3 * medium.sound_speed_m_s * medium.sound_speed_m_s;
        let cavity_spring = 9.0 * (1.0 / 6.0) * bulk;
        let effective = spring + cavity_spring;
        let decay = drag / (2.0 * inertia);
        let frequency = det::sqrt(effective / inertia - decay * decay);
        let phase = frequency * duration;
        let envelope = det::exp(-decay * duration);
        let q = force / effective * (1.0 - envelope * (det::cos(phase) + decay / frequency * det::sin(phase)));
        let v = force / (inertia * frequency) * envelope * det::sin(phase);
        let integral_q = (force * duration - inertia * v - drag * q) / effective;
        let pressure = -3.0 * bulk * q;
        let potential = -3.0 * bulk * integral_q;
        close(step.state.structural_displacement()[0], q, 2e-11);
        close(step.state.structural_velocity()[0], v, 2e-11);
        for &p in step.state.pressure_pa() { close(p, pressure, 3e-11); }
        for &psi in step.state.pressure_potential_pa_s() { close(psi, potential, 3e-11); }
        close(step.frame.fluid_force[0], 0.5 * pressure, 2e-11);
        close(step.frame.integrated_structural_velocity[0], q, 2e-11);
        for &integral in &step.frame.integrated_pressure_pa_s { close(integral, potential, 3e-11); }

        // Independent integral of drag*[F/(m*w)*exp(-a*t)*sin(w*t)]^2.
        let lambda = 2.0 * decay;
        let beta = 2.0 * frequency;
        let exponential = det::exp(-lambda * duration);
        let cosine_integral = (lambda + exponential
            * (-lambda * det::cos(beta * duration) + beta * det::sin(beta * duration)))
            / (lambda * lambda + beta * beta);
        let sine_squared_integral = 0.5 * ((1.0 - exponential) / lambda - cosine_integral);
        let viscous_loss = drag * (force / (inertia * frequency)).powi(2) * sine_squared_integral;
        let exchange = -0.5 * cavity_spring * q * q;
        close(step.frame.work.structural_input_work_j, force * q, 2e-12);
        close(step.frame.work.structural_dissipation_j, viscous_loss, 2e-12);
        close(step.frame.work.fluid_to_structure_work_j, exchange, 2e-12);
        close(step.frame.work.final_energy.fluid_compression_j, -exchange, 2e-12);
        close(step.frame.work.final_energy.structural_kinetic_j, 0.5 * inertia * v * v, 2e-12);
        close(step.frame.work.final_energy.structural_strain_j, 0.5 * spring * q * q, 2e-12);
        assert_eq!(step.frame.work.wall_dissipation_j, 0.0);
        assert_eq!(step.frame.work.acoustic_input_work_j, 0.0);
        assert!(step.frame.work.fluid_to_structure_work_j < 0.0);
        let product_of_means = -cavity_spring * integral_q * q / duration;
        assert!((product_of_means - exchange).abs() > 0.01 * exchange.abs());
        assert!(step.frame.work.balance_residual_j().abs() < 2e-11);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.structural_displacement(), &[0.0]);
    });
}

#[test]
fn free_consistent_mass_and_full_constant_pressure_mode_keep_their_exact_forced_limits() {
    let (complex, points) = single_tet();
    let medium = Medium { density_kg_m3: 1.1, sound_speed_m_s: 2.0 };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[2.0, 0.25, 0.25, 1.2], 2);
    let stiffness = dense(&[0.0; 4], 2);
    let force = [0.8, -0.3];
    let determinant = 2.0 * 1.2 - 0.25 * 0.25;
    let acceleration = [(1.2 * force[0] - 0.25 * force[1]) / determinant,
        (2.0 * force[1] - 0.25 * force[0]) / determinant];
    let displacement = [0.3, -0.2];
    let velocity = [0.7, 0.1];
    let pressure_rate = 0.3;
    let duration = 0.4;
    let flow: Vec<_> = (0..4).map(|i| (0..4).map(|j| {
        pressure_rate * interior.mass().get(i, j)
    }).sum()).collect();
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior, StructuralMatrices { mass: &mass, stiffness: &stiffness, damping: None },
            &[Vec::new(), Vec::new()], &[], controls(), cx,
        ).unwrap();
        let prepared = model.prepare_exact_time(duration, StructuralInteriorTimeLimits::default(), cx).unwrap();
        let mut state = prepared.state(1.25, &displacement, &velocity, &[0.4; 4], &[0.25; 4], cx).unwrap();
        let mut workspace = prepared.workspace();
        let frame = prepared.step_into(&mut state, &force, &flow, &mut workspace, cx).unwrap().clone();
        let expected_increment = (0..2).map(|i| velocity[i] * duration + 0.5 * acceleration[i] * duration * duration).collect::<Vec<_>>();
        for i in 0..2 {
            close(state.structural_displacement()[i], displacement[i] + expected_increment[i], 3e-12);
            close(state.structural_velocity()[i], velocity[i] + acceleration[i] * duration, 3e-12);
        }
        for &p in state.pressure_pa() { close(p, 0.25 + pressure_rate * duration, 3e-12); }
        for &psi in state.pressure_potential_pa_s() {
            close(psi, 0.4 + 0.25 * duration + 0.5 * pressure_rate * duration * duration, 3e-12);
        }
        let structural_work = force.iter().zip(&expected_increment).map(|(f, dq)| f * dq).sum::<f64>();
        let acoustic_work = flow.iter().sum::<f64>() * (0.25 * duration + 0.5 * pressure_rate * duration * duration);
        close(frame.work.structural_input_work_j, structural_work, 3e-12);
        close(frame.work.acoustic_input_work_j, acoustic_work, 3e-12);
        assert_eq!(frame.work.structural_dissipation_j, 0.0);
        assert_eq!(frame.work.wall_dissipation_j, 0.0);
        assert_eq!(frame.work.fluid_to_structure_work_j, 0.0);
        assert_eq!(state.time_s(), 1.25 + duration);

        // A refused force shape cannot advance accepted body, pressure or clock.
        // Reusing the same scratch after refusal reproduces a fresh candidate.
        let saved = state.clone();
        assert!(matches!(prepared.step_into(&mut state, &force[..1], &flow, &mut workspace, cx),
            Err(StructuralInteriorTimeError::Invalid(_))));
        assert_eq!(state.time_s().to_bits(), saved.time_s().to_bits());
        assert_eq!(state.displacement, saved.displacement);
        assert_eq!(state.velocity, saved.velocity);
        let expected = prepared.step(&saved, &force, &flow, cx).unwrap();
        prepared.step_into(&mut state, &force, &flow, &mut workspace, cx).unwrap();
        assert_eq!(state.time_s().to_bits(), expected.state.time_s().to_bits());
        assert_eq!(state.displacement, expected.state.displacement);
        assert_eq!(state.velocity, expected.state.velocity);
    });
}

#[test]
fn original_off_diagonal_inertia_and_gyroscopic_motion_are_retained_without_fake_stiffness() {
    let (complex, points) = single_tet();
    let medium = Medium { density_kg_m3: 1.2, sound_speed_m_s: 2.0 };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[2.0, 0.3, 0.3, 1.1], 2);
    let stiffness = dense(&[0.0; 4], 2);
    let gyro = 0.8;
    let damping = dense(&[0.0, gyro, -gyro, 0.0], 2);
    let duration = 1.3;
    let displacement = [0.2, -0.1];
    let velocity = [0.4, -0.2];
    let determinant = 2.0 * 1.1 - 0.3 * 0.3;
    let omega = gyro / det::sqrt(determinant);
    let av = [
        (0.3 * gyro * velocity[0] - 1.1 * gyro * velocity[1]) / determinant,
        (2.0 * gyro * velocity[0] - 0.3 * gyro * velocity[1]) / determinant,
    ];
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior, StructuralMatrices { mass: &mass, stiffness: &stiffness, damping: Some(&damping) },
            &[Vec::new(), Vec::new()], &[], controls(), cx,
        ).unwrap();
        let prepared = model.prepare_exact_time(duration, StructuralInteriorTimeLimits::default(), cx).unwrap();
        let initial = prepared.state(0.0, &displacement, &velocity, &[0.0; 4], &[0.0; 4], cx).unwrap();
        let step = prepared.step(&initial, &[0.0; 2], &[0.0; 4], cx).unwrap();
        for i in 0..2 {
            let v = det::cos(omega * duration) * velocity[i] + det::sin(omega * duration) / omega * av[i];
            let q = displacement[i] + det::sin(omega * duration) / omega * velocity[i]
                + (1.0 - det::cos(omega * duration)) / (omega * omega) * av[i];
            close(step.state.structural_velocity()[i], v, 2e-11);
            close(step.state.structural_displacement()[i], q, 2e-11);
        }
        close(step.frame.work.final_energy.structural_j(), step.frame.work.initial.structural_j(), 2e-12);
        assert_eq!(step.frame.work.structural_dissipation_j, 0.0);
        assert_eq!(step.frame.work.wall_dissipation_j, 0.0);
        assert_eq!(step.frame.work.structural_input_work_j, 0.0);
        assert_eq!(step.frame.work.acoustic_input_work_j, 0.0);
        assert!(step.state.pressure_pa().iter().all(|p| p.abs() < 1e-14));
    });
}

#[test]
fn original_heterogeneous_p2_fields_and_separate_work_integrals_compose_on_one_clock() {
    let (complex, points) = single_tet();
    let medium = Medium { density_kg_m3: 1.2, sound_speed_m_s: 3.0 };
    let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
    let interior = HeterogeneousHelmholtzFe::new(
        &geometry, 5,
        &mut |_, p: [f64; 3]| Ok::<_, AcousticFeError>(Medium {
            density_kg_m3: 1.1 + 0.2 * p[0],
            sound_speed_m_s: 3.0 + 0.4 * p[1],
        }), &mut || false,
    ).unwrap();
    let mass = dense(&[2.0, 0.4, 0.4, 1.1], 2);
    let stiffness = dense(&[70.0, -13.0, -13.0, 40.0], 2);
    let damping = dense(&[0.7, 0.3, -0.1, 0.5], 2);
    let traces = [vec![([0, 1, 2], [0.7, -0.3, 0.9])],
        vec![([2, 0, 1], [1.1, -0.2, 0.8])]];
    let wall = ResistiveBoundary { faces: vec![[0, 1, 3]], admittance_m_pa_s: 0.08 };
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior, StructuralMatrices { mass: &mass, stiffness: &stiffness, damping: Some(&damping) },
            &traces, &[wall], controls(), cx,
        ).unwrap();
        let full = model.prepare_exact_time(0.04, StructuralInteriorTimeLimits::default(), cx).unwrap();
        let half = model.prepare_exact_time(0.02, StructuralInteriorTimeLimits::default(), cx).unwrap();
        let potential: Vec<_> = (0..geometry.ndof()).map(|i| 0.005 * (i + 1) as f64).collect();
        let pressure: Vec<_> = (0..geometry.ndof()).map(|i| 0.02 * (i % 3) as f64 - 0.015).collect();
        let initial = full.state(0.0, &[0.01, -0.02], &[0.04, -0.03], &potential, &pressure, cx).unwrap();
        let complex_flow = model.project_inward_velocity(&[([1, 2, 3],
            [C64::new(0.1, 0.0), C64::new(-0.04, 0.0), C64::new(0.02, 0.0)])]).unwrap();
        let flow: Vec<_> = complex_flow.iter().map(|v| v.re).collect();
        let force = [0.7, -0.3];
        let once = full.step(&initial, &force, &flow, cx).unwrap();
        let first = half.step(&initial, &force, &flow, cx).unwrap();
        let second = half.step(&first.state, &force, &flow, cx).unwrap();
        assert_eq!(once.state.pressure_pa().len(), geometry.ndof());
        assert!(geometry.ndof() > points.len());
        for (a, b) in once.state.displacement.iter().zip(&second.state.displacement) { close(*a, *b, 3e-10); }
        for (a, b) in once.state.velocity.iter().zip(&second.state.velocity) { close(*a, *b, 3e-10); }
        for (whole, a, b) in [
            (once.frame.work.structural_input_work_j, first.frame.work.structural_input_work_j, second.frame.work.structural_input_work_j),
            (once.frame.work.acoustic_input_work_j, first.frame.work.acoustic_input_work_j, second.frame.work.acoustic_input_work_j),
            (once.frame.work.structural_dissipation_j, first.frame.work.structural_dissipation_j, second.frame.work.structural_dissipation_j),
            (once.frame.work.wall_dissipation_j, first.frame.work.wall_dissipation_j, second.frame.work.wall_dissipation_j),
            (once.frame.work.fluid_to_structure_work_j, first.frame.work.fluid_to_structure_work_j, second.frame.work.fluid_to_structure_work_j),
        ] { close(whole, a + b, 3e-11); }
        for (whole, (a, b)) in once.frame.integrated_pressure_pa_s.iter()
            .zip(first.frame.integrated_pressure_pa_s.iter().zip(&second.frame.integrated_pressure_pa_s)) {
            close(*whole, *a + *b, 3e-10);
        }
        assert!(once.frame.work.structural_dissipation_j > 0.0);
        assert!(once.frame.work.wall_dissipation_j > 0.0);
        assert!(once.frame.work.fluid_to_structure_work_j.abs() > 1e-10);
        assert!(once.frame.work.structural_balance_residual_j().abs() < 3e-11);
        assert!(once.frame.work.acoustic_balance_residual_j().abs() < 3e-11);
        assert!(once.frame.work.balance_residual_j().abs() < 3e-11);
        assert_eq!(once.state.time_s(), 0.04);
        assert_eq!(second.state.time_s(), 0.04);
    });
}
