//! Native time histories against closed-form Maxwell rates and independently
//! integrated original member, viscosity and mutual boundary equations.

use super::*;
use fs_phs::PortHamiltonian;
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_solid::frame2d::coupled::transient::FrameTransientControls;
use fs_solid::frame2d::relaxation::transient::FrameRelaxationInitialHistory;

const BROADBAND: (f64, f64) = (0.0, f64::MAX);

fn transient_controls() -> FrameTransientControls {
    FrameTransientControls {
        exact: DenseExponentialBudget {
            max_dimension: 128,
            max_dense_bytes: 64 << 20,
            max_multiply_terms: 1 << 42,
        },
        maximum_relative_balance_error: 1e-10,
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * expected.abs().max(1e-7),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

#[test]
fn all_three_maxwell_rates_and_forced_work_use_one_exact_native_history() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, selections) = frame(&[&[(0.5, 2.0 / 3.0)]], BROADBAND, true, (1.0, 0.0));
        let prepared = frame
            .prepare_relaxing_resonance(
                Rayleigh { a0: 4.5, a1: 0.0 },
                &selections,
                None,
                controls(),
                cx,
            )
            .unwrap();
        // M=2A, K=8A, C=9A, delta_K=.5A, 1/tau=1.5. The complete
        // characteristic polynomial is (s+1)(s+2)(s+3), not an isolated arm
        // decay followed by a held-velocity mechanical correction.
        let mass = 2.0 * area();
        let stiffness = 8.0 * area();
        let loss = 9.0 * area();
        let root = (0.5 * area()).sqrt();
        for duration in [0.03, 1.1, 5.0] {
            for rate in [-1.0_f64, -2.0, -3.0] {
                for equilibrium in [0.0, 0.013] {
                    let amplitude = 0.02;
                    let initial_z = root * rate / (rate + 1.5) * amplitude;
                    let initial = [
                        equilibrium + amplitude,
                        mass * rate * amplitude,
                        initial_z,
                        0.0,
                        0.0,
                    ];
                    let mut runtime = prepared
                        .prepare_transient(duration, transient_controls(), cx)
                        .unwrap();
                    runtime.restore_state(0.4, &initial, cx).unwrap();
                    let initial_energy = runtime.report().motion.total_energy_j;
                    let mut force = [0.0; 6];
                    force[3] = stiffness * equilibrium;
                    // Loads applied to stationary supports must do zero work.
                    force[0] = 13.0;
                    runtime.step(&force, cx).unwrap();
                    let decay = (rate * duration).exp();
                    let square_integral = -(2.0 * rate * duration).exp_m1() / (-2.0 * rate);
                    let report = runtime.report();
                    near(report.motion.time_s, 0.4 + duration);
                    near(
                        runtime.canonical_state()[0],
                        equilibrium + amplitude * decay,
                    );
                    near(
                        runtime.canonical_state()[1],
                        mass * rate * amplitude * decay,
                    );
                    near(runtime.canonical_state()[2], initial_z * decay);
                    near(report.motion.velocity[3], rate * amplitude * decay);
                    near(
                        report.material_dissipation_j,
                        1.5 * initial_z * initial_z * square_integral,
                    );
                    near(
                        report.motion.frame_dissipation_j,
                        loss * (rate * amplitude).powi(2) * square_integral,
                    );
                    near(
                        report.motion.applied_work_j,
                        force[3] * amplitude * (decay - 1.0),
                    );
                    near(
                        report.motion.total_energy_j - initial_energy,
                        report.motion.applied_work_j
                            - report.material_dissipation_j
                            - report.motion.frame_dissipation_j,
                    );
                    near(
                        report.material_memory_energy_j,
                        0.5 * (initial_z * decay).powi(2),
                    );
                    near(
                        report.motion.attached_dissipation_j,
                        report.material_dissipation_j,
                    );
                    assert_eq!(report.boundary_dissipation_j, 0.0);
                    assert_eq!(report.motion.stationary_support_work_j, 0.0);
                    assert!(report.motion.maximum_relative_equation_residual < 1e-10);
                }
            }
        }
    });
}

#[test]
fn constrained_dark_memory_relaxes_and_initial_history_remains_explicit() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, selections) = frame(
            &[&[(0.0, 0.1), (0.5, 2.0 / 3.0)]],
            BROADBAND,
            true,
            (1.0, 0.0),
        );
        let prepared = frame
            .prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx)
            .unwrap();
        assert_eq!(prepared.members()[0].source.spectra()[0].terms.len(), 2);
        assert_eq!(prepared.coordinates()[0].source_arm, 1);
        let mut runtime = prepared
            .prepare_transient(0.7, transient_controls(), cx)
            .unwrap();
        let mut displacement = [0.0; 6];
        displacement[3] = 0.02;
        runtime
            .initialize(
                0.0,
                &displacement,
                &[0.0; 6],
                FrameRelaxationInitialHistory::Relaxed,
                &[],
                cx,
            )
            .unwrap();
        assert_eq!(runtime.report().material_memory_energy_j, 0.0);
        runtime
            .initialize(
                0.0,
                &displacement,
                &[0.0; 6],
                FrameRelaxationInitialHistory::Unrelaxed,
                &[],
                cx,
            )
            .unwrap();
        near(runtime.report().memory[0].elastic_extension_m, 0.02);
        near(runtime.report().memory[0].viscous_extension_m, 0.0);
        near(
            runtime.report().material_memory_energy_j,
            0.5 * (0.5 * area()) * 0.02_f64.powi(2),
        );

        // A constrained bending strain can carry old internal stress even with
        // zero current nodal strain. It dissipates, and is never deleted as a
        // mode invisible to the free mechanical coordinate.
        let dark = &prepared.coordinates()[1];
        assert_eq!(dark.channel, FrameStrainChannel::SymmetricBending);
        let elastic = [0.0, 0.2 / dark.stiffness_n_m.sqrt(), 0.0];
        runtime
            .initialize(
                0.0,
                &[0.0; 6],
                &[0.0; 6],
                FrameRelaxationInitialHistory::ElasticExtensions(&elastic),
                &[],
                cx,
            )
            .unwrap();
        runtime.step(&[0.0; 6], cx).unwrap();
        let decay = (-1.5_f64 * 0.7).exp();
        near(runtime.canonical_state()[dark.native_state], 0.2 * decay);
        near(
            runtime.report().material_dissipation_j,
            0.02 * (1.0 - decay * decay),
        );
        near(
            runtime.report().material_memory_energy_j,
            0.02 * decay * decay,
        );
        assert_eq!(runtime.report().motion.displacement, vec![0.0; 6]);
        assert_eq!(runtime.report().motion.frame_dissipation_j, 0.0);
        assert_eq!(runtime.report().motion.stationary_support_work_j, 0.0);
    });
}

// Independent original two-member axial equations. z is the physical Maxwell
// energy coordinate. x_b are the unchanged two mutually coupled boundary
// states; the last three coordinates integrate the three physical heat rates.
fn member_derivative(x: &[f64; 10], force: [f64; 2]) -> [f64; 10] {
    let a = area();
    let root = (0.5 * a).sqrt();
    let effort = [1.4 * x[5] + 0.2 * x[6], 0.2 * x[5] + 0.9 * x[6]];
    let drag = 2.0 * a * (x[2] - x[3]);
    let momentum = [
        force[0] - 8.0 * a * (2.0 * x[0] - x[1]) - drag - root * x[4],
        force[1] - 8.0 * a * (x[1] - x[0]) + drag - (0.8 * effort[0] - 0.35 * effort[1]),
    ];
    [
        x[2],
        x[3],
        (2.0 * momentum[0] - momentum[1]) / (7.0 * a),
        (-momentum[0] + 4.0 * momentum[1]) / (7.0 * a),
        root * x[2] - 1.5 * x[4],
        -0.3 * effort[0] + 0.6 * effort[1] + 0.8 * x[3],
        -0.8 * effort[0] - 0.4 * effort[1] - 0.35 * x[3],
        2.0 * a * (x[2] - x[3]).powi(2),
        1.5 * x[4] * x[4],
        0.3 * effort[0] * effort[0] + 0.2 * effort[0] * effort[1] + 0.4 * effort[1] * effort[1],
    ]
}

fn independent_member_history(mut x: [f64; 10], force: [f64; 2], duration: f64) -> [f64; 10] {
    // Test-only fourth-order reference: no production generator, split, modal
    // eigenvector, storage gradient or heat Gramian is reused by this oracle.
    let dt = duration / 2048.0;
    for _ in 0..2048 {
        let a = member_derivative(&x, force);
        let b = member_derivative(&core::array::from_fn(|i| x[i] + 0.5 * dt * a[i]), force);
        let c = member_derivative(&core::array::from_fn(|i| x[i] + 0.5 * dt * b[i]), force);
        let d = member_derivative(&core::array::from_fn(|i| x[i] + dt * c[i]), force);
        for i in 0..x.len() {
            x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn full_consistent_mass_nonproportional_loss_and_mutual_boundary_keep_all_cross_terms() {
    with_cx(&CancelGate::new(), |cx| {
        let (mut frame, selections) = frame(&[&[(0.5, 2.0 / 3.0)]], BROADBAND, true, (1.0, 0.0));
        frame.add_node(2.0, 0.0);
        frame.fix(2, [false, true, true]);
        frame
            .add_material_viscous_element(1, 2, 5, super::super::source(2.0, 0.0, BROADBAND))
            .unwrap();
        let boundary = PortHamiltonian::new(
            2,
            1,
            vec![0.0, 0.7, -0.7, 0.0],
            vec![0.3, 0.1, 0.1, 0.4],
            vec![0.8, -0.35],
            Box::new(QuadraticStorage::new(vec![1.4, 0.2, 0.2, 0.9], 2).unwrap()),
        )
        .unwrap();
        let mut map = [0.0; 9];
        map[6] = 1.0;
        let prepared = frame
            .prepare_relaxing_resonance(
                NO_RAYLEIGH,
                &selections,
                Some(FrameDynamicPortLoad {
                    force_map: &map,
                    dynamics: FramePortDynamics::Ode(&boundary),
                }),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(prepared.coupled().free_dofs(), [3, 6]);
        assert_eq!(prepared.memory_state_range(), 4..7);
        assert_eq!(prepared.boundary_state_range(), 7..9);
        let initial = [0.01, -0.006, 0.02, -0.01, 0.01, 0.04, -0.03, 0.0, 0.0, 0.0];
        let forcing = [0.003, -0.002];
        let expected = independent_member_history(initial, forcing, 0.4);
        let mut displacement = [0.0; 9];
        displacement[3] = initial[0];
        displacement[6] = initial[1];
        let mut velocity = [0.0; 9];
        velocity[3] = initial[2];
        velocity[6] = initial[3];
        let elastic = [initial[4] / (0.5 * area()).sqrt(), 0.0, 0.0];
        let mut force = [0.0; 9];
        force[3] = forcing[0];
        force[6] = forcing[1];
        let mut endpoints = Vec::new();
        for intervals in [1, 4] {
            let mut runtime = prepared
                .prepare_transient(0.4 / intervals as f64, transient_controls(), cx)
                .unwrap();
            runtime
                .initialize(
                    0.0,
                    &displacement,
                    &velocity,
                    FrameRelaxationInitialHistory::ElasticExtensions(&elastic),
                    &initial[5..7],
                    cx,
                )
                .unwrap();
            let before = runtime.report().motion.total_energy_j;
            let mut heat = [0.0; 3];
            let mut work = 0.0;
            for _ in 0..intervals {
                let report = runtime.step(&force, cx).unwrap();
                heat[0] += report.motion.frame_dissipation_j;
                heat[1] += report.material_dissipation_j;
                heat[2] += report.boundary_dissipation_j;
                work += report.motion.applied_work_j;
                near(
                    report.motion.attached_dissipation_j,
                    report.material_dissipation_j + report.boundary_dissipation_j,
                );
            }
            let report = runtime.report();
            for (actual, expected) in [
                report.motion.displacement[3],
                report.motion.displacement[6],
                report.motion.velocity[3],
                report.motion.velocity[6],
                runtime.canonical_state()[4],
                runtime.canonical_state()[7],
                runtime.canonical_state()[8],
            ]
            .into_iter()
            .zip(&expected[..7])
            {
                near(actual, *expected);
            }
            for i in 0..3 {
                near(heat[i], expected[7 + i]);
            }
            near(
                work,
                forcing[0] * (expected[0] - initial[0]) + forcing[1] * (expected[1] - initial[1]),
            );
            near(
                report.motion.total_energy_j - before + heat.iter().sum::<f64>(),
                work,
            );
            endpoints.push(runtime.canonical_state().to_vec());
        }
        for (&a, &b) in endpoints[0].iter().zip(&endpoints[1]) {
            near(a, b);
        }
    });
}

#[test]
fn refusals_preserve_accepted_motion_history_and_retry() {
    with_cx(&CancelGate::new(), |cx| {
        let (frame, selections) = frame(&[&[(0.5, 2.0 / 3.0)]], BROADBAND, true, (1.0, 0.0));
        let prepared = frame
            .prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx)
            .unwrap();
        let mut runtime = prepared
            .prepare_transient(0.05, transient_controls(), cx)
            .unwrap();
        runtime
            .restore_state(0.3, &[0.01, -0.02, 0.03, 0.04, 0.05], cx)
            .unwrap();
        let accepted = runtime.canonical_state().to_vec();
        let accepted_motion = runtime.report().motion.clone();
        let accepted_memory: Vec<_> = runtime
            .report()
            .memory
            .iter()
            .map(|m| (m.elastic_extension_m, m.stored_energy_j))
            .collect();
        assert!(runtime.restore_state(-1.0, &[0.1; 5], cx).is_err());
        let mut force = [0.0; 6];
        force[3] = f64::NAN;
        assert!(runtime.step(&force, cx).is_err());
        let mut moved_support = [0.0; 6];
        moved_support[0] = 0.1;
        assert!(
            runtime
                .initialize(
                    1.0,
                    &moved_support,
                    &[0.0; 6],
                    FrameRelaxationInitialHistory::Unrelaxed,
                    &[],
                    cx
                )
                .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        with_cx(&gate, |cancelled| {
            assert!(matches!(
                runtime.step(&[0.0; 6], cancelled),
                Err(FrameRelaxationError::Native(FrameCoupledError::Cancelled))
            ));
        });
        assert_eq!(runtime.canonical_state(), accepted);
        assert_eq!(runtime.report().motion, &accepted_motion);
        assert_eq!(
            runtime
                .report()
                .memory
                .iter()
                .map(|m| (m.elastic_extension_m, m.stored_energy_j))
                .collect::<Vec<_>>(),
            accepted_memory
        );
        let mut retry = prepared
            .prepare_transient(0.05, transient_controls(), cx)
            .unwrap();
        retry.restore_state(0.3, &accepted, cx).unwrap();
        runtime.step(&[0.0; 6], cx).unwrap();
        retry.step(&[0.0; 6], cx).unwrap();
        assert_eq!(runtime.canonical_state(), retry.canonical_state());
        assert_eq!(runtime.report().motion, retry.report().motion);
    });
}

#[test]
fn finite_material_bands_and_algebraic_boundary_constraints_refuse_time_extrapolation() {
    with_cx(&CancelGate::new(), |cx| {
        let (limited, selections) = frame(&[&[(0.5, 2.0 / 3.0)]], (0.5, 3.0), true, (1.0, 0.0));
        let limited = limited
            .prepare_relaxing_resonance(NO_RAYLEIGH, &selections, None, controls(), cx)
            .unwrap();
        assert!(matches!(
            limited.prepare_transient(0.1, transient_controls(), cx),
            Err(FrameRelaxationError::Native(FrameCoupledError::Invalid(_)))
        ));
        let (frame, selections) = frame(&[&[(0.5, 2.0 / 3.0)]], BROADBAND, true, (1.0, 0.0));
        let boundary = DescriptorPortHamiltonian::new(
            2,
            1,
            1,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0; 4],
            vec![0.0, 1.0],
            Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
        )
        .unwrap();
        let map = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let constrained = frame
            .prepare_relaxing_resonance(
                NO_RAYLEIGH,
                &selections,
                Some(FrameDynamicPortLoad {
                    force_map: &map,
                    dynamics: FramePortDynamics::Descriptor(&boundary),
                }),
                controls(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            constrained.prepare_transient(0.1, transient_controls(), cx),
            Err(FrameRelaxationError::Native(FrameCoupledError::Invalid(_)))
        ));
        assert!(constrained.coupled().factor(1.0, cx).is_ok());
    });
}
