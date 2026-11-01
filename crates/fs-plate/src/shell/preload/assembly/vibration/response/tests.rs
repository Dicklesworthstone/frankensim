use super::super::super::equilibrate_shell_assembly;
use super::super::super::tests::{Part, budget as tie_budget, options as preload_options, seam};
use super::super::tests::{dense_reference, modal_options, output_budget, plan, solve, with_cx};
use super::super::{WindowPlan, eigen_windows};
use super::*;
use crate::shell::junction::{GlueLine, TieKind, TieOffset};
use crate::shell::preload::tests::reference_newton;
use fs_exec::CancelGate;

fn budget() -> AssemblyResponseBudget {
    AssemblyResponseBudget {
        maximum_projection_terms: 100_000_000,
        maximum_projection_bytes: 2 * 1024 * 1024,
        linear: LinearSplitBudget {
            max_modes: 16,
            max_omega_dt: 1.,
        },
    }
}

fn loaded_part() -> Part {
    let mut part = Part::rectangle(0, 2, 2, 0., true);
    for (node, force) in [(2, 25.), (5, 50.), (8, 25.)] {
        part.loads[6 * node] = force;
    }
    part.loads[6 * 8 + 2] = -0.01;
    part
}

fn relative(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * scale.max(f64::MIN_POSITIVE),
        "{a} != {b}; scale {scale}"
    );
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

#[test]
fn loaded_equilibrium_remains_exactly_at_rest_and_capacity_is_explicit() {
    let part = loaded_part();
    let body = solve(&[part.input()]);
    assert!(body.displacements[0].iter().any(|v| v.abs() > 1e-7));
    let windows = plan(&body, 4);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration =
            eigen_windows(&body, &windows, &modal_options(), output_budget(), cx).unwrap();
        let mut exact_cap = budget();
        exact_cap.linear.max_modes = 8;
        let mut response = AssemblyResponse::prepare(&vibration, 1e-4, exact_cap, cx).unwrap();
        let initial = bits(response.canonical_state());
        let zero = vec![0.; part.loads.len()];
        for _ in 0..8 {
            let frame = response.step(&[&zero], cx).unwrap();
            assert_eq!(frame.displacements, body.displacements);
            assert!(
                frame
                    .velocities
                    .iter()
                    .chain(&frame.mean_velocities)
                    .flatten()
                    .all(|v| *v == 0.)
            );
            assert_eq!(frame.perturbation_energy_j, 0.);
            assert_eq!(frame.energy_change_j, 0.);
            assert_eq!(frame.supplied_work_j, 0.);
            assert_eq!(frame.glue_dissipation_j, 0.);
            assert_eq!(bits(response.canonical_state()), initial);
        }
        let mut below = exact_cap;
        below.linear.max_modes -= 1;
        assert!(AssemblyResponse::prepare(&vibration, 1e-4, below, cx).is_err());
        below = exact_cap;
        below.maximum_projection_terms = 0;
        assert!(AssemblyResponse::prepare(&vibration, 1e-4, below, cx).is_err());
        below = exact_cap;
        below.maximum_projection_bytes = 1;
        assert!(AssemblyResponse::prepare(&vibration, 1e-4, below, cx).is_err());
        for duration in [0., -1., f64::NAN] {
            assert!(AssemblyResponse::prepare(&vibration, duration, exact_cap, cx).is_err());
        }
        assert!(response.restore_state(&[f64::NAN; 8], cx).is_err());
        assert!(response.initialize(&[0.; 3], &[0.; 4], cx).is_err());
        assert!(response.step(&[&[]], cx).is_err());
        assert_eq!(bits(response.canonical_state()), initial);
    });
}

#[test]
fn independent_forced_oscillator_matches_at_large_angular_steps() {
    let part = loaded_part();
    let body = solve(&[part.input()]);
    let reference = dense_reference(&body);
    let windows = WindowPlan {
        windows: vec![(0., 0.5 * (reference[0].lambda + reference[1].lambda))],
        threads: 1,
        maximum_windows: 1,
        maximum_bytes: 64 * 1024 * 1024,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    };
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration =
            eigen_windows(&body, &windows, &modal_options(), output_budget(), cx).unwrap();
        assert_eq!(vibration.mode_count(), 1);
        let probe = AssemblyResponse::prepare(&vibration, 1e-4, budget(), cx).unwrap();
        let (mass, stiffness) = (probe.mass_matrix()[0], probe.stiffness_matrix()[0]);
        let omega = (stiffness / mass).sqrt();
        let mut force = vec![0.; part.loads.len()];
        force[6 * 8 + 2] = 0.02;
        let modal_force = vibration.project_loads(&[&force], cx).unwrap()[0];
        let (q0, v0) = (2e-5, 3e-3);
        let equilibrium = modal_force / stiffness;
        let (sin, cos) = 2.5_f64.sin_cos(); // independent analytic oscillator oracle
        let expected_q = equilibrium + (q0 - equilibrium) * cos + v0 / omega * sin;
        let expected_v = -(q0 - equilibrium) * omega * sin + v0 * cos;
        let qscale = q0.abs().max(equilibrium.abs()).max((v0 / omega).abs());
        let vscale = omega * qscale + v0.abs();
        for steps in [1, 2, 4] {
            let duration = 2.5 / (omega * steps as f64);
            let mut response =
                AssemblyResponse::prepare(&vibration, duration, budget(), cx).unwrap();
            response.initialize(&[q0], &[v0], cx).unwrap();
            let mut supplied = 0.;
            for _ in 0..steps {
                let frame = response.step(&[&force], cx).unwrap();
                supplied += frame.supplied_work_j;
                assert_eq!(frame.glue_dissipation_j, 0.);
            }
            relative(response.canonical_state()[0], expected_q, qscale, 2e-11);
            relative(
                response.canonical_state()[1] / mass,
                expected_v,
                vscale,
                2e-11,
            );
            relative(
                supplied,
                modal_force * (expected_q - q0),
                modal_force.abs() * qscale,
                2e-11,
            );
        }
    });
}

// Independent physical-equation RK4 reference: M*qdd = F-K*q-C*qd.
// The extra two coordinates integrate glue heat and held-force work directly.
fn rk4(
    mass: &LdltFactor,
    k: &[f64],
    c: &[f64],
    force: &[f64],
    initial: &[f64],
    duration: f64,
    steps: usize,
) -> Vec<f64> {
    let n = force.len();
    let derivative = |state: &[f64]| {
        let (q, v) = (&state[..n], &state[n..2 * n]);
        let mut rhs = force.to_vec();
        let mut heat = 0.;
        for i in 0..n {
            for j in 0..n {
                rhs[i] -= k[i * n + j] * q[j] + c[i * n + j] * v[j];
                heat += v[i] * c[i * n + j] * v[j];
            }
        }
        let mut result = v.to_vec();
        result.extend(mass.solve(&rhs));
        result.push(heat);
        result.push(force.iter().zip(v).map(|(a, b)| a * b).sum());
        result
    };
    let mut state = initial.to_vec();
    state.resize(2 * n + 2, 0.);
    let dt = duration / steps as f64;
    for _ in 0..steps {
        let k1 = derivative(&state);
        let shifted = |slope: &[f64], scale: f64| {
            state
                .iter()
                .zip(slope)
                .map(|(a, b)| a + scale * dt * b)
                .collect::<Vec<_>>()
        };
        let k2 = derivative(&shifted(&k1, 0.5));
        let k3 = derivative(&shifted(&k2, 0.5));
        let k4 = derivative(&shifted(&k3, 1.));
        for i in 0..state.len() {
            state[i] += dt * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]) / 6.;
        }
    }
    state
}

fn state_error(a: &[f64], b: &[f64], n: usize) -> f64 {
    let mut worst = 0.0_f64;
    for part in [0..n, n..2 * n] {
        let scale = a[part.clone()]
            .iter()
            .chain(&b[part.clone()])
            .fold(0.0_f64, |s, v| s.max(v.abs()))
            .max(f64::MIN_POSITIVE);
        for i in part {
            worst = worst.max((a[i] - b[i]).abs() / scale);
        }
    }
    worst
}

#[test]
fn full_pencil_and_coupled_glue_match_original_equations_and_heat() {
    let mut left = Part::rectangle(0, 2, 2, 0., true);
    let right = Part::rectangle(2, 2, 2, 0., true);
    left.loads[6 * 5 + 2] = -0.01;
    let mut tie = seam(2, TieOffset::Midsurface);
    tie.kind = TieKind::Compliant(GlueLine {
        normal_n_m3: 1e5,
        shear_n_m3: [1e5, 2e5],
        viscosity_n_s_m3: Some([1e3, 2e3, 3e3]),
        width_m: 0.01,
    });
    let body = equilibrate_shell_assembly(
        &[left.input(), right.input()],
        &[tie],
        &[],
        tie_budget(),
        preload_options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let windows = plan(&body, 6);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration =
            eigen_windows(&body, &windows, &modal_options(), output_budget(), cx).unwrap();
        let n = vibration.mode_count();
        let duration = 0.002;
        let modes: Vec<_> = vibration.modes().map(|mode| mode.phi.as_slice()).collect();
        let mut polls = 0;
        let source_rows = body
            .assembly
            .project_glue_damping_source_rows(&modes, 2 * 1024 * 1024, 100_000_000, || {
                polls += 1;
                false
            })
            .unwrap();
        assert_eq!(
            source_rows.len(),
            body.assembly.glue_damping_source_row_count() * n
        );
        assert!(body.assembly.glue_damping_source_row_count() > 2 * n);
        for (bytes, terms) in [(0, u64::MAX), (usize::MAX, 0)] {
            assert!(
                body.assembly
                    .project_glue_damping_source_rows(&modes, bytes, terms, || false)
                    .is_err()
            );
        }
        let mut cancelled_polls = 0;
        assert!(matches!(
            body.assembly.project_glue_damping_source_rows(
                &modes,
                2 * 1024 * 1024,
                100_000_000,
                || {
                    cancelled_polls += 1;
                    cancelled_polls == polls
                },
            ),
            Err(crate::PlateError::Cancelled)
        ));
        let mut response = AssemblyResponse::prepare(&vibration, duration, budget(), cx).unwrap();
        assert_eq!(
            bits(response.damping_matrix()),
            bits(vibration.glue_damping())
        );
        assert!((0..n).any(|i| (0..i).any(|j| response.damping_matrix()[i * n + j].abs() > 1e-8)));
        // Finite eigen residuals and orthogonality error are retained, not
        // silently overwritten by a diagonal lambda/identity pencil.
        assert!((0..n).any(|i| (0..i).any(|j| response.stiffness_matrix()[i * n + j] != 0.)));
        assert!((0..n).any(|i| {
            (0..n).any(|j| response.mass_matrix()[i * n + j] != if i == j { 1. } else { 0. })
        }));
        let mass = positive_factor(
            response.mass_matrix(),
            n,
            budget().maximum_projection_bytes,
            cx,
        )
        .unwrap();
        let q: Vec<_> = (0..n).map(|i| (i + 1) as f64 * 1e-5).collect();
        let v: Vec<_> = (0..n)
            .map(|i| (i + 1) as f64 * if i % 2 == 0 { 0.005 } else { -0.005 })
            .collect();
        let source_power: f64 = source_rows
            .chunks_exact(n)
            .map(|row| row.iter().zip(&v).map(|(a, b)| a * b).sum::<f64>().powi(2))
            .sum();
        let retained_velocity = vibration.retained_perturbation(&v, cx).unwrap();
        let original_power = body
            .assembly
            .glue_dissipation_rate(&retained_velocity)
            .unwrap();
        relative(source_power, original_power, original_power, 1e-12);
        response.initialize(&q, &v, cx).unwrap();
        let mut initial = q.clone();
        initial.extend(&v);
        let mut forces = [vec![0.; left.loads.len()], vec![0.; right.loads.len()]];
        forces[0][6 * 8 + 2] = 0.02;
        forces[1][6 * 8 + 2] = -0.01;
        let loads = [&forces[0][..], &forces[1][..]];
        let force = vibration.project_loads(&loads, cx).unwrap();
        let frame = response.step(&loads, cx).unwrap();
        let mut actual = response.canonical_state()[..n].to_vec();
        actual.extend(mass.solve(&response.canonical_state()[n..]));
        let references: Vec<_> = [256, 512, 1024]
            .into_iter()
            .map(|steps| {
                rk4(
                    &mass,
                    response.stiffness_matrix(),
                    response.damping_matrix(),
                    &force,
                    &initial,
                    duration,
                    steps,
                )
            })
            .collect();
        let errors: Vec<_> = references
            .iter()
            .map(|r| state_error(&actual, r, n))
            .collect();
        assert!(
            errors[1] < 0.2 * errors[0] && errors[2] < 0.2 * errors[1],
            "RK4 ladder {errors:?}"
        );
        assert!(errors[2] < 2e-8, "{errors:?}");
        let reference = &references[2];
        relative(
            frame.glue_dissipation_j,
            reference[2 * n],
            frame.perturbation_energy_j,
            2e-8,
        );
        relative(
            frame.supplied_work_j,
            reference[2 * n + 1],
            frame.perturbation_energy_j,
            2e-8,
        );
        assert!(frame.glue_dissipation_j > 0.);
        assert!(frame.energy_balance_residual_j.abs() < 1e-12 * frame.perturbation_energy_j);
        let mut diagonal = response.damping_matrix().to_vec();
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    diagonal[i * n + j] = 0.;
                }
            }
        }
        let falsifier = rk4(
            &mass,
            response.stiffness_matrix(),
            &diagonal,
            &force,
            &initial,
            duration,
            1024,
        );
        assert!(
            state_error(&actual, &falsifier, n) > 1e-6,
            "cross-mode damping must affect physical motion"
        );
        // Independently reconstruct energy in the original physical FE pencil.
        let mut physical_q = vec![0.; body.assembly.model.free];
        let mut physical_v = physical_q.clone();
        for (i, mode) in vibration.modes().enumerate() {
            for j in 0..physical_q.len() {
                physical_q[j] += mode.phi[j] * actual[i];
                physical_v[j] += mode.phi[j] * actual[n + i];
            }
        }
        let mut kq = physical_q.clone();
        let mut mv = physical_v.clone();
        body.assembly.model.k.spmv(&physical_q, &mut kq);
        body.assembly.model.m.spmv(&physical_v, &mut mv);
        let original = 0.5
            * (physical_q.iter().zip(kq).map(|(a, b)| a * b).sum::<f64>()
                + physical_v.iter().zip(mv).map(|(a, b)| a * b).sum::<f64>());
        relative(frame.perturbation_energy_j, original, original, 2e-11);
    });
}

#[test]
fn nonmatching_slave_forces_exchange_exact_physical_interval_work() {
    let left = Part::rectangle(0, 2, 2, 0., true);
    let mut right = Part::rectangle(2, 2, 3, 0., false);
    right.loads[6 * 8 + 2] = -0.1;
    let body = equilibrate_shell_assembly(
        &[left.input(), right.input()],
        &[seam(3, TieOffset::Midsurface)],
        &[],
        tie_budget(),
        preload_options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let windows = plan(&body, 4);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let vibration =
            eigen_windows(&body, &windows, &modal_options(), output_budget(), cx).unwrap();
        let n = vibration.mode_count();
        let mut response = AssemblyResponse::prepare(&vibration, 1e-4, budget(), cx).unwrap();
        response
            .initialize(&vec![1e-5; n], &vec![0.01; n], cx)
            .unwrap();
        let before = vibration
            .displacements(&response.canonical_state()[..n], cx)
            .unwrap();
        let mut forces: Vec<Vec<f64>> = [&left, &right]
            .iter()
            .map(|part| {
                (0..part.loads.len())
                    .map(|i| ((i * 7 % 11) as f64 - 5.) * 0.001)
                    .collect()
            })
            .collect();
        forces[1][6 * 3 + 2] += 0.1; // load on an eliminated nonmatching slave
        let frame = response
            .step(&forces.iter().map(Vec::as_slice).collect::<Vec<_>>(), cx)
            .unwrap();
        let physical_work: f64 = forces
            .iter()
            .flatten()
            .zip(frame.displacements.iter().flatten())
            .zip(before.iter().flatten())
            .map(|((f, b), a)| f * (b - a))
            .sum();
        let port_work: f64 = forces
            .iter()
            .flatten()
            .zip(frame.mean_velocities.iter().flatten())
            .map(|(f, v)| f * v * response.duration_s())
            .sum();
        relative(
            frame.supplied_work_j,
            physical_work,
            physical_work.abs(),
            2e-11,
        );
        relative(frame.supplied_work_j, port_work, port_work.abs(), 2e-12);
        for ((a, b), mean) in before
            .iter()
            .flatten()
            .zip(frame.displacements.iter().flatten())
            .zip(frame.mean_velocities.iter().flatten())
        {
            relative(
                b - a,
                mean * response.duration_s(),
                (b - a).abs().max(1e-10),
                2e-11,
            );
        }
    });
}

#[test]
fn final_cancellation_is_atomic_and_checkpoint_retry_is_bitwise() {
    let part = loaded_part();
    let body = solve(&[part.input()]);
    let windows = plan(&body, 4);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let vibration =
            eigen_windows(&body, &windows, &modal_options(), output_budget(), cx).unwrap();
        let mut response = AssemblyResponse::prepare(&vibration, 1e-4, budget(), cx).unwrap();
        response.initialize(&[1e-5; 4], &[0.01; 4], cx).unwrap();
        let initial = response.canonical_state().to_vec();
        let mut force = vec![0.; part.loads.len()];
        force[6 * 8 + 2] = 0.02;
        let mut polls = 0;
        let accepted = response
            .step_with_checkpoint(&[&force], cx, &mut || {
                polls += 1;
                checkpoint(cx)
            })
            .unwrap();
        let expected = bits(response.canonical_state());
        response.restore_state(&initial, cx).unwrap();
        let mut calls = 0;
        let error = response
            .step_with_checkpoint(&[&force], cx, &mut || {
                calls += 1;
                if calls == polls {
                    gate.request();
                }
                checkpoint(cx)
            })
            .unwrap_err();
        assert!(error.contains("cancelled"));
        assert_eq!(calls, polls);
        assert_eq!(bits(response.canonical_state()), bits(&initial));
        assert!(AssemblyResponse::prepare(&vibration, 1e-4, budget(), cx).is_err());
        with_cx(&CancelGate::new_clock_free(), |fresh| {
            let retry = response.step(&[&force], fresh).unwrap();
            assert_eq!(bits(response.canonical_state()), expected);
            assert_eq!(retry.displacements, accepted.displacements);
            assert_eq!(
                retry.supplied_work_j.to_bits(),
                accepted.supplied_work_j.to_bits()
            );
            let checkpoint = response.canonical_state().to_vec();
            response.step(&[&force], fresh).unwrap();
            let resumed = bits(response.canonical_state());
            response.restore_state(&checkpoint, fresh).unwrap();
            response.step(&[&force], fresh).unwrap();
            assert_eq!(bits(response.canonical_state()), resumed);
        });
    });
}
