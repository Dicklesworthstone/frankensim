//! Original P2 motion, nonlinear Maxwell strain and actual source-time checks.

use super::*;
use fs_solid::hyper3::mixed::MixedUpError;
use fs_solid::hyper3::mixed::visco::dynamics::{
    MixedViscoAdaptiveControls, MixedViscoAdaptiveError, MixedViscoPairLoads,
};

fn accuracy() -> MixedViscoAdaptiveControls {
    MixedViscoAdaptiveControls {
        absolute_displacement_m: 1.0,
        absolute_velocity_m_s: 1.0,
        absolute_green_strain: 1.0,
        absolute_inverse_viscous_metric: 1.0,
        absolute_cauchy_stress_pa: 1e6,
        relative: 0.0,
        minimum_pair_span_s: 1e-7,
        maximum_pair_span_s: 0.05,
        maximum_attempts: 16,
    }
}

fn loads(prescribed_m: Vec<(usize, f64)>) -> MixedViscoPairLoads {
    MixedViscoPairLoads {
        prescribed_m,
        nodal_forces_n: vec![],
        body_force_n_m3: [0.0; 3],
        pressure_gauge: None,
    }
}

fn norm(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0.0, f64::hypot)
}

fn oscillator_step(q: f64, v: f64, dt: f64) -> (f64, f64) {
    let mass = 8.0 * DENSITY * VOLUME / 105.0;
    let stiffness = 16.0 * MU * VOLUME / 5.0;
    let next = (mass * q + mass * dt * v) / (mass + dt * dt * stiffness);
    (next, (next - q) / dt)
}

#[test]
fn native_quadratic_shear_accepts_only_fine_history_and_reduces_physical_time_error() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, false);
        let edge = TET10_EDGES
            .iter()
            .position(|e| e.contains(&2) && e.contains(&3))
            .unwrap();
        let node = elements[0][4 + edge];
        let dof = 3 * node;
        let amplitude = 0.02;
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d != dof)
            .map(|d| (d, 0.0))
            .collect();
        let placed: Vec<_> = (0..3 * nodes.len())
            .map(|d| (d, if d == dof { amplitude } else { 0.0 }))
            .collect();
        let mass = 8.0 * DENSITY * VOLUME / 105.0;
        let stiffness = 16.0 * MU * VOLUME / 5.0;
        let omega = (stiffness / mass).sqrt();
        let target = 0.04;
        let mut final_errors = Vec::new();
        let mut accepted_counts = Vec::new();
        for tolerance in [1e-4, 2.5e-5] {
            let initial = solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::OriginalMaterial,
                cx,
            )
            .advance(
                0.0,
                &placed,
                &[],
                [0.0; 3],
                None,
                settings(),
                controls().material,
                cx,
            )
            .unwrap()
            .solid;
            let mut motion = MixedViscoTetMotion::try_new(
                initial,
                &[DENSITY],
                &vec![[0.0; 3]; nodes.len()],
                MASS_BYTES,
                cx,
            )
            .unwrap();
            let mut control = accuracy();
            control.absolute_displacement_m = tolerance;
            let (mut q, mut v) = (amplitude, 0.0);
            let mut proposed = 0.02;
            let mut numerical_loss = 0.0;
            let mut pairs = 0;
            let mut rejected = 0;
            while motion.solid().time_s() < target {
                assert!(pairs < 64);
                let start = motion.solid().time_s();
                let pair = motion
                    .advance_adaptive_pair_toward(
                        target,
                        proposed,
                        |_| Ok(loads(fixed.clone())),
                        settings(),
                        controls(),
                        control,
                        cx,
                    )
                    .unwrap();
                assert!(pair.estimate.maximum_ratio() <= 1.0);
                assert_eq!(pair.native_solves, 3 * (pair.rejected_attempts + 1));
                near(
                    pair.accepted_pair_span_s,
                    pair.second.motion.solid().time_s() - start,
                    0.0,
                );
                let dt1 = pair.first.motion.solid().time_s() - start;
                let dt2 = pair.second.motion.solid().time_s() - pair.first.motion.solid().time_s();
                let first = oscillator_step(q, v, dt1);
                let second = oscillator_step(first.0, first.1, dt2);
                let coarse = oscillator_step(q, v, dt1 + dt2);
                near(
                    pair.first.motion.solid().displacement_m()[node][0],
                    first.0,
                    8e-10,
                );
                near(
                    pair.second.motion.solid().displacement_m()[node][0],
                    second.0,
                    8e-10,
                );
                near(pair.second.motion.velocity_mps()[node][0], second.1, 3e-8);
                near(
                    pair.estimate.displacement_ratio,
                    (second.0 - coarse.0).abs() / tolerance,
                    2e-5,
                );
                work_identity(&pair.first);
                work_identity(&pair.second);
                assert_eq!(pair.first.material_energy.dissipated_j, 0.0);
                assert_eq!(pair.second.material_energy.dissipated_j, 0.0);
                near(
                    pair.endpoint_work_defect_magnitude_j,
                    0.5 * stiffness * ((first.0 - q).powi(2) + (second.0 - first.0).powi(2)),
                    2e-10,
                );
                near(
                    pair.numerical_kinetic_loss_j,
                    0.5 * mass * ((first.1 - v).powi(2) + (second.1 - first.1).powi(2)),
                    2e-10,
                );
                numerical_loss +=
                    pair.endpoint_work_defect_magnitude_j + pair.numerical_kinetic_loss_j;
                rejected += pair.rejected_attempts;
                proposed = pair.suggested_pair_span_s;
                motion = pair.second.motion;
                (q, v) = second;
                pairs += 1;
            }
            assert!(rejected > 0);
            near(motion.solid().time_s(), target, 0.0);
            near(
                0.5 * stiffness * amplitude * amplitude - 0.5 * (stiffness * q * q + mass * v * v),
                numerical_loss,
                3e-9,
            );
            final_errors.push(
                (q / amplitude - (omega * target).cos())
                    .hypot(v / (amplitude * omega) + (omega * target).sin()),
            );
            accepted_counts.push(pairs);
        }
        assert!(accepted_counts[1] > accepted_counts[0]);
        assert!(final_errors[1] < 0.7 * final_errors[0], "{final_errors:?}");
    });
}

// Independent scalar stationarity for an isochoric uniaxial NH Maxwell branch:
// a-a_trial + dt/(3 tau)*(exp(2a)-exp(-a)) = 0. No native material update enters.
fn relaxed_log_strain(trial: f64, dt: f64, tau: f64) -> f64 {
    let mut a = trial;
    for _ in 0..30 {
        let residual = a - trial + dt / (3.0 * tau) * ((2.0 * a).exp() - (-a).exp());
        let slope = 1.0 + dt / (3.0 * tau) * (2.0 * (2.0 * a).exp() + (-a).exp());
        a -= residual / slope;
    }
    a
}

#[test]
fn held_finite_strain_refines_original_maxwell_metrics_and_physical_pressure_stress() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        let stretch = 1.15_f64;
        let dilation = 1.03_f64;
        let j = dilation.powi(3);
        let pressure = 10.0 * MU * (j - 1.0);
        let branch_mu = 3.0 * MU;
        let tau = 0.07;
        let initial_a = stretch.ln();
        let mut invariant = None;
        for (cosine, sine) in [(1.0, 0.0), (0.6, 0.8)] {
            let f = [
                cosine * dilation * stretch,
                -sine * dilation / stretch.sqrt(),
                0.0,
                sine * dilation * stretch,
                cosine * dilation / stretch.sqrt(),
                0.0,
                0.0,
                0.0,
                dilation / stretch.sqrt(),
            ];
            let placed: Vec<_> = nodes
                .iter()
                .enumerate()
                .flat_map(|(n, x)| {
                    (0..3).map(move |i| {
                        let value = (0..3).map(|k| f[3 * i + k] * x[k]).sum::<f64>() - x[i];
                        (3 * n + i, value)
                    })
                })
                .collect();
            let initial = solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::OriginalMaterial,
                cx,
            )
            .advance(
                0.0,
                &placed,
                &[],
                [0.0; 3],
                None,
                settings(),
                controls().material,
                cx,
            )
            .unwrap()
            .solid;
            let motion = MixedViscoTetMotion::try_new(
                initial,
                &[DENSITY],
                &vec![[0.0; 3]; nodes.len()],
                MASS_BYTES,
                cx,
            )
            .unwrap();
            let before = motion.solid().quadrature_states().to_vec();
            let mut control = accuracy();
            control.absolute_inverse_viscous_metric = 5e-5;
            control.absolute_cauchy_stress_pa = 1.0;
            control.relative = 1e-3;
            let pair = motion
                .advance_adaptive_pair_toward(
                    0.1,
                    0.05,
                    |_| Ok(loads(placed.clone())),
                    settings(),
                    controls(),
                    control,
                    cx,
                )
                .unwrap();
            assert!(pair.rejected_attempts > 0);
            assert!(pair.estimate.inverse_viscous_metric_ratio > 0.05);
            assert!(pair.estimate.maximum_ratio() <= 1.0);
            let dt1 = pair.first.motion.solid().time_s();
            let dt2 = pair.second.motion.solid().time_s() - dt1;
            let coarse_a = relaxed_log_strain(initial_a, dt1 + dt2, tau);
            let first_a = relaxed_log_strain(initial_a, dt1, tau);
            let fine_a = relaxed_log_strain(first_a, dt2, tau);
            let cv = |a: f64| {
                [
                    (2.0 * a).exp() / stretch.powi(2),
                    (-a).exp() * stretch,
                    (-a).exp() * stretch,
                ]
            };
            let coarse_cv = cv(coarse_a);
            let fine_cv = cv(fine_a);
            let metric_change = |a: [f64; 3]| norm(a.into_iter().map(|x| x - 1.0));
            let metric_error = norm(coarse_cv.iter().zip(fine_cv).map(|(a, b)| a - b));
            near(
                pair.estimate.inverse_viscous_metric_ratio,
                metric_error
                    / (control.absolute_inverse_viscous_metric
                        + control.relative * metric_change(coarse_cv).max(metric_change(fine_cv))),
                2e-8,
            );
            let cauchy = |a: f64| {
                let tau_eq = MU * (stretch.powi(2) - 1.0 / stretch);
                let tau_branch = branch_mu * ((2.0 * a).exp() - (-a).exp());
                let value = (tau_eq + tau_branch) / (3.0 * j);
                [pressure + 2.0 * value, pressure - value, pressure - value]
            };
            let cs = cauchy(coarse_a);
            let fs = cauchy(fine_a);
            let stress_error = norm(cs.iter().zip(fs).map(|(a, b)| a - b));
            near(
                pair.estimate.cauchy_stress_ratio,
                stress_error
                    / (control.absolute_cauchy_stress_pa
                        + control.relative * norm(cs.into_iter()).max(norm(fs.into_iter()))),
                2e-9,
            );
            for state in pair.second.motion.solid().quadrature_states() {
                for i in 0..3 {
                    near(state.inverse_viscous_cauchy_green[0][i], fine_cv[i], 2e-10);
                }
                for i in 3..6 {
                    near(state.inverse_viscous_cauchy_green[0][i], 0.0, 2e-10);
                }
            }
            for p in pair.second.motion.solid().pressure_pa() {
                near(*p, pressure, 2e-8);
            }
            let branch_energy =
                |a: f64| 0.5 * branch_mu * ((2.0 * a).exp() + 2.0 * (-a).exp() - 3.0);
            let loss = VOLUME * (branch_energy(initial_a) - branch_energy(fine_a));
            near(
                pair.first.material_energy.dissipated_j + pair.second.material_energy.dissipated_j,
                loss,
                2e-10,
            );
            near(
                motion.solid().stored_energy_j() - pair.second.motion.solid().stored_energy_j(),
                loss,
                2e-10,
            );
            near(pair.numerical_kinetic_loss_j, 0.0, 2e-12);
            work_identity(&pair.first);
            work_identity(&pair.second);
            assert_eq!(motion.solid().quadrature_states(), before.as_slice());
            assert_eq!(motion.solid().time_s(), 0.0);
            let observed = (
                pair.accepted_pair_span_s,
                pair.estimate.inverse_viscous_metric_ratio,
                pair.estimate.cauchy_stress_ratio,
            );
            if let Some(previous) = invariant {
                let (span, metric, stress): (f64, f64, f64) = previous;
                near(observed.0, span, 2e-10);
                near(observed.1, metric, 2e-8);
                near(observed.2, stress, 2e-9);
            }
            invariant = Some(observed);
        }
    });
}

#[test]
fn explicit_time_sources_and_velocity_error_control_follow_actual_retry_endpoints() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(1e300, true);
        let initial_v = [0.1, -0.2, 0.05];
        let acceleration = |t: f64| [0.3 + 0.7 * t, -0.2 + 0.4 * t, 0.1 - 0.5 * t];
        let motion = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::ExactIncompressible,
                cx,
            ),
            &[DENSITY],
            &vec![initial_v; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let mut sampled = Vec::new();
        let mut control = accuracy();
        control.absolute_velocity_m_s = 1e-6;
        let pair = motion
            .advance_adaptive_pair_toward(
                0.018,
                0.05,
                |time| {
                    sampled.push(time);
                    let mut value = loads(vec![]);
                    value.body_force_n_m3 = acceleration(time).map(|a| DENSITY * a);
                    Ok(value)
                },
                settings(),
                controls(),
                control,
                cx,
            )
            .unwrap();
        assert!(pair.rejected_attempts > 0);
        assert_eq!(sampled.len(), 2 * (pair.rejected_attempts + 1));
        for sample in sampled.chunks_exact(2) {
            near(sample[0], 0.5 * sample[1], 0.0);
            assert!(sample[1] <= 0.018);
        }
        for samples in sampled.chunks_exact(2).collect::<Vec<_>>().windows(2) {
            assert!(samples[1][1] < samples[0][1]);
        }
        let middle = pair.first.motion.solid().time_s();
        let end = pair.second.motion.solid().time_s();
        assert!(end < 0.018); // A call advances one accepted pair, not the entire target.
        assert_eq!(&sampled[sampled.len() - 2..], &[middle, end]);
        let a1 = acceleration(middle);
        let a2 = acceleration(end);
        let first_v: [f64; 3] = core::array::from_fn(|i| initial_v[i] + middle * a1[i]);
        let fine_v: [f64; 3] = core::array::from_fn(|i| first_v[i] + (end - middle) * a2[i]);
        let fine_u: [f64; 3] =
            core::array::from_fn(|i| middle * first_v[i] + (end - middle) * fine_v[i]);
        let coarse_v: [f64; 3] = core::array::from_fn(|i| initial_v[i] + end * a2[i]);
        for (u, v) in pair
            .second
            .motion
            .solid()
            .displacement_m()
            .iter()
            .zip(pair.second.motion.velocity_mps())
        {
            for i in 0..3 {
                near(u[i], fine_u[i], 3e-10);
                near(v[i], fine_v[i], 3e-9);
            }
        }
        near(
            pair.estimate.velocity_ratio,
            norm(coarse_v.iter().zip(fine_v).map(|(c, f)| c - f)) / control.absolute_velocity_m_s,
            3e-3,
        );
        work_identity(&pair.first);
        work_identity(&pair.second);
        assert_eq!(motion.solid().time_s(), 0.0);
    });
}

#[test]
fn pair_work_budget_sources_and_attempt_limits_cannot_publish_partial_motion() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(1e300, true);
        let motion = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::ExactIncompressible,
                cx,
            ),
            &[DENSITY],
            &vec![[0.0; 3]; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let before = (
            motion.solid().displacement_m().to_vec(),
            motion.solid().pressure_pa().to_vec(),
            motion.solid().quadrature_states().to_vec(),
            motion.velocity_mps().to_vec(),
        );
        // a(mid)=2, a(end)=1: full loss=.0002 J; fine losses=.0002 and .00005 J.
        // Every native interval passes .00022 J, but the accepted pair must refuse.
        let mut work = controls();
        work.maximum_numerical_kinetic_loss_j = 0.00022;
        let failed = motion.advance_adaptive_pair_toward(
            0.02,
            0.02,
            |time| {
                let mut value = loads(vec![]);
                value.body_force_n_m3[0] = DENSITY * (3.0 - 100.0 * time);
                Ok(value)
            },
            settings(),
            work,
            accuracy(),
            cx,
        );
        match failed {
            Err(MixedViscoAdaptiveError::Native(MixedViscoTetError::EnergyGate {
                quantity,
                value_j,
                ..
            })) => {
                assert_eq!(quantity, "adaptive pair numerical kinetic loss");
                near(value_j, 0.00025, 2e-10);
            }
            other => panic!("expected pair-wide work refusal, got {other:?}"),
        }
        let source_error = MixedViscoTetError::Native(MixedUpError::InvalidInput {
            what: "declared source unavailable".into(),
        });
        let refused = motion
            .advance_adaptive_pair_toward(
                0.02,
                0.02,
                |_| Err(source_error.clone()),
                settings(),
                controls(),
                accuracy(),
                cx,
            )
            .unwrap_err();
        assert_eq!(refused, MixedViscoAdaptiveError::Native(source_error));
        let mut control = accuracy();
        control.maximum_attempts = 1;
        control.absolute_velocity_m_s = 1e-15;
        let rejected = motion
            .advance_adaptive_pair_toward(
                0.02,
                0.02,
                |time| {
                    let mut value = loads(vec![]);
                    value.body_force_n_m3[0] = DENSITY * time;
                    Ok(value)
                },
                settings(),
                controls(),
                control,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            rejected,
            MixedViscoAdaptiveError::Rejected { attempts: 1, .. }
        ));
        for attempts in [0, 33] {
            control.maximum_attempts = attempts;
            assert!(matches!(
                motion.advance_adaptive_pair_toward(
                    0.02,
                    0.02,
                    |_| panic!("invalid controls must precede source evaluation"),
                    settings(),
                    controls(),
                    control,
                    cx,
                ),
                Err(MixedViscoAdaptiveError::Invalid(_))
            ));
        }
        assert_eq!(
            (
                motion.solid().displacement_m().to_vec(),
                motion.solid().pressure_pa().to_vec(),
                motion.solid().quadrature_states().to_vec(),
                motion.velocity_mps().to_vec()
            ),
            before,
        );
        assert_eq!(motion.solid().time_s(), 0.0);
    });
}

#[test]
fn cancellation_after_source_sampling_leaves_all_accepted_histories_and_time_unchanged() {
    let cancel = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &cancel,
            arena,
            StreamKey {
                seed: 61,
                kernel_id: 22,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let (nodes, elements) = mesh(&cx);
        let material = law(100.0 * MU, true);
        let motion = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::OriginalMaterial,
                &cx,
            ),
            &[DENSITY],
            &vec![[0.0; 3]; nodes.len()],
            MASS_BYTES,
            &cx,
        )
        .unwrap();
        let before = motion.solid().quadrature_states().to_vec();
        let mut queries = 0;
        let result = motion.advance_adaptive_pair_toward(
            0.02,
            0.02,
            |_| {
                queries += 1;
                if queries == 2 {
                    cancel.request();
                }
                Ok(loads(vec![]))
            },
            settings(),
            controls(),
            accuracy(),
            &cx,
        );
        assert!(matches!(
            result,
            Err(MixedViscoAdaptiveError::Native(MixedViscoTetError::Native(
                MixedUpError::Cancelled
            )))
        ));
        assert_eq!(queries, 2);
        assert_eq!(motion.solid().quadrature_states(), before.as_slice());
        assert_eq!(motion.solid().time_s(), 0.0);
        assert!(
            motion
                .solid()
                .displacement_m()
                .iter()
                .flatten()
                .all(|u| *u == 0.0)
        );
        assert!(motion.velocity_mps().iter().flatten().all(|v| *v == 0.0));
    });
}
