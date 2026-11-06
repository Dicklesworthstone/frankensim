//! Adaptive midpoint paths retain native work and irreversible half-step history.

use super::*;
use fs_solid::brittle::solid::dynamics::midpoint::time::BrittleTetMidpointTimeControls;
use fs_solid::brittle::solid::dynamics::time::{
    BrittleTetTimeError, BrittleTetTimeRefinementCause,
};

fn time_controls(duration_s: f64) -> BrittleTetMidpointTimeControls {
    BrittleTetMidpointTimeControls {
        initial_interval_s: duration_s,
        minimum_interval_s: duration_s / 4096.0,
        maximum_interval_s: duration_s,
        maximum_intervals: 256,
        maximum_trials: 2048,
        displacement_tolerance_m: 1e-5,
        velocity_tolerance_mps: 1e-5,
        damage_tolerance: 1e-4,
        fracture_dissipation_tolerance_j: 1e-10,
        relative_tolerance: 0.0,
        maximum_damage_increment: 0.08,
        work: controls(),
    }
}

fn x_force_interval(start_n: f64, end_n: f64) -> Vec<(usize, [f64; 2])> {
    (0..4)
        .map(|node| (3 * node, [start_n / 4.0, end_n / 4.0]))
        .collect()
}

#[test]
fn tighter_midpoint_time_error_improves_oscillator_phase_without_kinetic_damping() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let mut velocity = [[0.0; 3]; 4];
        velocity[1][0] = 0.002;
        let old = BrittleTetMotion::try_from_solid(&solid, &velocity, 0.0, 1024, cx).unwrap();
        let supports: Vec<_> = prescribed(0.0)
            .into_iter()
            .filter(|&(dof, _)| dof != 3)
            .collect();
        let duration = 0.05_f64;
        let loose_controls = BrittleTetMidpointTimeControls {
            displacement_tolerance_m: 1e-3,
            velocity_tolerance_mps: 0.1,
            ..time_controls(duration)
        };
        let loose = old
            .advance_midpoint_time(duration, &supports, &[], settings(), loose_controls, cx)
            .unwrap();
        let tight = old
            .advance_midpoint_time(
                duration,
                &supports,
                &[],
                settings(),
                BrittleTetMidpointTimeControls {
                    displacement_tolerance_m: 1e-8,
                    velocity_tolerance_mps: 1e-6,
                    ..loose_controls
                },
                cx,
            )
            .unwrap();

        // The only free DOF has the consistent diagonal mass rho*V/10,
        // not the row-sum lumped mass rho*V/4; its stiffness is E*V/L^2.
        let oscillator_mass = DENSITY * VOLUME / 10.0;
        let omega = (YOUNG * VOLUME / oscillator_mass).sqrt();
        let exact_u = velocity[1][0] * (omega * duration).sin() / omega;
        let exact_v = velocity[1][0] * (omega * duration).cos();
        let phase_error = |motion: &BrittleTetMotion<'_>| {
            let du = motion.solid().solution().displacement_m[1][0] - exact_u;
            let dv = (motion.velocity_mps()[1][0] - exact_v) / omega;
            du.hypot(dv)
        };
        assert_eq!(loose.records.len(), 1);
        assert!(tight.records.len() > loose.records.len());
        assert!(tight.rejected_intervals > 0);
        assert!(phase_error(&tight.motion) < 0.15 * phase_error(&loose.motion));
        let initial_energy = 0.5 * oscillator_mass * velocity[1][0].powi(2);
        near(old.kinetic_energy_j(), initial_energy, 2e-18);
        for path in [&loose, &tight] {
            assert_eq!(path.energy.inertia.numerical_loss_j, 0.0);
            assert_eq!(path.energy.fracture_dissipation_j, 0.0);
            assert_eq!(path.energy.external_mean_work_j(), 0.0);
            let mut previous = old.time_s();
            for record in &path.records {
                assert_eq!(record.start_time_s, previous);
                assert!(record.start_time_s < record.midpoint_time_s);
                assert!(record.midpoint_time_s < record.end_time_s);
                assert!(record.estimate.normalized_error <= 1.0);
                for half in &record.energy {
                    assert_eq!(half.inertia.numerical_loss_j, 0.0);
                    near(
                        half.inertia.kinetic_after_j + half.stored_after_j,
                        initial_energy,
                        2e-12,
                    );
                }
                previous = record.end_time_s;
            }
            assert_eq!(previous, duration);
            assert_eq!(path.motion.time_s(), duration);
            assert_eq!(path.motion.accepted_steps(), 2 * path.records.len());
            assert!(path.trials >= 3 * path.records.len());
            assert_eq!(path.motion.reference_mass_kg(), old.reference_mass_kg());
            assert_eq!(
                path.motion.solid().solution().element_states[0].damage(),
                0.0
            );
            near(
                path.motion.kinetic_energy_j() + path.motion.solid().stored_energy_j(),
                initial_energy,
                2e-12,
            );
        }
        assert_eq!(old.velocity_mps(), &velocity);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.accepted_steps(), 0);
    });
}

#[test]
fn affine_forces_use_each_actual_interval_mean_and_supports_reach_absolute_targets() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let velocity = [[0.03, -0.02, 0.01]; 4];
        let initial_time = 1.25;
        let old =
            BrittleTetMotion::try_from_solid(&solid, &velocity, initial_time, 1024, cx).unwrap();
        let duration = 0.1_f64;
        let force_start = 2.0;
        let force_end = 6.0;
        let forces = x_force_interval(force_start, force_end);
        let step = old
            .advance_midpoint_time(
                duration,
                &[],
                &forces,
                settings(),
                BrittleTetMidpointTimeControls {
                    displacement_tolerance_m: 2e-6,
                    ..time_controls(duration)
                },
                cx,
            )
            .unwrap();
        assert!(step.records.len() > 1);
        assert!(step.rejected_intervals > 0);
        let mass = old.reference_mass_kg();
        let acceleration_start = force_start / mass;
        let jerk = (force_end - force_start) / (mass * duration);
        let exact_velocity =
            velocity[0][0] + acceleration_start * duration + 0.5 * jerk * duration.powi(2);
        let exact_position = velocity[0][0] * duration
            + 0.5 * acceleration_start * duration.powi(2)
            + jerk * duration.powi(3) / 6.0;
        // Mean acceleration integrates affine force exactly in velocity.
        // Each native midpoint position has excess jerk*h^3/12, so this
        // formula also detects sampling loads at substep ENDPOINTS.
        let sum_h_cubed: f64 = step
            .records
            .iter()
            .map(|record| {
                (record.midpoint_time_s - record.start_time_s).powi(3)
                    + (record.end_time_s - record.midpoint_time_s).powi(3)
            })
            .sum();
        let discrete_position = exact_position + jerk * sum_h_cubed / 12.0;
        for (u, v) in step
            .motion
            .solid()
            .solution()
            .displacement_m
            .iter()
            .zip(step.motion.velocity_mps())
        {
            near(u[0], discrete_position, 2e-10);
            near(v[0], exact_velocity, 2e-10);
            for axis in 1..3 {
                near(u[axis], velocity[0][axis] * duration, 2e-10);
                near(v[axis], velocity[0][axis], 2e-10);
            }
        }
        near(
            step.energy.dead_load_mean_work_j,
            0.5 * mass * (exact_velocity.powi(2) - velocity[0][0].powi(2)),
            2e-10,
        );
        assert_eq!(step.energy.inertia.numerical_loss_j, 0.0);
        assert_eq!(step.energy.support_mean_work_j, 0.0);
        assert_eq!(step.motion.time_s(), initial_time + duration);
        near(step.motion.solid().stored_energy_j(), 0.0, 2e-12);
        near(step.motion.solid().fracture_dissipation_j(), 0.0, 2e-12);

        // A compatible imposed rigid velocity has no kinetic change.
        // Complete mean reactions therefore cancel the applied mean work.
        let travel = 0.02;
        let supported_velocity = [[travel / duration, 0.0, 0.0]; 4];
        let supported =
            BrittleTetMotion::try_from_solid(&solid, &supported_velocity, initial_time, 8, cx)
                .unwrap();
        let targets: Vec<_> = (0..12)
            .map(|dof| (dof, if dof % 3 == 0 { travel } else { 0.0 }))
            .collect();
        let moved = supported
            .advance_midpoint_time(
                duration,
                &targets,
                &forces,
                settings(),
                time_controls(duration),
                cx,
            )
            .unwrap();
        assert_eq!(moved.records.len(), 1);
        assert_eq!(moved.trials, 3);
        assert_eq!(moved.motion.accepted_steps(), 2);
        for (u, v) in moved
            .motion
            .solid()
            .solution()
            .displacement_m
            .iter()
            .zip(moved.motion.velocity_mps())
        {
            assert_eq!(*u, [travel, 0.0, 0.0]);
            near(v[0], travel / duration, 2e-12);
            assert_eq!(v[1..], [0.0; 2]);
        }
        let mean_force_work = travel * 0.5 * (force_start + force_end);
        near(moved.energy.dead_load_mean_work_j, mean_force_work, 2e-12);
        near(moved.energy.support_mean_work_j, -mean_force_work, 2e-12);
        near(moved.energy.inertia.kinetic_change_j, 0.0, 2e-12);
        assert_eq!(moved.energy.inertia.numerical_loss_j, 0.0);
        assert_eq!(
            moved.motion.reference_mass_kg(),
            supported.reference_mass_kg()
        );
    });
}

#[test]
fn tensile_damage_limits_keep_native_half_steps_work_and_discarded_trial_telemetry() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let duration = 0.1;
        let extension = 0.0045_f64;
        let mut velocity = [[0.0; 3]; 4];
        velocity[1][0] = extension / duration;
        let old = BrittleTetMotion::try_from_solid(&solid, &velocity, 0.0, 1024, cx).unwrap();
        let targets = prescribed(extension);
        let force_start = 0.3;
        let force_end = 0.9;
        let forces = [(3, [force_start, force_end])];
        let controls = time_controls(duration);
        let accepted = old
            .advance_midpoint_time(duration, &targets, &forces, settings(), controls, cx)
            .unwrap();
        assert!(accepted.records.len() > 1);
        assert!(accepted.rejected_intervals > 0);
        assert!(accepted.records[0].end_time_s < duration);

        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * GC / (STRENGTH * sources[0].law().band_width_m());
        let stress = STRENGTH * (failure - extension) / (failure - onset);
        let work = VOLUME
            * (0.5 * STRENGTH * onset
                + STRENGTH / (failure - onset)
                    * (failure * (extension - onset) - 0.5 * (extension.powi(2) - onset.powi(2))));
        let stored = 0.5 * stress * extension * VOLUME;
        near(accepted.energy.material_work_j, work, 2e-12);
        near(accepted.energy.mean_internal_work_j, work, 2e-12);
        near(accepted.motion.solid().stored_energy_j(), stored, 2e-12);
        near(accepted.energy.fracture_dissipation_j, work - stored, 2e-12);
        near(
            accepted.motion.solid().solution().element_states[0].damage(),
            1.0 - stress / (YOUNG * extension),
            2e-12,
        );
        near(
            accepted.motion.solid().solution().element_stress[0][0],
            stress,
            2e-11,
        );
        near(
            accepted.energy.dead_load_mean_work_j,
            extension * 0.5 * (force_start + force_end),
            2e-12,
        );
        near(accepted.energy.external_mean_work_j(), work, 2e-11);
        near(accepted.energy.inertia.kinetic_change_j, 0.0, 2e-12);
        assert_eq!(accepted.energy.inertia.numerical_loss_j, 0.0);

        // Native replay of only accepted half steps must reproduce the exact
        // physical endpoint and every mean-work receipt, without extrapolation.
        let mut replay: Option<BrittleTetMotion<'_>> = None;
        let mut replay_work = 0.0;
        let mut absolute_mean_error = 0.0;
        let mut accepted_tangents = 0;
        let mut previous_damage = 0.0;
        for record in &accepted.records {
            assert!(record.estimate.normalized_error <= 1.0);
            assert!(record.maximum_damage_increment <= controls.maximum_damage_increment);
            for (half, end_time) in [record.midpoint_time_s, record.end_time_s]
                .into_iter()
                .enumerate()
            {
                let state = replay.as_ref().unwrap_or(&old);
                let fraction = end_time / duration;
                let target = if end_time == duration {
                    extension
                } else {
                    fraction * extension
                };
                let mean_fraction = 0.5 * (state.time_s() / duration) + 0.5 * fraction;
                let mean_force =
                    mean_fraction.mul_add(force_end, (1.0 - mean_fraction) * force_start);
                let step = state
                    .advance_midpoint(
                        end_time - state.time_s(),
                        &prescribed(target),
                        &[(3, mean_force)],
                        settings(),
                        controls.work,
                        cx,
                    )
                    .unwrap();
                assert_eq!(step.energy, record.energy[half]);
                let damage = step.motion.solid().solution().element_states[0].damage();
                assert!(damage >= previous_damage);
                assert!(damage - previous_damage <= controls.maximum_damage_increment);
                previous_damage = damage;
                replay_work += step.energy.material_work_j;
                absolute_mean_error += step.energy.mean_work_error_j.abs();
                accepted_tangents += step.numerical_tangent_evaluations;
                replay = Some(step.motion);
            }
        }
        let replay = replay.unwrap();
        assert_eq!(
            replay.solid().solution(),
            accepted.motion.solid().solution()
        );
        assert_eq!(replay.velocity_mps(), accepted.motion.velocity_mps());
        near(replay_work, accepted.energy.material_work_j, 2e-15);
        near(
            absolute_mean_error,
            accepted.absolute_mean_work_error_j,
            2e-15,
        );
        assert!(accepted.absolute_mean_work_error_j <= controls.work.material_energy_tolerance_j);

        // The initial rejected comparison completes one coarse and two fine
        // solves. Their fallback costs must remain in total trial telemetry.
        let coarse = old
            .advance_midpoint(
                duration,
                &targets,
                &[(3, 0.5 * (force_start + force_end))],
                settings(),
                controls.work,
                cx,
            )
            .unwrap();
        let mut discarded_tangents = coarse.numerical_tangent_evaluations;
        let mut rejected_fine: Option<BrittleTetMotion<'_>> = None;
        for fraction in [0.5, 1.0] {
            let state = rejected_fine.as_ref().unwrap_or(&old);
            let end_time = fraction * duration;
            let mean_fraction = 0.5 * (state.time_s() / duration) + 0.5 * fraction;
            let mean_force = mean_fraction.mul_add(force_end, (1.0 - mean_fraction) * force_start);
            let step = state
                .advance_midpoint(
                    end_time - state.time_s(),
                    &prescribed(fraction * extension),
                    &[(3, mean_force)],
                    settings(),
                    controls.work,
                    cx,
                )
                .unwrap();
            discarded_tangents += step.numerical_tangent_evaluations;
            rejected_fine = Some(step.motion);
        }
        assert!(discarded_tangents > 0);
        assert!(accepted.numerical_tangent_evaluations >= accepted_tangents + discarded_tangents);
        assert_eq!(accepted.motion.time_s(), duration);
        assert_eq!(accepted.motion.accepted_steps(), 2 * accepted.records.len());
        assert_eq!(accepted.motion.reference_mass_kg(), old.reference_mass_kg());
        assert!(core::ptr::eq(
            accepted.motion.solid().materials().as_ptr(),
            sources.as_ptr(),
        ));
        assert_eq!(old.solid().solution().element_states[0].damage(), 0.0);
        assert_eq!(old.solid().fracture_dissipation_j(), 0.0);
        assert_eq!(old.velocity_mps(), &velocity);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.accepted_steps(), 0);
    });
}

#[test]
fn incomplete_midpoint_paths_preserve_motion_history_budget_and_native_refusals() {
    with_cx(|cx, gate| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let old = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 1.0, 1024, cx).unwrap();
        let snapshot = old.solid().solution().clone();
        let duration = 0.1;
        let forces = x_force_interval(2.0, 6.0);
        let controls = BrittleTetMidpointTimeControls {
            displacement_tolerance_m: 0.01,
            velocity_tolerance_mps: 0.01,
            ..time_controls(duration)
        };
        for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                old.advance_midpoint_time(invalid, &[], &forces, settings(), controls, cx),
                Err(BrittleTetTimeError::Invalid(_))
            ));
        }
        assert!(matches!(
            old.advance_midpoint_time(duration, &[(12, 0.0)], &forces, settings(), controls, cx,),
            Err(BrittleTetTimeError::Invalid(_))
        ));
        let limited = BrittleTetMidpointTimeControls {
            initial_interval_s: 0.5 * duration,
            maximum_interval_s: 0.5 * duration,
            maximum_trials: 3,
            ..controls
        };
        // A pair was accepted internally, but the whole requested time path
        // remains incomplete. No partial clock, velocity or history escapes.
        assert!(matches!(
            old.advance_midpoint_time(duration, &[], &forces, settings(), limited, cx),
            Err(BrittleTetTimeError::WorkLimit {
                accepted_intervals: 1,
                trials: 3,
            })
        ));
        let unresolved = BrittleTetMidpointTimeControls {
            minimum_interval_s: duration,
            displacement_tolerance_m: 1e-12,
            ..controls
        };
        match old.advance_midpoint_time(duration, &[], &forces, settings(), unresolved, cx) {
            Err(BrittleTetTimeError::MinimumInterval {
                cause: BrittleTetTimeRefinementCause::ErrorEstimate(error),
                ..
            }) => assert!(error.normalized_error > 1.0),
            other => panic!("expected unresolved midpoint time error, got {other:?}"),
        }
        let one_step =
            BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 1.0, 1, cx).unwrap();
        assert!(matches!(
            one_step.advance_midpoint_time(duration, &[], &forces, settings(), controls, cx),
            Err(BrittleTetTimeError::Endpoint(
                BrittleTetDynamicError::StepBudget { .. }
            ))
        ));

        // An unavailable analytic constitutive tangent is the same native
        // refusal at either entry point; refinement must not relabel it.
        let mut strict = controls;
        strict.work.difference_step_strain = None;
        let strict_targets = prescribed(0.0002);
        let native_refusal = old
            .advance_midpoint(
                (old.time_s() + duration) - old.time_s(),
                &strict_targets,
                &[],
                settings(),
                strict.work,
                cx,
            )
            .unwrap_err();
        assert_eq!(
            old.advance_midpoint_time(duration, &strict_targets, &[], settings(), strict, cx)
                .unwrap_err(),
            BrittleTetTimeError::Endpoint(native_refusal),
        );
        gate.request();
        assert!(matches!(
            old.advance_midpoint_time(duration, &[], &forces, settings(), controls, cx),
            Err(BrittleTetTimeError::Endpoint(
                BrittleTetDynamicError::Solid(BrittleTetError::Cancelled)
            ))
        ));
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.velocity_mps(), &[[0.0; 3]; 4]);
        assert_eq!(old.time_s(), 1.0);
        assert_eq!(old.accepted_steps(), 0);
        assert_eq!(old.remaining_steps(), 1024);
        assert_eq!(old.kinetic_energy_j(), 0.0);
        assert_eq!(old.solid().fracture_dissipation_j(), 0.0);
        assert_eq!(one_step.accepted_steps(), 0);
        let retry_gate = CancelGate::new_clock_free();
        let retry = Cx::new(
            &retry_gate,
            cx.arena(),
            cx.stream_key(),
            cx.budget(),
            cx.mode(),
        );
        let completed = old
            .advance_midpoint_time(duration, &[], &forces, settings(), controls, &retry)
            .unwrap();
        assert_eq!(completed.motion.time_s(), 1.0 + duration);
        assert_eq!(completed.motion.accepted_steps(), 2);
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.time_s(), 1.0);
    });
}
