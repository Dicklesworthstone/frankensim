//! Native step doubling controls physical motion without extrapolating damage.
use super::*;
use fs_solid::brittle::solid::dynamics::time::{
    BrittleTetTimeControls, BrittleTetTimeError, BrittleTetTimeRefinementCause,
};

fn time_controls(duration_s: f64) -> BrittleTetTimeControls {
    BrittleTetTimeControls {
        initial_interval_s: duration_s,
        minimum_interval_s: duration_s / 4096.0,
        maximum_interval_s: duration_s,
        maximum_intervals: 256,
        maximum_trials: 2048,
        displacement_tolerance_m: 1e-7,
        velocity_tolerance_mps: 1e-6,
        damage_tolerance: 0.02,
        fracture_dissipation_tolerance_j: 1e-8,
        relative_tolerance: 0.0,
        maximum_damage_increment: 0.08,
        work: dynamic_controls(),
    }
}

fn x_force_interval(start_n: f64, end_n: f64) -> Vec<(usize, [f64; 2])> {
    (0..4)
        .map(|node| (3 * node, [start_n / 4.0, end_n / 4.0]))
        .collect()
}

#[test]
fn tighter_time_error_reduces_free_motion_error_and_numerical_kinetic_loss() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = dynamic_card(2000.0);
        let materials = [dynamic_material(&source, &geometry)];
        let solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let velocity = [[0.2, -0.1, 0.05]; 4];
        let old = BrittleTetMotion::try_from_solid(&solid, &velocity, 0.0, 1024, cx).unwrap();
        let duration = 1e-5_f64;
        let total_force = 2.0;
        let forces = x_force_interval(total_force, total_force);
        let loose_controls = time_controls(duration);
        let tight_controls = BrittleTetTimeControls {
            displacement_tolerance_m: 5e-10,
            ..loose_controls
        };
        let loose = old
            .advance_time(
                duration,
                &[],
                &forces,
                IncrementSettings::default(),
                loose_controls,
                cx,
            )
            .unwrap();
        let tight = old
            .advance_time(
                duration,
                &[],
                &forces,
                IncrementSettings::default(),
                tight_controls,
                cx,
            )
            .unwrap();
        let acceleration = total_force / old.reference_mass_kg();
        let exact_position = velocity[0][0] * duration + 0.5 * acceleration * duration.powi(2);
        let error = |motion: &BrittleTetMotion<'_>| {
            motion
                .solid()
                .solution()
                .displacement_m
                .iter()
                .map(|u| (u[0] - exact_position).abs())
                .fold(0.0_f64, f64::max)
        };
        assert_eq!(loose.records.len(), 1);
        assert!(tight.records.len() > loose.records.len());
        assert!(tight.rejected_intervals > 0);
        assert!(error(&tight.motion) < 0.4 * error(&loose.motion));
        assert!(
            tight.energy.inertia.numerical_loss_j < 0.4 * loose.energy.inertia.numerical_loss_j
        );

        // For arbitrary accepted BE timesteps, the exact discrete free-body
        // solution is u = v0*T + a*(T^2 + sum(h^2))/2. Its excess displacement
        // over the analytic trajectory and numerical kinetic loss both shrink.
        let mut sum_h_squared = 0.0;
        let mut previous_time = 0.0;
        for record in &tight.records {
            assert_eq!(record.start_time_s, previous_time);
            assert!(record.start_time_s < record.midpoint_time_s);
            assert!(record.midpoint_time_s < record.end_time_s);
            assert!(record.estimate.normalized_error <= 1.0);
            sum_h_squared += (record.midpoint_time_s - record.start_time_s).powi(2)
                + (record.end_time_s - record.midpoint_time_s).powi(2);
            previous_time = record.end_time_s;
        }
        assert_eq!(previous_time, duration);
        let discrete_position = exact_position + 0.5 * acceleration * sum_h_squared;
        for (u, v) in tight
            .motion
            .solid()
            .solution()
            .displacement_m
            .iter()
            .zip(tight.motion.velocity_mps())
        {
            close(u[0], discrete_position, 2e-8);
            close(v[0], velocity[0][0] + acceleration * duration, 2e-8);
            for axis in 1..3 {
                close(u[axis], velocity[0][axis] * duration, 2e-8);
                close(v[axis], velocity[0][axis], 2e-8);
            }
        }
        close(
            tight.energy.inertia.numerical_loss_j,
            0.5 * old.reference_mass_kg() * acceleration.powi(2) * sum_h_squared,
            2e-6,
        );
        assert_eq!(tight.motion.time_s(), duration);
        assert_eq!(tight.motion.accepted_steps(), 2 * tight.records.len());
        assert!(tight.trials >= 3 * tight.records.len());
        assert_eq!(
            tight.motion.solid().solution().element_states[0].damage(),
            0.0
        );
        assert!(tight.energy.fracture_dissipation_j.abs() < 1e-12);
        assert!(tight.energy.total_balance_error_j.abs() < 1e-9);
        assert_eq!(tight.motion.reference_mass_kg(), old.reference_mass_kg());
        assert_eq!(old.velocity_mps(), &velocity);
        assert_eq!(old.solid().solution().displacement_m, vec![[0.0; 3]; 4]);
        assert_eq!(old.accepted_steps(), 0);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn damage_limited_time_path_retains_native_half_steps_and_exact_fracture_work() {
    with_cx(false, |cx| {
        let length = 0.03;
        let nodes = points(length);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = dynamic_card(2000.0);
        let materials = [dynamic_material(&source, &geometry)];
        let solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let duration = 1e-4;
        let final_strain = 0.001;
        let final_displacement = length * final_strain;
        let targets: Vec<_> = (0..12)
            .map(|dof| (dof, if dof == 3 { final_displacement } else { 0.0 }))
            .collect();
        // Match the imposed constant strain rate at the initial endpoint, so
        // the test's fracture loss is independent of a support velocity impulse.
        let mut velocity = [[0.0; 3]; 4];
        velocity[1][0] = final_displacement / duration;
        let old = BrittleTetMotion::try_from_solid(&solid, &velocity, 0.0, 1024, cx).unwrap();
        let controls = time_controls(duration);
        let accepted = old
            .advance_time(
                duration,
                &targets,
                &[],
                IncrementSettings::default(),
                controls,
                cx,
            )
            .unwrap();
        assert!(accepted.records.len() > 1);
        assert!(accepted.rejected_intervals > 0);

        let volume = geometry.element_volumes_m3()[0];
        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * FRACTURE_WORK / (STRENGTH * materials[0].law().band_width_m());
        let stress = STRENGTH * (failure - final_strain) / (failure - onset);
        let work = volume
            * (0.5 * STRENGTH * onset
                + STRENGTH / (failure - onset)
                    * (failure * (final_strain - onset)
                        - 0.5 * (final_strain.powi(2) - onset.powi(2))));
        let storage = 0.5 * stress * final_strain * volume;
        close(accepted.energy.material_work_j, work, 2e-8);
        close(accepted.motion.solid().stored_energy_j(), storage, 2e-8);
        close(accepted.energy.fracture_dissipation_j, work - storage, 2e-8);
        close(
            accepted.motion.solid().solution().element_states[0].damage(),
            1.0 - stress / (YOUNG * final_strain),
            2e-8,
        );
        assert!(accepted.energy.fracture_dissipation_j > 0.0);

        // Replay only returned physical half-step endpoints through the public
        // native owner. Extrapolated or coarse histories cannot match this path.
        let mut replay: Option<BrittleTetMotion<'_>> = None;
        let mut native_work = 0.0;
        let mut absolute_defect = 0.0;
        for record in &accepted.records {
            assert!(record.maximum_damage_increment <= controls.maximum_damage_increment);
            assert!(record.estimate.normalized_error <= 1.0);
            for (half, time) in [record.midpoint_time_s, record.end_time_s]
                .into_iter()
                .enumerate()
            {
                let state = replay.as_ref().unwrap_or(&old);
                let fraction = time / duration;
                let prescribed: Vec<_> = targets
                    .iter()
                    .map(|&(dof, target)| {
                        (
                            dof,
                            if time == duration {
                                target
                            } else {
                                fraction * target
                            },
                        )
                    })
                    .collect();
                let step = state
                    .advance(
                        time - state.time_s(),
                        &prescribed,
                        &[],
                        IncrementSettings::default(),
                        dynamic_controls(),
                        cx,
                    )
                    .unwrap();
                assert_eq!(step.energy, record.energy[half]);
                native_work += step.energy.material_work_j;
                absolute_defect += step.energy.constitutive_endpoint_defect_j.abs();
                replay = Some(step.motion);
            }
        }
        let replay = replay.unwrap();
        assert_eq!(
            replay.solid().solution(),
            accepted.motion.solid().solution()
        );
        assert_eq!(replay.velocity_mps(), accepted.motion.velocity_mps());
        assert_eq!(replay.time_s(), duration);
        close(native_work, accepted.energy.material_work_j, 1e-12);
        close(
            absolute_defect,
            accepted.absolute_constitutive_defect_j,
            1e-12,
        );
        assert!(
            accepted.absolute_constitutive_defect_j
                <= controls.work.work.maximum_endpoint_work_defect_j
        );
        assert!(
            accepted.energy.inertia.numerical_loss_j
                <= controls.work.maximum_numerical_kinetic_loss_j
        );
        assert!(accepted.energy.inertia.numerical_loss_j < 1e-15);
        assert!(accepted.energy.inertia.kinetic_change_j.abs() < 1e-12);
        assert!(accepted.energy.total_balance_error_j.abs() < 1e-9);
        assert_eq!(accepted.motion.accepted_steps(), 2 * accepted.records.len());
        assert_eq!(old.solid().solution().element_states[0].damage(), 0.0);
        assert_eq!(old.solid().fracture_dissipation_j(), 0.0);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn moving_supports_reach_exact_targets_and_keep_applied_and_kinetic_work_separate() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = dynamic_card(2000.0);
        let materials = [dynamic_material(&source, &geometry)];
        let solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let old = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 8, cx).unwrap();
        let duration = 1e-4_f64;
        let travel = 1e-5;
        let supports: Vec<_> = (0..12)
            .map(|dof| (dof, if dof % 3 == 0 { travel } else { 0.0 }))
            .collect();
        let force_start = 2.0;
        let force_end = 4.0;
        let step = old
            .advance_time(
                duration,
                &supports,
                &x_force_interval(force_start, force_end),
                IncrementSettings::default(),
                time_controls(duration),
                cx,
            )
            .unwrap();
        assert_eq!(step.records.len(), 1);
        assert_eq!(step.trials, 3);
        assert_eq!(step.motion.accepted_steps(), 2);
        assert_eq!(step.motion.time_s(), duration);
        for (u, v) in step
            .motion
            .solid()
            .solution()
            .displacement_m
            .iter()
            .zip(step.motion.velocity_mps())
        {
            assert_eq!(*u, [travel, 0.0, 0.0]);
            close(v[0], travel / duration, 2e-12);
            assert_eq!(v[1..], [0.0; 2]);
        }
        let kinetic = 0.5 * old.reference_mass_kg() * (travel / duration).powi(2);
        let midpoint_force = 0.5 * (force_start + force_end);
        let dead_work = 0.5 * travel * (midpoint_force + force_end);
        close(step.energy.inertia.kinetic_after_j, kinetic, 2e-12);
        close(step.energy.inertia.kinetic_change_j, kinetic, 2e-12);
        close(step.energy.inertia.numerical_loss_j, kinetic, 2e-12);
        close(step.energy.inertia.inertial_work_j, 2.0 * kinetic, 2e-12);
        close(step.energy.dead_load_endpoint_work_j, dead_work, 2e-12);
        close(
            step.energy.support_endpoint_work_j,
            2.0 * kinetic - dead_work,
            2e-12,
        );
        assert!(step.energy.support_endpoint_work_j < 0.0);
        assert!(step.energy.internal_endpoint_work_j.abs() < 1e-12);
        assert!(step.energy.fracture_dissipation_j.abs() < 1e-12);
        assert!(step.energy.total_balance_error_j.abs() < 1e-12);
        assert_eq!(step.motion.reference_mass_kg(), old.reference_mass_kg());
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn incomplete_time_requests_leave_the_original_motion_and_history_unchanged() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = dynamic_card(2000.0);
        let materials = [dynamic_material(&source, &geometry)];
        let solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let old = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 1024, cx).unwrap();
        let snapshot = old.solid().solution().clone();
        let duration = 1e-5;
        let forces = x_force_interval(2.0, 2.0);
        let controls = time_controls(duration);
        for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                old.advance_time(
                    invalid,
                    &[],
                    &forces,
                    IncrementSettings::default(),
                    controls,
                    cx
                ),
                Err(BrittleTetTimeError::Invalid(_))
            ));
        }
        assert!(matches!(
            old.advance_time(
                duration,
                &[(12, 0.0)],
                &forces,
                IncrementSettings::default(),
                controls,
                cx,
            ),
            Err(BrittleTetTimeError::Invalid(_))
        ));
        let limited = BrittleTetTimeControls {
            initial_interval_s: 0.5 * duration,
            maximum_interval_s: 0.5 * duration,
            maximum_trials: 3,
            ..controls
        };
        // One complete native pair is solved internally, but the requested
        // duration is incomplete. That partial motion must not escape.
        assert!(matches!(
            old.advance_time(
                duration,
                &[],
                &forces,
                IncrementSettings::default(),
                limited,
                cx
            ),
            Err(BrittleTetTimeError::WorkLimit {
                accepted_intervals: 1,
                trials: 3
            })
        ));
        let unresolved = BrittleTetTimeControls {
            minimum_interval_s: duration,
            displacement_tolerance_m: 1e-12,
            ..controls
        };
        match old.advance_time(
            duration,
            &[],
            &forces,
            IncrementSettings::default(),
            unresolved,
            cx,
        ) {
            Err(BrittleTetTimeError::MinimumInterval {
                cause: BrittleTetTimeRefinementCause::ErrorEstimate(error),
                ..
            }) => assert!(error.normalized_error > 1.0),
            other => panic!("expected resolved time-error refusal, got {other:?}"),
        }
        let one_step =
            BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 1, cx).unwrap();
        assert!(matches!(
            one_step.advance_time(
                duration,
                &[],
                &forces,
                IncrementSettings::default(),
                controls,
                cx,
            ),
            Err(BrittleTetTimeError::Endpoint(
                BrittleTetDynamicError::StepBudget { .. }
            ))
        ));
        with_cx(true, |cancelled| {
            assert!(
                old.advance_time(
                    duration,
                    &[],
                    &forces,
                    IncrementSettings::default(),
                    controls,
                    cancelled,
                )
                .is_err()
            );
        });
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.velocity_mps(), &[[0.0; 3]; 4]);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.accepted_steps(), 0);
        assert_eq!(old.kinetic_energy_j(), 0.0);
        assert_eq!(old.solid().fracture_dissipation_j(), 0.0);
        assert_eq!(one_step.accepted_steps(), 0);
        let accepted = old
            .advance_time(
                duration,
                &[],
                &forces,
                IncrementSettings::default(),
                controls,
                cx,
            )
            .unwrap();
        assert_eq!(accepted.motion.time_s(), duration);
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.time_s(), 0.0);
    });
}
