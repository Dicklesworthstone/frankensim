use super::*;
use fs_solid::ductile::evolution::dynamics::{DuctileMotionControls, DuctileTetMotion};

fn motion_controls() -> DuctileMotionControls {
    DuctileMotionControls {
        continuation: controls(),
        maximum_numerical_kinetic_loss_j: 0.1,
    }
}

#[test]
fn accelerated_free_body_has_exact_mass_momentum_and_separate_numerical_loss() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m());
        let materials = [&material];
        // rho*V = 1 kg. Each P1 mass row sums to 1/4 kg; these four
        // forces therefore produce exactly uniform acceleration a_x=1 m/s2.
        let forces = [(0, 0.25), (3, 0.25), (6, 0.25), (9, 0.25)];
        let specimen = DuctileTetEvolution::new(
            &geometry,
            &materials,
            &bands,
            &[],
            &forces,
            TetAssemblyBudget::standard(),
            IncrementSettings::default(),
            cx,
        )
        .unwrap();
        let motion = DuctileTetMotion::new(&specimen, &[6.0], cx).unwrap();
        let virgin = [material.initial_state()];
        let step = motion
            .advance_from(
                &[0.0; 12],
                &[[0.0; 3]; 4],
                &virgin,
                [1.0, 1.0],
                conditions(0.1),
                motion_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(step.records.len(), 1);
        for node in 0..4 {
            near(step.mechanics.velocity_mps[node][0], 0.1);
            near(step.mechanics.solution.displacement_m[node][0], 0.01);
            for axis in 1..3 {
                near(step.mechanics.velocity_mps[node][axis], 0.0);
                near(step.mechanics.solution.displacement_m[node][axis], 0.0);
            }
        }
        near(step.energy.inertia.reference_mass_kg, 1.0);
        near(step.energy.inertia.kinetic_change_j, 0.005);
        near(step.energy.inertia.numerical_loss_j, 0.005);
        near(step.energy.external_work_j, 0.01);
        near(step.energy.material.external_work_j, 0.0);
        near(step.energy.material.plastic_dissipation_j, 0.0);
        near(step.energy.material.damage_release_j, 0.0);
        near(step.energy.balance_error_j, 0.0);
        assert_eq!(step.mechanics.solution.element_states[0].damage(), 0.0);

        // The same native motion is refined to satisfy a physically separate
        // kinetic numerical-loss budget. Refinement preserves total momentum.
        let refined = motion
            .advance_from(
                &[0.0; 12],
                &[[0.0; 3]; 4],
                &virgin,
                [1.0, 1.0],
                conditions(0.1),
                DuctileMotionControls {
                    maximum_numerical_kinetic_loss_j: 1e-4,
                    ..motion_controls()
                },
                cx,
            )
            .unwrap();
        assert!(refined.rejected_trials > 0 && refined.records.len() > 1);
        assert!(refined.energy.inertia.numerical_loss_j <= 1e-4);
        near(refined.energy.inertia.kinetic_change_j, 0.005);
        for velocity in &refined.mechanics.velocity_mps {
            near(velocity[0], 0.1);
        }
        near(refined.energy.material.damage_release_j, 0.0);
        near(refined.energy.balance_error_j, 0.0);
        assert_eq!(virgin[0].damage(), 0.0);
    });
}

#[test]
fn free_shear_degree_solves_momentum_with_active_ductile_softening() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap()];
        let calibration = model(bands[0].band_width_m());
        // This oracle is one backward-Euler interval. Explicitly admit its
        // damage increment, so subdivision cannot change the tested equation.
        let material = DuctileJ2::new(
            calibration.native().clone(),
            calibration.initiation_law().clone(),
            DuctileJ2Settings {
                maximum_damage_increment: 1.0,
                ..calibration.settings()
            },
        )
        .unwrap();
        let materials = [&material];
        let supports: Vec<_> = (0..12)
            .filter(|&dof| dof != 6)
            .map(|dof| (dof, 0.0))
            .collect();
        let duration_s = 0.01;
        let target_u = 0.012;
        let volume = bands[0].reference_volume_m3();
        let reference_mass = 6.0 * volume;
        let diagonal_mass = reference_mass / 10.0;
        let tau = shear_stress(0.5 * target_u, bands[0].band_width_m());
        // Only node 2's x displacement is free: eps_xy=u/2, f_int=V*tau,
        // and the original consistent reference mass gives M_66=rho*V/10.
        let force = diagonal_mass * target_u / (duration_s * duration_s) + volume * tau;
        let forces = [(6, force)];
        let specimen = DuctileTetEvolution::new(
            &geometry,
            &materials,
            &bands,
            &supports,
            &forces,
            TetAssemblyBudget::standard(),
            IncrementSettings::default(),
            cx,
        )
        .unwrap();
        let motion = DuctileTetMotion::new(&specimen, &[6.0], cx).unwrap();
        let step = motion
            .advance_from(
                &[0.0; 12],
                &[[0.0; 3]; 4],
                &[material.initial_state()],
                [0.0, 1.0],
                conditions(duration_s),
                motion_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(step.records.len(), 1);
        assert!(step.mechanics.solution.element_states[0].damage() > 0.5);
        near(step.mechanics.solution.displacement_m[2][0], target_u);
        near(step.mechanics.velocity_mps[2][0], target_u / duration_s);
        near(step.mechanics.solution.element_stress[0][3], tau);
        near(
            step.mechanics.solution.element_states[0]
                .native_state()
                .alpha,
            alpha_at(0.5 * target_u),
        );
        let kinetic = 0.5 * diagonal_mass * (target_u / duration_s).powi(2);
        near(step.energy.inertia.kinetic_change_j, kinetic);
        near(step.energy.inertia.numerical_loss_j, kinetic);
        near(step.energy.external_work_j, force * target_u);
        near(
            step.energy.material.external_work_j,
            volume * tau * target_u,
        );
        assert!(step.energy.material.plastic_dissipation_j > 0.0);
        // The endpoint update releases damage from old recoverable storage;
        // this virgin start has none. Damage itself is already active.
        near(step.energy.material.damage_release_j, 0.0);
        near(step.energy.balance_error_j, 0.0);

        // Constrained rows retain their consistent-mass coupling. Together
        // applied and support forces advance the full body's P1 momentum.
        let reactions = &step.mechanics.solution.increments[0].reactions_n;
        let net_x =
            supports.iter().zip(reactions).fold(
                force,
                |sum, (&(dof, _), &r)| {
                    if dof % 3 == 0 { sum + r } else { sum }
                },
            );
        near(
            net_x,
            0.25 * reference_mass * target_u / (duration_s * duration_s),
        );
    });
}

#[test]
fn imposed_motion_crosses_ductile_initiation_without_turning_inertia_into_heat() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m());
        let materials = [&material];
        let supports = prescribed(&NODES);
        // Loading a prescribed DOF changes its support reaction; applied and
        // support work together must still equal the native material work.
        let supported_force = [(6, 7.0)];
        let specimen = DuctileTetEvolution::new(
            &geometry,
            &materials,
            &bands,
            &supports,
            &supported_force,
            TetAssemblyBudget::standard(),
            IncrementSettings::default(),
            cx,
        )
        .unwrap();
        let motion = DuctileTetMotion::new(&specimen, &[6.0], cx).unwrap();
        let virgin = [material.initial_state()];
        // The prescribed affine path starts with its compatible velocity.
        // No impulsive start contaminates the constitutive work comparison.
        let velocity = [[0.0; 3], [0.0; 3], [0.008, 0.0, 0.0], [0.0; 3]];
        let step = motion
            .advance_from(
                &[0.0; 12],
                &velocity,
                &virgin,
                [0.0, 0.012],
                conditions(3.0),
                motion_controls(),
                cx,
            )
            .unwrap();
        assert!(step.rejected_trials > 0);
        assert!(step.mechanics.solution.element_states[0].damage() > 0.85);
        near(
            step.mechanics.solution.element_states[0]
                .native_state()
                .alpha,
            alpha_at(0.012),
        );
        near(
            step.mechanics.solution.element_stress[0][3],
            shear_stress(0.012, bands[0].band_width_m()),
        );
        near(step.energy.inertia.kinetic_change_j, 0.0);
        near(step.energy.inertia.numerical_loss_j, 0.0);
        assert!(step.energy.material.plastic_dissipation_j > 0.0);
        assert!(step.energy.material.damage_release_j > 0.0);
        near(
            step.energy.external_work_j,
            step.energy.material.stored_change_j
                + step.energy.material.plastic_dissipation_j
                + step.energy.material.damage_release_j
                + step.energy.material.numerical_loss_j
                + step.energy.inertia.kinetic_change_j
                + step.energy.inertia.numerical_loss_j,
        );
        near(step.energy.balance_error_j, 0.0);
        let mut old_damage = 0.0;
        let mut old_fraction = 0.0;
        let mut old_factor = 0.0;
        for record in &step.records {
            assert!(record.maximum_damage - old_damage <= DAMAGE_LIMIT + 1e-13);
            near(record.load_factor, record.fraction * 0.012);
            near(record.duration_s, 3.0 * (record.fraction - old_fraction));
            near(
                record.energy.material.external_work_j,
                2.0 * bands[0].reference_volume_m3()
                    * shear_stress(record.load_factor, bands[0].band_width_m())
                    * (record.load_factor - old_factor),
            );
            old_damage = record.maximum_damage;
            old_fraction = record.fraction;
            old_factor = record.load_factor;
        }
        assert_eq!(old_fraction, 1.0);
        near(step.records.iter().map(|r| r.duration_s).sum(), 3.0);
        assert_eq!(virgin[0].damage(), 0.0);
    });
}

#[test]
fn refused_or_cancelled_motion_keeps_full_state_and_replays_identically() {
    with_cx(|gate, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m());
        let materials = [&material];
        let supports = prescribed(&NODES);
        let specimen = DuctileTetEvolution::new(
            &geometry,
            &materials,
            &bands,
            &supports,
            &[],
            TetAssemblyBudget::standard(),
            IncrementSettings::default(),
            cx,
        )
        .unwrap();
        assert!(DuctileTetMotion::new(&specimen, &[0.0], cx).is_err());
        let motion = DuctileTetMotion::new(&specimen, &[6.0], cx).unwrap();
        let old_u = [0.0; 12];
        let old_v = [[0.0; 3], [0.0; 3], [0.008, 0.0, 0.0], [0.0; 3]];
        let old_states = [material.initial_state()];
        let saved = old_states.clone();
        // This rigid shift has the same valid virgin constitutive history,
        // but the declared support path starts at zero. Do not invent an
        // initial displacement jump or an associated inertial impulse.
        let shifted = [0.01, 0.0, 0.0].repeat(4);
        assert!(matches!(
            motion.advance_from(
                &shifted,
                &old_v,
                &old_states,
                [0.0, 0.012],
                conditions(3.0),
                motion_controls(),
                cx,
            ),
            Err(DuctileTetEvolutionError::Invalid(_))
        ));
        let mut bounded = motion_controls();
        bounded.continuation.maximum_accepted_steps = 1;
        assert!(matches!(
            motion.advance_from(
                &old_u,
                &old_v,
                &old_states,
                [0.0, 0.012],
                conditions(3.0),
                bounded,
                cx,
            ),
            Err(DuctileTetEvolutionError::WorkLimit {
                accepted_steps: 1,
                ..
            })
        ));
        assert_eq!(old_states, saved);
        let accepted = motion
            .advance_from(
                &old_u,
                &old_v,
                &old_states,
                [0.0, 0.012],
                conditions(3.0),
                motion_controls(),
                cx,
            )
            .unwrap();
        gate.request();
        assert!(matches!(
            motion.advance_from(
                &old_u,
                &old_v,
                &old_states,
                [0.0, 0.012],
                conditions(3.0),
                motion_controls(),
                cx,
            ),
            Err(DuctileTetEvolutionError::Cancelled)
        ));
        let restarted = with_cx(|_, fresh| {
            motion.advance_from(
                &old_u,
                &old_v,
                &old_states,
                [0.0, 0.012],
                conditions(3.0),
                motion_controls(),
                fresh,
            )
        })
        .unwrap();
        assert_eq!(accepted.mechanics.solution, restarted.mechanics.solution);
        assert_eq!(
            accepted.mechanics.velocity_mps,
            restarted.mechanics.velocity_mps
        );
        assert_eq!(accepted.records, restarted.records);
        assert_eq!(accepted.energy, restarted.energy);
        assert_eq!(old_states, saved);
        assert_eq!(old_u, [0.0; 12]);
    });
}

#[test]
fn native_deformation_and_retained_history_must_match_before_any_unloading() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m());
        let materials = [&material];
        let supports = prescribed(&NODES);
        let specimen = DuctileTetEvolution::new(
            &geometry,
            &materials,
            &bands,
            &supports,
            &[],
            TetAssemblyBudget::standard(),
            IncrementSettings::default(),
            cx,
        )
        .unwrap();
        let virgin = [material.initial_state()];
        let prefix = specimen
            .advance_from(
                &[0.0; 12],
                &virgin,
                [0.0, 0.003],
                conditions(1.0),
                controls(),
                cx,
            )
            .unwrap();
        let old_u = prefix.mechanics.displacement_m.as_flattened();
        let old_states = &prefix.mechanics.element_states;
        // Zero endpoint stress can otherwise make the work mismatch disappear:
        // unloading an already deformed body must not erase its initial state.
        assert!(matches!(
            specimen.advance_from(
                old_u,
                &virgin,
                [0.003, 0.0],
                conditions(1.0),
                controls(),
                cx,
            ),
            Err(DuctileTetEvolutionError::Invalid(_))
        ));
        let held = specimen
            .advance_from(
                old_u,
                old_states,
                [0.003, 0.003],
                conditions(1.0),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(held.mechanics.element_states, *old_states);
        assert_eq!(
            held.mechanics.displacement_m,
            prefix.mechanics.displacement_m
        );
        near(held.energy.external_work_j, 0.0);

        let motion = DuctileTetMotion::new(&specimen, &[6.0], cx).unwrap();
        assert!(matches!(
            motion.advance_from(
                old_u,
                &[[0.0; 3]; 4],
                &virgin,
                [0.003, 0.0],
                conditions(1.0),
                motion_controls(),
                cx,
            ),
            Err(DuctileTetEvolutionError::Invalid(_))
        ));
        let held_motion = motion
            .advance_from(
                old_u,
                &[[0.0; 3]; 4],
                old_states,
                [0.003, 0.003],
                conditions(1.0),
                motion_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(held_motion.mechanics.solution.element_states, *old_states);
        assert_eq!(
            held_motion.mechanics.solution.displacement_m,
            prefix.mechanics.displacement_m
        );
        assert_eq!(held_motion.mechanics.velocity_mps, [[0.0; 3]; 4]);
        near(held_motion.energy.external_work_j, 0.0);
        assert_eq!(virgin[0].damage(), 0.0);
    });
}
