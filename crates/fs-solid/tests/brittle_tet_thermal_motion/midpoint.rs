//! Original mean forces, thermal work and retained physical endpoint history.

use super::*;
use fs_solid::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use fs_solid::brittle::solid::expansion::dynamics::midpoint::ThermalBrittleTetMidpointEnergy;

fn midpoint_controls() -> BrittleTetMidpointControls {
    BrittleTetMidpointControls {
        material_energy_tolerance_j: 1e-10,
        equilibrium_work_tolerance_j: 1e-10,
        kinetic_work_tolerance_j: 1e-10,
        total_energy_tolerance_j: 1e-10,
        difference_step_strain: Some(1e-7),
    }
}

fn check_midpoint_balances(energy: &ThermalBrittleTetMidpointEnergy) {
    for defect in [
        energy.strain_split_error_j,
        energy.material_balance_error_j,
        energy.mean_work_error_j,
        energy.equilibrium_work_residual_j,
        energy.inertia.balance_error_j,
        energy.total_balance_error_j,
    ] {
        assert!(defect.abs() < 1e-10, "{energy:?}");
    }
    assert_eq!(energy.inertia.numerical_loss_j, 0.0);
}

#[test]
fn cooling_midpoint_uses_analytic_mean_force_and_no_numerical_kinetic_loss() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = card("midpoint-cooling", 4.5e-5, false);
        let reference = [source(&card, &geometry, 400.0, 400.0)];
        let cold = [reference[0].at_temperature(&card, 300.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let old =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 8, cx).unwrap();
        let dt = 0.01_f64;
        let step = old
            .advance_midpoint(
                &cold,
                dt,
                &axial_supports(),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();

        // Original bilinear opening work w(e) has mean stress
        // c*(failure - e/2 - onset*failure/(2e)).
        // With u=e+free and v_next=2u/dt, native centered momentum is
        // (2m/dt^2)*u + V*w(e)/e = 0. Its positive quadratic root
        // independently determines the actual free displacement.
        let free = -4.5e-5 * 100.0;
        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * GC / (STRENGTH * (1.0 / 3.0));
        let slope = STRENGTH / (failure - onset);
        let mass = DENSITY * VOLUME / 10.0;
        let inertia = 2.0 * mass / dt.powi(2);
        let a = inertia - 0.5 * VOLUME * slope;
        let b = inertia * free + VOLUME * slope * failure;
        let c = -0.5 * VOLUME * slope * onset * failure;
        let strain = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
        let u = strain + free;
        let velocity = 2.0 * u / dt;
        let stress = slope * (failure - strain);
        let opening_work = 0.5 * STRENGTH * onset + 0.5 * (STRENGTH + stress) * (strain - onset);
        let mean_x = opening_work / strain;
        let transverse_storage = YOUNG * free * free;
        let stored = VOLUME * (0.5 * stress * strain + transverse_storage);
        let fracture = VOLUME * (opening_work - 0.5 * stress * strain);
        let mechanical_work = VOLUME * mean_x * u;
        let eigenstrain_work = -VOLUME * free * (mean_x - YOUNG * free);
        let kinetic = 0.5 * mass * velocity * velocity;
        assert!(strain > onset && strain < failure);
        assert!((mean_x - stress).abs() > 0.3);
        near(
            step.motion.solid().solution().displacement_m[1][0],
            u,
            2e-12,
        );
        near(step.motion.velocity_mps()[1][0], velocity, 2e-10);
        near(
            step.motion.solid().solution().element_strain[0][0],
            u,
            2e-12,
        );
        near(
            step.motion.solid().solution().element_states[0].strain()[0],
            strain,
            2e-12,
        );
        // The published solid contains PHYSICAL endpoint stress, not the
        // different mean force response that drove the momentum equation.
        near(
            step.motion.solid().solution().element_stress[0][0],
            stress,
            2e-10,
        );
        near(
            step.motion.solid().solution().element_states[0].damage(),
            1.0 - stress / (YOUNG * strain),
            2e-11,
        );
        near(step.motion.solid().stored_energy_j(), stored, 2e-12);
        near(
            step.motion.solid().fracture_dissipation_j(),
            fracture,
            2e-12,
        );
        near(step.energy.fracture_dissipation_j, fracture, 2e-12);
        near(
            step.energy.constitutive_strain_work_j,
            VOLUME * (opening_work + transverse_storage),
            2e-12,
        );
        near(step.energy.mechanical_work_j, mechanical_work, 2e-12);
        near(step.energy.mean_internal_work_j, mechanical_work, 2e-12);
        near(step.energy.eigenstrain_work_j, eigenstrain_work, 2e-12);
        near(step.motion.kinetic_energy_j(), kinetic, 2e-12);
        near(step.energy.inertia.inertial_work_j, kinetic, 2e-12);
        near(step.energy.inertia.kinetic_change_j, kinetic, 2e-12);
        near(step.energy.parameter_work_j, 0.0, 2e-14);
        near(step.energy.parameter_fracture_dissipation_j, 0.0, 2e-14);
        assert_eq!(step.energy.external_mean_work_j(), 0.0);
        near(mechanical_work, -kinetic, 2e-12);
        near(eigenstrain_work, stored + fracture + kinetic, 2e-12);
        check_midpoint_balances(&step.energy);
        near(step.motion.reference_mass_kg(), DENSITY * VOLUME, 2e-15);
        assert_eq!(step.motion.reference_density_kg_m3(), &[DENSITY]);
        assert_eq!(step.motion.accepted_steps(), 1);
        assert_eq!(step.motion.time_s(), dt);
        assert_eq!(old.solid().solution().element_states[0].damage(), 0.0);
        assert_eq!(old.velocity_mps(), &[[0.0; 3]; 4]);
    });
}

#[test]
fn fixed_temperature_midpoint_vibration_retains_energy_and_its_independent_phase() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = card("midpoint-vibration", 1e-5, false);
        let materials = [source(&card, &geometry, 300.0, 300.0)];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &materials, budget(), cx).unwrap();
        let mut velocity = [[0.0; 3]; 4];
        velocity[1][0] = 0.002;
        let mut motion =
            ThermalBrittleTetMotion::try_from_solid(&solid, &velocity, 0.0, 48, cx).unwrap();
        let mass = DENSITY * VOLUME / 10.0;
        let omega = (YOUNG * VOLUME / mass).sqrt();
        let initial_energy = 0.5 * mass * velocity[1][0].powi(2);
        let dt = 0.01_f64;
        let phase_per_step = 2.0 * (0.5 * omega * dt).atan();
        let supports = axial_supports();
        let mut peak_storage = 0.0_f64;
        for index in 1..=48 {
            let step = motion
                .advance_midpoint(
                    &materials,
                    dt,
                    &supports,
                    &[],
                    settings(),
                    midpoint_controls(),
                    cx,
                )
                .unwrap();
            let phase = f64::from(index) * phase_per_step;
            near(
                step.motion.solid().solution().displacement_m[1][0],
                velocity[1][0] * phase.sin() / omega,
                2e-11,
            );
            near(
                step.motion.velocity_mps()[1][0],
                velocity[1][0] * phase.cos(),
                2e-10,
            );
            near(
                step.motion.kinetic_energy_j() + step.motion.solid().stored_energy_j(),
                initial_energy,
                2e-12,
            );
            near(
                step.motion.kinetic_energy_j(),
                quadrature_kinetic(step.motion.velocity_mps()),
                2e-13,
            );
            assert_eq!(step.energy.fracture_dissipation_j, 0.0);
            assert_eq!(step.energy.eigenstrain_work_j, 0.0);
            near(step.energy.parameter_work_j, 0.0, 2e-14);
            assert_eq!(step.energy.external_mean_work_j(), 0.0);
            assert_eq!(
                step.motion.solid().solution().element_states[0].damage(),
                0.0
            );
            check_midpoint_balances(&step.energy);
            peak_storage = peak_storage.max(step.motion.solid().stored_energy_j());
            motion = step.motion;
        }
        assert!(peak_storage > 0.5 * initial_energy);
        assert_eq!(motion.accepted_steps(), 48);
        assert_eq!(motion.remaining_steps(), 0);
        near(motion.time_s(), 0.48, 2e-15);
        near(motion.reference_mass_kg(), DENSITY * VOLUME, 2e-15);
    });
}

#[test]
fn thermal_midpoint_support_work_uses_original_inventory_and_centered_velocity() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = card("midpoint-expansion", 1e-5, true);
        let reference = [source(&card, &geometry, 300.0, 300.0)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let dt = 0.1;
        let free = 1e-5 * 100.0;
        let targets = affine_targets(free);
        let compatible: [[f64; 3]; 4] =
            core::array::from_fn(|node| core::array::from_fn(|axis| free * NODES[node][axis] / dt));
        let coast = ThermalBrittleTetMotion::try_from_solid(&solid, &compatible, 0.0, 2, cx)
            .unwrap()
            .advance_midpoint(&hot, dt, &targets, &[], settings(), midpoint_controls(), cx)
            .unwrap();
        let driven = ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 2, cx)
            .unwrap()
            .advance_midpoint(&hot, dt, &targets, &[], settings(), midpoint_controls(), cx)
            .unwrap();
        near(
            hot[0].band().reference_density_kg_m3().unwrap(),
            0.5 * DENSITY,
            1e-14,
        );
        let coast_energy = quadrature_kinetic(&compatible);
        for (step, velocity_factor) in [(&coast, 1.0), (&driven, 2.0)] {
            for (node, (u, v)) in step
                .motion
                .solid()
                .solution()
                .displacement_m
                .iter()
                .zip(step.motion.velocity_mps())
                .enumerate()
            {
                for axis in 0..3 {
                    near(u[axis], free * NODES[node][axis], 2e-15);
                    near(v[axis], velocity_factor * compatible[node][axis], 2e-14);
                }
            }
            near(
                step.motion.kinetic_energy_j(),
                velocity_factor.powi(2) * coast_energy,
                2e-14,
            );
            near(step.motion.reference_mass_kg(), DENSITY * VOLUME, 2e-15);
            assert_eq!(step.motion.reference_density_kg_m3(), &[DENSITY]);
            assert_eq!(
                step.motion.solid().solution().element_states[0].damage(),
                0.0
            );
            assert_eq!(step.motion.solid().fracture_dissipation_j(), 0.0);
            for component in 0..6 {
                near(
                    step.motion.solid().solution().element_states[0].strain()[component],
                    0.0,
                    2e-15,
                );
                near(
                    step.motion.solid().solution().element_stress[0][component],
                    0.0,
                    2e-12,
                );
            }
            near(step.energy.mechanical_work_j, 0.0, 2e-14);
            near(step.energy.eigenstrain_work_j, 0.0, 2e-14);
            near(step.energy.parameter_work_j, 0.0, 2e-14);
            assert_eq!(step.energy.dead_load_mean_work_j, 0.0);
            check_midpoint_balances(&step.energy);
        }
        near(coast.energy.inertia.kinetic_change_j, 0.0, 2e-14);
        near(coast.energy.support_mean_work_j, 0.0, 2e-14);
        near(
            driven.energy.inertia.kinetic_change_j,
            4.0 * coast_energy,
            2e-14,
        );
        near(
            driven.energy.inertia.inertial_work_j,
            4.0 * coast_energy,
            2e-14,
        );
        near(driven.energy.support_mean_work_j, 4.0 * coast_energy, 2e-14);
    });
}

#[test]
fn thermal_midpoint_native_refusals_leave_original_history_and_clock_retryable() {
    with_cx(|cx, gate| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = card("midpoint-refusal", 1e-5, false);
        let materials = [source(&card, &geometry, 300.0, 300.0)];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &materials, budget(), cx).unwrap();
        let old =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 1.0, 1, cx).unwrap();
        let snapshot = old.solid().solution().clone();
        let targets: Vec<_> = (0..12)
            .map(|dof| (dof, if dof == 3 { 0.0002 } else { 0.0 }))
            .collect();
        let shifted_reference = [source(&card, &geometry, 350.0, 300.0)];
        assert!(matches!(
            old.advance_midpoint(
                &shifted_reference,
                0.1,
                &targets,
                &[],
                settings(),
                midpoint_controls(),
                cx,
            ),
            Err(BrittleTetDynamicError::Source { .. })
        ));
        // The virgin zero-length mechanical path lacks a classical native
        // mean derivative. Only the explicitly admitted fallback may resolve it.
        let strict = BrittleTetMidpointControls {
            difference_step_strain: None,
            ..midpoint_controls()
        };
        assert!(
            old.advance_midpoint(&materials, 0.1, &targets, &[], settings(), strict, cx,)
                .is_err()
        );
        let oversized = BrittleTetMidpointControls {
            difference_step_strain: Some(0.1),
            ..midpoint_controls()
        };
        assert!(
            old.advance_midpoint(&materials, 0.1, &targets, &[], settings(), oversized, cx,)
                .is_err()
        );
        let stalled = IncrementSettings {
            max_newton_iterations: 1,
            ..settings()
        };
        assert!(matches!(
            old.advance_midpoint(
                &materials,
                0.1,
                &targets,
                &[],
                stalled,
                midpoint_controls(),
                cx,
            ),
            Err(BrittleTetDynamicError::Solid(BrittleTetError::Solve(
                SmallStrainTetError::NewtonStalled { .. }
            )))
        ));
        gate.request();
        assert!(matches!(
            old.advance_midpoint(
                &materials,
                0.1,
                &targets,
                &[],
                settings(),
                midpoint_controls(),
                cx,
            ),
            Err(BrittleTetDynamicError::Solid(BrittleTetError::Cancelled))
        ));
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.velocity_mps(), &[[0.0; 3]; 4]);
        assert_eq!(old.kinetic_energy_j(), 0.0);
        assert_eq!(old.accepted_steps(), 0);
        assert_eq!(old.remaining_steps(), 1);
        assert_eq!(old.time_s(), 1.0);
        let retry_gate = CancelGate::new_clock_free();
        let retry = Cx::new(
            &retry_gate,
            cx.arena(),
            cx.stream_key(),
            cx.budget(),
            cx.mode(),
        );
        let accepted = old
            .advance_midpoint(
                &materials,
                0.1,
                &targets,
                &[],
                settings(),
                midpoint_controls(),
                &retry,
            )
            .unwrap();
        assert!(accepted.numerical_tangent_evaluations > 0);
        near(
            accepted.motion.solid().solution().element_stress[0][0],
            YOUNG * 0.0002,
            2e-12,
        );
        assert_eq!(accepted.motion.accepted_steps(), 1);
        assert_eq!(accepted.motion.remaining_steps(), 0);
        assert_eq!(accepted.motion.time_s(), 1.1);
        check_midpoint_balances(&accepted.energy);
        assert!(matches!(
            accepted.motion.advance_midpoint(
                &materials,
                0.1,
                &targets,
                &[],
                settings(),
                midpoint_controls(),
                &retry,
            ),
            Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: 1,
                maximum_steps: 1
            })
        ));
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.accepted_steps(), 0);
    });
}
