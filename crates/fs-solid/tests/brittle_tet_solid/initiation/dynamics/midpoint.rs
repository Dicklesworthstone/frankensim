//! G1: elastic oscillator phase/energy and source-directed midpoint fracture.

use super::*;
use fs_solid::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use fs_solid::brittle::solid::initiation::dynamics::midpoint::{
    TensileElasticTetMidpointStep, TensileTetMidpointProbe,
};

fn midpoint_controls() -> BrittleTetMidpointControls {
    BrittleTetMidpointControls {
        material_energy_tolerance_j: 1e-10,
        equilibrium_work_tolerance_j: 1e-10,
        kinetic_work_tolerance_j: 1e-10,
        total_energy_tolerance_j: 1e-10,
        difference_step_strain: Some(1e-7),
    }
}

fn midpoint_step<'a>(probe: TensileTetMidpointProbe<'a>) -> TensileElasticTetMidpointStep<'a> {
    match probe {
        TensileTetMidpointProbe::Elastic(step) => step,
        TensileTetMidpointProbe::Onset(_) => {
            panic!("expected a strictly elastic midpoint endpoint")
        }
    }
}

fn supports_y() -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 7)
        .map(|dof| (dof, 0.0))
        .collect()
}

#[test]
fn g1_elastic_midpoint_preserves_energy_and_matches_its_independent_warped_phase() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let density = 2000.0;
        let card = dynamic_source(density);
        let materials = [dynamic_material(&card, &geometry)];
        let virgin = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let length = 0.02;
        let volume = geometry.element_volumes_m3()[0];
        let mass = density * volume / 10.0;
        let stiffness = YOUNG * volume / (length * length);
        let q0 = 0.2 * (STRENGTH / YOUNG) * length;
        let v0 = 0.05;
        let supports = supports_y();
        let solid = elastic(
            virgin
                .probe(
                    &supports,
                    &[(7, stiffness * q0)],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let mut velocity = [[0.0; 3]; 4];
        velocity[2][1] = v0;
        let mut motion = TensileElasticTetMotion::try_from_solid(
            &solid,
            &velocity,
            0.0,
            32,
            initiation_settings(),
            cx,
        )
        .unwrap();
        let h = 1e-6;
        let omega = (stiffness / mass).sqrt();
        let phase_step = 2.0 * (0.5 * h * omega).atan();
        let initial_energy = 0.5 * stiffness * q0 * q0 + 0.5 * mass * v0 * v0;
        let before = motion.solution().clone();
        // The linear native solve can produce a candidate outside the source
        // strain domain. Late source inspection must publish no successor.
        assert!(
            motion
                .probe_midpoint(
                    h,
                    &supports,
                    &[(7, 1e9)],
                    IncrementSettings::default(),
                    midpoint_controls(),
                    None,
                    cx
                )
                .is_err()
        );
        assert_eq!(motion.solution(), &before);
        assert_eq!(motion.accepted_steps(), 0);
        let damped = elastic_step(
            motion
                .probe(
                    h,
                    &supports,
                    &[],
                    IncrementSettings::default(),
                    dynamic_controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        assert!(damped.energy.inertia.numerical_loss_j > 0.0);
        assert!(
            damped.motion.stored_energy_j() + damped.motion.kinetic_energy_j() < initial_energy
        );
        for index in 1..=32 {
            let old_strain = motion.solution().element_strain[0][1];
            let step = midpoint_step(
                motion
                    .probe_midpoint(
                        h,
                        &supports,
                        &[],
                        IncrementSettings::default(),
                        midpoint_controls(),
                        None,
                        cx,
                    )
                    .unwrap(),
            );
            let phase = f64::from(index) * phase_step;
            let q = q0 * phase.cos() + (v0 / omega) * phase.sin();
            let v = -omega * q0 * phase.sin() + v0 * phase.cos();
            close(step.motion.solution().displacement_m[2][1], q, 2e-8);
            close(step.motion.velocity_mps()[2][1], v, 2e-8);
            close(
                step.motion.solution().element_stress[0][1],
                YOUNG * q / length,
                2e-8,
            );
            close(
                step.mean_element_stress_pa[0][1],
                0.5 * YOUNG * (old_strain + q / length),
                2e-8,
            );
            close(
                step.motion.stored_energy_j() + step.motion.kinetic_energy_j(),
                initial_energy,
                2e-8,
            );
            assert_eq!(step.energy.fracture_dissipation_j, 0.0);
            assert_eq!(step.energy.inertia.numerical_loss_j, 0.0);
            assert!(step.energy.total_balance_error_j.abs() < 1e-12);
            assert_eq!(step.motion.reference_mass_kg(), motion.reference_mass_kg());
            motion = step.motion;
        }
        assert_eq!(motion.accepted_steps(), 32);
        assert_eq!(motion.remaining_steps(), 0);
        close(motion.time_s(), 32.0 * h, 1e-14);
        assert!(matches!(
            motion.probe_midpoint(
                h,
                &supports,
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                None,
                cx
            ),
            Err(TensileTetDynamicError::Motion(
                BrittleTetDynamicError::StepBudget { .. }
            ))
        ));
    });
}

#[test]
fn g1_source_midpoint_onset_retains_mean_loads_support_work_and_consumed_motion() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let density = 2000.0;
        let card = dynamic_source(density);
        let materials = [dynamic_material(&card, &geometry)];
        let solid = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let original = TensileElasticTetMotion::try_from_solid(
            &solid,
            &[[0.0; 3]; 4],
            2.0,
            2,
            initiation_settings(),
            cx,
        )
        .unwrap();
        let h = 1e-5;
        let length = 0.02;
        let volume = geometry.element_volumes_m3()[0];
        let area = volume / length;
        let mass = density * volume / 10.0;
        let stiffness = YOUNG * volume / (length * length);
        let q_peak = (STRENGTH / YOUNG) * length;
        let q0 = 0.4 * q_peak;
        let supports = supports_y();
        let first_force = (0.5 * stiffness + 2.0 * mass / (h * h)) * q0;
        let accepted = midpoint_step(
            original
                .probe_midpoint(
                    h,
                    &supports,
                    &[(7, first_force)],
                    IncrementSettings::default(),
                    midpoint_controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let motion = &accepted.motion;
        close(motion.solution().displacement_m[2][1], q0, 2e-9);
        close(motion.velocity_mps()[2][1], 2.0 * q0 / h, 2e-9);
        assert_eq!(motion.accepted_steps(), 1);
        assert_eq!(motion.remaining_steps(), 1);
        let before = motion.solution().clone();
        let velocity = motion.velocity_mps().to_vec();
        let target = 2.0 * q_peak;
        let band = TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap();
        let q_failure = 2.0 * FRACTURE_WORK * length / (STRENGTH * band.band_width_m());
        let peak_force = STRENGTH * area;
        // Independently integrate the bilinear tensile traction across onset.
        let work = 0.5 * stiffness * (q_peak * q_peak - q0 * q0)
            + peak_force / (q_failure - q_peak)
                * (q_failure * (target - q_peak) - 0.5 * (target * target - q_peak * q_peak));
        let delta = target - q0;
        let mean_force = work / delta + 2.0 * mass / (h * h) * (delta - h * velocity[2][1]);
        // Uniform transverse translation gives nonzero native support work.
        let translation = 2e-7;
        let targets: Vec<_> = supports
            .iter()
            .map(|&(dof, _)| (dof, if dof % 3 == 0 { translation } else { 0.0 }))
            .collect();
        let forces = [(7, 0.6 * mean_force), (0, 2.0), (7, 0.4 * mean_force)];
        let prediction = match motion
            .probe_midpoint(
                h,
                &targets,
                &forces,
                IncrementSettings::default(),
                midpoint_controls(),
                None,
                cx,
            )
            .unwrap()
        {
            TensileTetMidpointProbe::Onset(onset) => onset,
            TensileTetMidpointProbe::Elastic(_) => panic!("expected a tensile midpoint crossing"),
        };
        assert_eq!(prediction.initial_solution(), &before);
        assert_eq!(prediction.requested_dt_s(), h);
        assert_eq!(
            prediction.requested_interval_s(),
            [motion.time_s(), motion.time_s() + h]
        );
        assert_eq!(prediction.prescribed_end_m(), targets.as_slice());
        assert_eq!(prediction.mean_nodal_forces_n(), &forces);
        let normal = prediction.events()[0].as_ref().unwrap().native().normal();
        close(normal[1].abs(), 1.0, 2e-12);
        assert!(normal[0].abs() < 1e-12 && normal[2].abs() < 1e-12);
        let predictor_q = (mean_force + 2.0 * mass / (h * h) * (q0 + h * velocity[2][1])
            - 0.5 * stiffness * q0)
            / (2.0 * mass / (h * h) + 0.5 * stiffness);
        close(
            prediction.elastic_predictor().displacement_m[2][1],
            predictor_q,
            2e-9,
        );
        close(
            prediction.elastic_predictor().element_stress[0][1],
            YOUNG * predictor_q / length,
            2e-9,
        );
        assert!(predictor_q > q_peak);
        assert!(0.5 * (q0 + predictor_q) < q_peak);
        let bands = prediction.resolve_bands(cx).unwrap();
        let handoff = prediction.handoff(&bands, 1e-9, cx).unwrap();
        assert_eq!(
            handoff.motion.solid().solution().displacement_m,
            before.displacement_m
        );
        assert_eq!(handoff.motion.velocity_mps(), velocity.as_slice());
        assert_eq!(
            handoff.motion.reference_mass_kg(),
            motion.reference_mass_kg()
        );
        assert_eq!(handoff.motion.kinetic_energy_j(), motion.kinetic_energy_j());
        assert_eq!(handoff.motion.time_s(), motion.time_s());
        assert_eq!(handoff.motion.accepted_steps(), 1);
        assert_eq!(handoff.motion.remaining_steps(), 1);
        assert_eq!(handoff.motion.solid().fracture_dissipation_j(), 0.0);
        with_cx(true, |cancelled| {
            assert!(
                prediction
                    .replay(&bands, 1e-9, midpoint_controls(), cancelled)
                    .is_err()
            );
        });
        let replay = prediction
            .replay(&bands, 1e-9, midpoint_controls(), cx)
            .unwrap();
        close(
            replay.motion.solid().solution().displacement_m[2][1],
            target,
            2e-8,
        );
        close(
            replay.motion.velocity_mps()[2][1],
            2.0 * delta / h - velocity[2][1],
            2e-8,
        );
        close(replay.energy.material_work_j, work, 2e-8);
        close(replay.energy.mean_internal_work_j, work, 2e-8);
        close(
            replay.energy.dead_load_mean_work_j,
            mean_force * delta + 2.0 * translation,
            2e-8,
        );
        let translation_energy = 0.5 * (density * volume) * (2.0 * translation / h).powi(2);
        close(
            replay.energy.support_mean_work_j,
            translation_energy - 2.0 * translation,
            2e-8,
        );
        let stress = STRENGTH * (q_failure - target) / (q_failure - q_peak);
        let stored_after = 0.5 * volume * stress * (target / length);
        close(
            replay.energy.fracture_dissipation_j,
            work - (stored_after - motion.stored_energy_j()),
            2e-8,
        );
        assert!(replay.energy.fracture_dissipation_j > 0.0);
        assert_eq!(replay.energy.inertia.numerical_loss_j, 0.0);
        assert!(replay.energy.total_balance_error_j.abs() < 1e-9);
        assert_eq!(replay.motion.accepted_steps(), 2);
        assert_eq!(replay.motion.remaining_steps(), 0);
        assert_eq!(replay.motion.time_s(), motion.time_s() + h);
        assert_eq!(motion.solution(), &before);
        assert_eq!(motion.velocity_mps(), velocity.as_slice());
        let repeated = prediction
            .replay(&bands, 1e-9, midpoint_controls(), cx)
            .unwrap();
        assert_eq!(
            repeated.motion.solid().solution(),
            replay.motion.solid().solution()
        );
        assert_eq!(repeated.motion.velocity_mps(), replay.motion.velocity_mps());
    });
}
