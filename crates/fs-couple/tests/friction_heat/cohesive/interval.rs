//! G1/G3/G4: accepted deformable stop paths heat original finite receivers.
use super::*;
use fs_couple::friction_heat::cohesive::interval::*;
use fs_solid::cohesive::dynamics::interval::IntervalConfig;

fn interval_controls() -> IntervalConfig {
    IntervalConfig {
        minimum_step_s: 1e-12,
        reversal_time_tolerance_s: 1e-11,
        maximum_trials: 160,
        maximum_substeps: 40,
        maximum_total_residual_evaluations: 200_000,
        maximum_total_newton_iterations: 6000,
    }
}

fn coulomb_controls() -> CoulombConfig {
    CoulombConfig {
        traction_scale_pa_per_m: 1000.0,
        slip_tolerance_m: 1e-11,
        velocity_tolerance_m_s: 1e-8,
    }
}

#[test]
fn unaligned_deformable_stopping_heat_survives_static_recapture_and_breakaway() {
    let source = material(false);
    let elastic = Scalar(0.0);
    let mass = Scalar(0.1);
    let joint = system(&elastic, &source, &[7], true);
    let incoming = 0.037;
    let motion = compressed(&joint, 0.005)
        .into_dynamics(&mass, vec![incoming], 0.0, dynamics())
        .unwrap()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::Sliding)],
            coulomb_controls(),
            &mut || false,
        )
        .unwrap();
    let thermal = network(
        1,
        [1.0, 3.0],
        [2.0, 9.0],
        [0.5, 2.0],
        [0.0; 2],
        [40.0, 10.0],
    );
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [2.0, 9.0]).unwrap();
    let contacts = [CohesiveHeatContact {
        patch_id: 7,
        route: &route,
        temperature_guess_k: 320.0,
    }];
    let duration = 0.04;
    let offsets = [[-0.001, 0.0, 0.0]];
    let plain = motion
        .advance_interval_controlled(&[0.0], &offsets, duration, interval_controls(), &mut || {
            false
        })
        .unwrap();
    let heated = propose_deformable_cohesive_interval_heat(
        &motion,
        &[0.0],
        &offsets,
        duration,
        interval_controls(),
        &contacts,
        &thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids, |_| [2.0, 9.0])),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(
        heated.mechanical.motion.state().points(),
        plain.motion.state().points()
    );
    assert_eq!(
        heated.mechanical.motion.state().displacements_m(),
        plain.motion.state().displacements_m()
    );
    assert_eq!(
        heated.mechanical.motion.velocity_m_s(),
        plain.motion.velocity_m_s()
    );
    assert_eq!(heated.mechanical.patch_dissipation, plain.patch_dissipation);
    assert!(heated.mechanical.attempted_substeps > heated.mechanical.accepted_substeps);
    assert!(
        heated.mechanical.motion.velocity_m_s()[0].abs()
            < coulomb_controls().velocity_tolerance_m_s
    );
    let heat = 0.5 * mass.0 * incoming * incoming;
    close(heated.heat.dissipated_j, heat, 1e-12);
    close(
        heated.mechanical.motion.state().displacements_m()[0],
        heat / 0.6,
        1e-11,
    );
    assert_eq!(heated.fracture_work_j, 0.0);
    assert_eq!(heated.mechanical.viscous_dissipated_j, 0.0);
    assert_eq!(
        heated
            .mechanical
            .motion
            .damaged_area_coulomb_states()
            .unwrap()[0]
            .occupancy,
        CoulombRegime::Sticking
    );
    let mut matrix = factor(&thermal, duration, [4.0, 6.0], 0.0);
    let (energy, face) = endpoint(&thermal, &mut matrix, duration, heat / duration);
    for (i, cell) in thermal.cells().values().enumerate() {
        close(
            heated.thermal.cells()[&cell.id].specific_h_j_kg,
            cell.specific_h_j_kg + energy[i] / cell.mass_kg,
            2e-10,
        );
    }
    close(
        heated.heat.interfaces[0].interface_temperature_k,
        face,
        2e-10,
    );
    assert!(heated.heat.interfaces[0].body_heat_j[0] < 0.0);
    close(heated.combined_energy_residual_j, 0.0, 2e-10);
    close(heated.summed_energy_residual_j, 0.0, 2e-10);
    assert_eq!(heated.mechanical.motion.time_s(), heated.thermal.time_s());
    let held = propose_deformable_cohesive_interval_heat(
        &heated.mechanical.motion,
        &[0.8],
        &offsets,
        0.013,
        interval_controls(),
        &contacts,
        &heated.thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids, |_| [2.0, 9.0])),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(held.heat.dissipated_j, 0.0); // .6 kinetic < applied .8 < 1.0 static.
    assert_eq!(held.mechanical.motion.time_s(), duration + 0.013);
    assert_eq!(held.mechanical.motion.time_s(), held.thermal.time_s());
    let sliding = propose_deformable_cohesive_interval_heat(
        &held.mechanical.motion,
        &[1.2],
        &offsets,
        0.01,
        interval_controls(),
        &contacts,
        &held.thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(inputs(ids, |_| [2.0, 9.0])),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(sliding.heat.dissipated_j > 0.0);
    assert_eq!(
        sliding
            .mechanical
            .motion
            .damaged_area_coulomb_states()
            .unwrap()[0]
            .occupancy,
        CoulombRegime::Sliding
    );
    close(
        sliding.heat.dissipated_j,
        sliding.mechanical.friction_dissipated_j,
        1e-15,
    );
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(motion.velocity_m_s(), &[incoming]);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn prescribed_multisite_intervals_keep_surface_creation_and_shared_cell_heat_separate() {
    let source = material(true);
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7, 9], false);
    let start = [-0.001, 0.0005, 0.0];
    let motion = joint
        .bonded_state()
        .solve_over_time_controlled(&[0.0], &[start; 2], 1.0, dynamics().solve, &mut || false)
        .unwrap()
        .state
        .into_dynamics(&mass, vec![0.0], 0.0, dynamics())
        .unwrap()
        .with_damaged_area_friction_controlled(&[Some(1e-4); 2], &mut || false)
        .unwrap();
    let thermal = three_cells();
    let routes = [
        FrictionHeatRoute::from_link(&thermal, 20, [1.0, 2.0]).unwrap(),
        FrictionHeatRoute::from_link(&thermal, 10, [2.0, 3.0]).unwrap(),
    ];
    let contacts = [
        CohesiveHeatContact {
            patch_id: 7,
            route: &routes[0],
            temperature_guess_k: 315.0,
        },
        CohesiveHeatContact {
            patch_id: 9,
            route: &routes[1],
            temperature_guess_k: 325.0,
        },
    ];
    let duration = 0.04;
    let ends: [f64; 2] = [0.003, 0.002];
    let heated = propose_deformable_cohesive_interval_heat(
        &motion,
        &[0.0],
        &ends.map(|q| [-0.001, q, 0.0]),
        duration,
        interval_controls(),
        &contacts,
        &thermal,
        controls(3, 2),
        &[0.0; 3],
        |_, ids, _| {
            Ok::<_, &'static str>(inputs(
                ids,
                |id| if id == 10 { [2.0, 3.0] } else { [1.0, 2.0] },
            ))
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(heated.mechanical.accepted_substeps, 2);
    assert_eq!(heated.patch_ids, [9, 7]);
    let mut expected_heat = 0.0;
    let mut expected_fracture = 0.0;
    for (index, end) in ends.into_iter().enumerate() {
        // Integrate the actual surviving fraction: 1 to .001, then
        // 1-d(q) = -.333333... + .001333333.../q to the endpoint.
        let surviving = 0.0005 - (end - 0.001) / 3.0 + (0.004 / 3.0) * (end / 0.001).ln();
        let friction = 0.6 * (end - start[1] - surviving);
        let viscosity = 4.0 * (end - start[1]) / duration * surviving;
        let damage = (4.0 / 3.0) * (1.0 - 0.001 / end);
        let work = 0.5 * 1000.0 * (0.001_f64.powi(2) - start[1].powi(2))
            + 0.5 * (1.0 + (0.004 - end) / 0.003) * (end - 0.001);
        let stored = 0.5 * 1000.0 * ((1.0 - damage) * end * end - start[1] * start[1]);
        let loss = heated.mechanical.patch_dissipation[index];
        close(loss.friction_j, friction, 2e-12);
        close(loss.viscous_j, viscosity, 2e-12);
        close(loss.fracture_j, work - stored, 2e-12);
        close(
            heated.heat.interfaces[1 - index].dissipated_j,
            friction + viscosity,
            2e-12,
        );
        expected_heat += friction + viscosity;
        expected_fracture += work - stored;
    }
    close(heated.fracture_work_j, expected_fracture, 3e-12);
    close(heated.heat.dissipated_j, expected_heat, 3e-12);
    close(
        heated.thermal.energy_j() - thermal.energy_j(),
        expected_heat,
        2e-10,
    );
    assert!(expected_fracture > 0.0);
    close(heated.combined_energy_residual_j, 0.0, 2e-10);
    close(heated.summed_energy_residual_j, 0.0, 2e-10);
    assert_eq!(heated.mechanical.maximum_normal_jump_m, [-0.001; 2]);
}

#[test]
fn late_thermal_refusal_after_deformable_stopping_preserves_original_motion_and_heat() {
    let source = material(false);
    let elastic = Scalar(0.0);
    let mass = Scalar(0.1);
    let joint = system(&elastic, &source, &[7], true);
    let motion = compressed(&joint, 0.005)
        .into_dynamics(&mass, vec![0.037], 0.0, dynamics())
        .unwrap()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::Sliding)],
            coulomb_controls(),
            &mut || false,
        )
        .unwrap();
    let thermal = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [1.0; 2]).unwrap();
    let contacts = [CohesiveHeatContact {
        patch_id: 7,
        route: &route,
        temperature_guess_k: 300.0,
    }];
    let points = motion.state().points().to_vec();
    let friction = motion.damaged_area_coulomb_states().unwrap().to_vec();
    let cells = thermal.cells().clone();
    for cancel in [false, true] {
        let gate = CancelGate::new_clock_free();
        let calls = Cell::new(0);
        let refused = propose_deformable_cohesive_interval_heat(
            &motion,
            &[0.0],
            &[[-0.001, 0.0, 0.0]],
            0.04,
            interval_controls(),
            &contacts,
            &thermal,
            controls(2, 1),
            &[0.0; 2],
            |_, ids, _| {
                calls.set(calls.get() + 1);
                if cancel {
                    gate.request();
                    Ok(inputs(ids, |_| [1.0; 2]))
                } else {
                    Err("late thermal source")
                }
            },
            &gate,
        );
        assert!(refused.is_err());
        assert!(calls.get() > 0);
        assert_eq!(motion.state().points(), points);
        assert_eq!(motion.damaged_area_coulomb_states().unwrap(), friction);
        assert_eq!(motion.state().displacements_m(), &[0.0]);
        assert_eq!(motion.velocity_m_s(), &[0.037]);
        assert_eq!(motion.steps(), 0);
        assert_eq!(motion.time_s(), thermal.time_s());
        assert_eq!(thermal.cells(), &cells);
    }
}

#[test]
fn deformable_opening_inside_an_interval_refuses_despite_final_reclosure() {
    let source = material(false);
    let elastic = Scalar(0.0);
    let mass = Scalar(10.0);
    let patch = CohesivePatch::from_interface(
        7,
        1.0,
        &source,
        [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
        vec![JumpTerm {
            dof: 0,
            local: [1.0, 0.0, 0.0],
        }],
    )
    .unwrap();
    let joint = CohesiveSystem::new(
        &elastic,
        vec![patch],
        Budget {
            max_dofs: 1,
            max_patches: 1,
            max_jump_terms: 1,
            max_krylov_basis_entries: 128,
        },
    )
    .unwrap();
    let offsets = [[-0.001, 0.0, 0.0]];
    let motion = joint
        .bonded_state()
        .solve_controlled(&[-2.0], &offsets, dynamics().solve, &mut || false)
        .unwrap()
        .state
        .into_dynamics(&mass, vec![0.04], 0.0, dynamics())
        .unwrap();
    let duration = 0.4;
    let plain = motion
        .advance_interval_controlled(
            &[-2.0],
            &offsets,
            duration,
            interval_controls(),
            &mut || false,
        )
        .unwrap();
    assert!(plain.maximum_normal_jump_m[0] > 0.0);
    assert!(plain.motion.state().points()[0].jump_m()[0] <= 0.0);
    let thermal = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [1.0; 2]).unwrap();
    let calls = Cell::new(0);
    let refused = propose_deformable_cohesive_interval_heat(
        &motion,
        &[-2.0],
        &offsets,
        duration,
        interval_controls(),
        &[CohesiveHeatContact {
            patch_id: 7,
            route: &route,
            temperature_guess_k: 300.0,
        }],
        &thermal,
        controls(2, 1),
        &[0.0; 2],
        |_, ids, _| {
            calls.set(calls.get() + 1);
            Ok::<_, &'static str>(inputs(ids, |_| [1.0; 2]))
        },
        &CancelGate::new_clock_free(),
    );
    assert!(matches!(refused, Err(CohesiveHeatError::Invalid(_))));
    assert_eq!(calls.get(), 0);
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(motion.velocity_m_s(), &[0.04]);
}
