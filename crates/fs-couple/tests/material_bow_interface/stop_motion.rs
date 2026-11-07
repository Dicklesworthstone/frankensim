//! End-to-end support motion on the original sourced string/twist/body fixture.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::{
    MultiBowPerformance, motion::{ObstacleMotionLimits, StopMotionKnot, StopMotionProgram},
};

fn motion_program() -> StopMotionProgram {
    StopMotionProgram {
        knots: [(12, -3e-6), (36, -3e-6), (72, 0.0)].into_iter()
            .map(|(sample, offset)| StopMotionKnot { sample, offsets_m: vec![offset] }).collect(),
        limits: ObstacleMotionLimits { maximum_abs_offset_m: 0.001, maximum_speed_m_s: 1.0 },
        max_knots: 16, max_coordinate_values: 32,
    }
}
fn offset(endpoint: u64) -> f64 {
    if endpoint < 12 { -3e-6*(endpoint as f64/12.0) }
    else if endpoint <= 36 { -3e-6 }
    else if endpoint < 72 { -3e-6+3e-6*((endpoint-36) as f64/36.0) }
    else { 0.0 }
}
fn performance(cx: &Cx<'_>, gate: &CancelGate) -> MultiBowPerformance {
    stopped(cx, gate).into_scheduled(schedule(), vec!["first".into(), "second".into()],
        MultiBowPerformanceConfig { samples: 96, max_block: 64, max_work: 10000, max_events: 64 }).unwrap()
}

#[test]
fn motion_bows_and_load_release_use_one_pressure_history_independent_of_blocks() {
    with_cx(|cx, gate| {
        let mut direct = stopped(cx, gate).with_stop_motion(motion_program().limits, gate).unwrap();
        let load = forced(&direct, gate); let zero = vec![0.0; load.len()];
        let mut expected = Vec::new(); let mut moved_work = false;
        for sample in 0..96 {
            let frame = direct.step_with_stop_offsets(if sample < 64 { &load } else { &zero },
                &drivers(), &[offset(sample+1)]).unwrap().clone();
            moved_work |= frame.obstacle_work_j != 0.0;
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            expected.push(frame.observer_pressure_pa.to_bits());
        }
        assert!(moved_work);
        for partition in [1, 37, 64] {
            let mut performance = performance(cx, gate)
                .with_string_forces(program(ForceInitialization::RetainState), gate).unwrap()
                .with_stop_motion(motion_program(), gate).unwrap();
            assert_eq!(performance.stop_motion_program(), Some(&motion_program()));
            let mut output = vec![0.0; 96];
            for block in output.chunks_mut(partition) { performance.block(block).unwrap(); }
            assert_eq!(output.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), expected);
            assert_eq!(performance.bank().system().obstacle_offsets_m().unwrap(), &[0.0]);
            assert_eq!(performance.bank().source().stops().len(), 1);
            assert_eq!(performance.bank().bow_interface(0), Some(&interface(0.4, 10.0)));
            assert!(performance.bank().system().total_energy_j().unwrap() > 0.0);
        }
    });
}

#[test]
fn geometry_motion_changes_actual_pressure_without_retuning_material_modes() {
    with_cx(|cx, gate| {
        let mut fixed = performance(cx, gate);
        let mut moving = performance(cx, gate).with_stop_motion(motion_program(), gate).unwrap();
        for (a,b) in fixed.bank().system().components().iter().zip(moving.bank().system().components()) {
            assert_eq!(a.modes(), b.modes());
        }
        let mut a = vec![0.0; 96]; let mut b = a.clone();
        for block in a.chunks_mut(64) { fixed.block(block).unwrap(); }
        for block in b.chunks_mut(37) { moving.block(block).unwrap(); }
        assert!(a.iter().zip(&b).any(|(x,y)| (x-y).abs() > 1e-15));
        assert_ne!(fixed.bank().system().components()[0].states(), moving.bank().system().components()[0].states());
    });
}

#[test]
fn contact_preload_is_identical_on_either_side_of_motion_binding() {
    with_cx(|cx, gate| {
        let mut before = performance(cx, gate)
            .with_string_forces(program(ForceInitialization::StaticPreload), gate).unwrap()
            .with_stop_motion(motion_program(), gate).unwrap();
        let mut after = performance(cx, gate).with_stop_motion(motion_program(), gate).unwrap()
            .with_string_forces(program(ForceInitialization::StaticPreload), gate).unwrap();
        assert_eq!(before.preload_energy_j(), after.preload_energy_j());
        assert_eq!(before.samples_rendered(), 0);
        assert_eq!(after.bank().system().obstacle_offsets_m().unwrap(), &[0.0]);
        for (a,b) in before.bank().system().components().iter().zip(after.bank().system().components()) {
            assert_eq!(a.states(), b.states());
        }
        let mut a = vec![0.0; 96]; let mut b = a.clone();
        for block in a.chunks_mut(64) { before.block(block).unwrap(); }
        for block in b.chunks_mut(37) { after.block(block).unwrap(); }
        assert_eq!(a.iter().map(|x|x.to_bits()).collect::<Vec<_>>(), b.iter().map(|x|x.to_bits()).collect::<Vec<_>>());
    });
}

#[test]
fn malformed_motion_and_invalid_callbacks_never_advance_geometry() {
    with_cx(|cx, gate| {
        for which in 0..8 {
            let mut motion = motion_program();
            match which {
                0 => motion.knots[0].sample = 0,
                1 => motion.knots[1].sample = 12,
                2 => motion.knots[0].offsets_m.clear(),
                3 => motion.knots[0].offsets_m[0] = f64::NAN,
                4 => motion.knots[0].offsets_m[0] = -0.1,
                5 => motion.limits.maximum_speed_m_s = 1e-8,
                6 => motion.knots[2].sample = 97,
                _ => motion.max_coordinate_values = 2,
            }
            assert!(performance(cx, gate).with_stop_motion(motion, gate).is_err());
        }
        let mut performance = performance(cx, gate).with_stop_motion(motion_program(), gate).unwrap();
        assert!(performance.block(&mut [0.0; 65]).is_err());
        assert_eq!(performance.samples_rendered(), 0);
        assert_eq!(performance.bank().system().obstacle_offsets_m().unwrap(), &[0.0]);
        performance.block(&mut [0.0; 12]).unwrap();
        assert_eq!(performance.bank().system().obstacle_offsets_m().unwrap(), &[-3e-6]);
        assert!(performance.with_stop_motion(motion_program(), gate).is_err());
    });
}


/// Keep the real material/body/source path while condensing its omitted modes.
fn residual_stopped(cx: &Cx<'_>, gate: &CancelGate) -> MultiBowStringBank {
    let binding = interface(0.4, 10.0);
    build(cx, gate, binding.clone())
        .with_residual_flexibility(gate).unwrap()
        .with_additional_bows(vec![MaterialBowContact {
            location: StringBowLocation { string: 0, station_fraction: 0.31 },
            law: binding.card().friction_law().clone(),
            friction: options().friction,
        }], MultiFrictionConfig {
            max_ports: 4, max_sweeps: 128, max_setup_terms: 100_000,
        }, gate).unwrap()
        .with_bow_interface(1, binding).unwrap()
        .with_string_stops(vec![stop()], limits(), gate).unwrap()
}

#[test]
fn condensed_material_stops_preserve_motion_work_and_scheduled_pressure() {
    with_cx(|cx, gate| {
        let mut direct = residual_stopped(cx, gate)
            .with_stop_motion(motion_program().limits, gate).unwrap();
        let compliance = direct.system().obstacle_residual_compliance().unwrap()[0];
        assert!(compliance > 0.0);
        let column = direct.system().obstacle_normal_column(0).unwrap().to_vec();
        let source = direct.source().stops()[0].law.clone();
        let load = forced(&direct, gate);
        let zero = vec![0.0; load.len()];
        let (mut previous_reaction, mut previous_offset) = (0.0, 0.0);
        let mut expected = Vec::new();
        let mut supplied_work = false;
        let mut residual_storage = false;
        for sample in 0..96 {
            let endpoint = offset(sample + 1);
            let frame = direct.step_with_stop_offsets(
                if sample < 64 { &load } else { &zero }, &drivers(), &[endpoint],
            ).unwrap().clone();
            let normal = direct.system().last_obstacle_frames().unwrap()[0];
            let retained_closure: f64 = direct.system().components().iter()
                .flat_map(|model| model.states())
                .zip(&column)
                .map(|(state, shape)| shape * state.displacement_m_sqrt_kg)
                .sum();
            let closure = retained_closure - compliance * normal.normal_force_n - endpoint;
            let penetration = (closure - source.gaps()[0]).max(0.0);
            let stored = 0.5 * source.stiffness() * penetration * penetration;
            assert!((normal.penetration_m - penetration).abs() < 1e-12);
            assert!((normal.stored_energy_j - stored).abs() < 1e-13);
            let reaction_change = normal.normal_force_n - previous_reaction;
            assert!((frame.obstacle_residual_energy_j
                - 0.5 * compliance * normal.normal_force_n.powi(2)).abs() < 1e-13);
            assert!((frame.obstacle_residual_loss_j
                - 0.5 * compliance * reaction_change.powi(2)).abs() < 1e-13);
            assert!((frame.obstacle_work_j
                + normal.normal_force_n * (endpoint - previous_offset)).abs() < 1e-13);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            supplied_work |= frame.obstacle_work_j != 0.0;
            residual_storage |= frame.obstacle_residual_energy_j > 0.0;
            expected.push(frame.observer_pressure_pa.to_bits());
            (previous_reaction, previous_offset) = (normal.normal_force_n, endpoint);
        }
        assert!(supplied_work && residual_storage);
        assert!(direct.system().components()[1].states().iter()
            .any(|state| state.velocity_m_sqrt_kg_per_s != 0.0));

        for partition in [1, 37, 64] {
            let mut performance = residual_stopped(cx, gate)
                .into_scheduled(schedule(), vec!["first".into(), "second".into()],
                    MultiBowPerformanceConfig {
                        samples: 96, max_block: 64, max_work: 10_000, max_events: 64,
                    }).unwrap()
                .with_string_forces(program(ForceInitialization::RetainState), gate).unwrap()
                .with_stop_motion(motion_program(), gate).unwrap();
            let mut output = vec![0.0; 96];
            for block in output.chunks_mut(partition) { performance.block(block).unwrap(); }
            assert_eq!(output.iter().map(|value| value.to_bits()).collect::<Vec<_>>(), expected);
            assert_eq!(performance.bank().system().obstacle_offsets_m().unwrap(), &[0.0]);
            assert_eq!(performance.bank().system().obstacle_residual_compliance().unwrap(), &[compliance]);
            for port in 0..2 {
                assert_eq!(performance.bank().bow_interface(port), Some(&interface(0.4, 10.0)));
            }
        }
    });
}
