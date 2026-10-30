//! Normal load, moving material crossings and support work share one solve.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    performance::bow_stations::{BowStationKnot, BowStationMotionLimits, BowStationProgram},
    ribbon::UniformBowRibbon,
};

#[test]
fn normal_support_acceleration_matches_the_continuum_and_width_partition() {
    // G1: independent forced sine oscillator; rigid normal support travel must
    // accelerate the source mass even while the string contact is separated.
    with_cx(|cx, gate| {
        let samples = 96;
        let rate = 48000.0;
        let mut previous_energy: Option<f64> = None;
        for cells in [1, 3, 5] {
            let fiber = fibers(cells, 2, 1300.0, 3e9, 0.0, 40.0, rate as u32, gate);
            let omega = fiber.cell_model().modes()[2].angular_frequency_rad_s;
            let gamma = 2.0 * (2.0 * fiber.mass_kg() / cells as f64).sqrt()
                / core::f64::consts::PI;
            let bank = base(cx, gate, 7800.0).with_tensioned_fiber_ribbon(
                fiber, 0.4, ribbon_limits(), gate).unwrap();
            let first = bank.source().fiber_ribbon().unwrap().component_start();
            let config = MultiBowPerformanceConfig { samples, ..performance_config() };
            let mut run = bank.into_scheduled_ribbon(fiber_load_schedule(0.0), "bow".into(), config).unwrap()
                .with_fiber_support_program(support_program(&[(samples, 0.0)]), gate).unwrap()
                .with_fiber_normal_support_program(support_program(&[(samples, samples as f64 / rate)]), gate).unwrap();
            let mut work = 0.0;
            for step in 1..=samples {
                let mut output = [0.0];
                run.block(&mut output).unwrap();
                assert_eq!(output, [0.0]);
                let t = step as f64 / rate;
                let state = run.bank().system().components()[first].states()[2];
                let q = -gamma * (1.0 - (omega * t).cos()) / omega.powi(2);
                let v = -gamma * (omega * t).sin() / omega;
                assert!((state.displacement_m_sqrt_kg - q).abs() < 1e-17);
                assert!((state.velocity_m_sqrt_kg_per_s - v).abs() < 1e-13);
                assert_eq!(run.bank().system().components()[first].states()[0], ModalAcousticState::default());
                assert_eq!(run.bank().system().components()[first].states()[3], ModalAcousticState::default());
                let frame = run.fiber_support_frame().unwrap();
                let normal = frame.normal.unwrap();
                assert!((normal.translation_m - 0.5 * t * t).abs() < 1e-19);
                assert_eq!(normal.contact_translation_work_j, 0.0);
                work += frame.support_work_j;
                assert!((frame.laboratory_mechanical_energy_j - work).abs() < 1e-14);
            }
            let energy = run.fiber_support_frame().unwrap().laboratory_mechanical_energy_j;
            assert!(energy > 0.0);
            if let Some(previous) = previous_energy {
                assert!((energy - previous).abs() < 1e-18);
            }
            previous_energy = Some(energy);
        }
    });
}

fn pressing_program() -> FiberSupportProgram {
    support_program(&[(64, -0.02), (128, 0.0), (192, 0.02), (256, 0.0)])
}

#[test]
fn normal_support_press_and_release_solve_load_and_close_laboratory_work() {
    with_cx(|cx, gate| {
        let make = |tension| normal_fiber_bank(cx, gate, tension, 5e-6)
            .into_scheduled_ribbon(fiber_load_schedule(0.0), "bow".into(), performance_config()).unwrap()
            .with_fiber_support_program(support_program(&[(64, 0.01), (384, 0.01)]), gate).unwrap()
            .with_fiber_normal_support_program(pressing_program(), gate).unwrap();
        let mut actual = make(40.0);
        let mut chunks = make(40.0);
        let mut taut = make(160.0);
        let (mut loaded, mut separated, mut balance) = (false, false, 0.0);
        let mut samples = Vec::new();
        for _ in 0..513 {
            let mut output = [0.0];
            actual.block(&mut output).unwrap();
            samples.push(output[0]);
            let frame = actual.bank().system().last_frame().unwrap();
            let support = actual.fiber_support_frame().unwrap();
            let normal = support.normal.unwrap();
            let reactions = actual.bank().system().last_obstacle_frames().unwrap();
            for (port, reaction) in frame.ports.iter().zip(reactions) {
                assert_eq!(port.normal_force_n, reaction.normal_force_n);
                loaded |= reaction.normal_force_n > 0.0;
                separated |= reaction.normal_force_n == 0.0;
                if reaction.normal_force_n == 0.0 { assert_eq!(port.traction_n, 0.0); }
            }
            balance += frame.external_work_j - support.inertial_coordinate_work_j
                + support.support_work_j + frame.obstacle_work_j
                - support.normal_attachment_work_j - normal.contact_translation_work_j
                - frame.network_dissipation_j - frame.friction_dissipation_j
                - frame.obstacle_dissipation_j - frame.obstacle_residual_loss_j;
            assert!((support.laboratory_mechanical_energy_j - balance).abs() < 2e-9);
            for offset in actual.bank().system().obstacle_offsets_m().unwrap() {
                assert_eq!(*offset, normal.translation_m);
            }
        }
        let mut other = vec![0.0; 513];
        let mut changed = other.clone();
        for block in other.chunks_mut(37) { chunks.block(block).unwrap(); }
        for block in changed.chunks_mut(128) { taut.block(block).unwrap(); }
        assert!(loaded && separated);
        assert_eq!(samples, other);
        assert_eq!(states(actual.bank()), states(chunks.bank()));
        assert_eq!(actual.fiber_support_frame(), chunks.fiber_support_frame());
        assert!(samples.iter().any(|v| v.abs() > 1e-12));
        assert!(samples.iter().zip(changed).any(|(a, b)| (a-b).abs() > 1e-12));
    });
}

#[test]
fn refused_normal_translation_preserves_support_contact_and_material_history() {
    with_cx(|cx, gate| {
        let mut under_budget = pressing_program();
        under_budget.max_update_terms = 1;
        assert!(normal_fiber_bank(cx, gate, 40.0, 5e-6)
            .with_fiber_support_motion(0.0, support_limits(), gate).unwrap()
            .with_fiber_normal_support_program(under_budget, 513, gate).is_err());
        let make = || normal_fiber_bank(cx, gate, 40.0, 5e-6)
            .with_fiber_support_motion(0.0, support_limits(), gate).unwrap()
            .with_fiber_normal_support_program(pressing_program(), 513, gate).unwrap();
        let mut actual = make();
        let mut expected = make();
        let forces = vec![0.0; actual.system().mode_count()];
        let drives = vec![drive(0.0, 0.0); 3];
        assert_eq!(actual.step(&forces, &drives).unwrap(), expected.step(&forces, &drives).unwrap());
        let before = states(&actual);
        let support = actual.fiber_support_frame().copied();
        let normal = actual.fiber_normal_support_kinematics();
        let offsets = actual.system().obstacle_offsets_m().unwrap().to_vec();
        // Two prescribed support translations cannot own one physical boundary.
        assert!(actual.step_with_stop_offsets(&forces, &drives, &[1e-7; 3]).is_err());
        assert_eq!(states(&actual), before);
        assert_eq!(actual.fiber_support_frame().copied(), support);
        assert_eq!(actual.fiber_normal_support_kinematics(), normal);
        assert_eq!(actual.system().obstacle_offsets_m().unwrap(), offsets);
        assert_eq!(actual.step(&forces, &drives).unwrap(), expected.step(&forces, &drives).unwrap());
        assert_eq!(states(&actual), states(&expected));
        assert_eq!(actual.fiber_support_frame(), expected.fiber_support_frame());
    });
}

#[test]
fn a_noncanonical_normal_row_cannot_acquire_invented_material_geometry() {
    with_cx(|cx, gate| {
        let altered = normal_fiber_bank_with_receiver_scale(cx, gate, 40.0, -1e-6, 1.01);
        // An authored row is usable as a fixed supplied receiver, but its
        // scale does not establish a station on the retained fiber continuum.
        assert!(
            altered
                .with_fiber_support_motion(0.0, support_limits(), gate)
                .is_err()
        );
    });
}

#[test]
fn a_stationary_support_preserves_the_original_solved_normal_contact_bitwise() {
    with_cx(|cx, gate| {
        let mut original = normal_fiber_bank(cx, gate, 40.0, -1e-6);
        let mut moving = normal_fiber_bank(cx, gate, 40.0, -1e-6)
            .with_fiber_support_motion(0.0, support_limits(), gate)
            .unwrap();
        let zero = vec![0.0; original.system().mode_count()];
        let drives = vec![drive(0.0, 0.0); 3];
        for _ in 0..64 {
            assert_eq!(
                original.step(&zero, &drives).unwrap(),
                moving.step(&zero, &drives).unwrap()
            );
            assert_eq!(states(&original), states(&moving));
            assert_eq!(
                original.system().last_obstacle_frames(),
                moving.system().last_obstacle_frames()
            );
            assert_eq!(
                moving
                    .fiber_support_frame()
                    .unwrap()
                    .normal_attachment_work_j,
                0.0
            );
        }
    });
}

#[test]
fn accelerating_hair_transports_the_normal_basis_and_closes_laboratory_work() {
    with_cx(|cx, gate| {
        let mut moving = normal_fiber_bank(cx, gate, 40.0, -1e-6)
            .with_fiber_support_motion(0.0, support_limits(), gate)
            .unwrap();
        let initial = moving.system().total_energy_j().unwrap();
        let fiber = moving.source().fiber_ribbon().unwrap().clone();
        let first = fiber.component_start();
        let per_cell =
            fiber.model().resolution().axial_modes + fiber.model().resolution().transverse_modes;
        let first_mode: usize = moving.system().components()[..first]
            .iter()
            .map(|c| c.modes().len())
            .sum();
        let scale = (2.0
            / (fiber.model().cell_line_density_kg_m() * fiber.model().geometry().length_m))
            .sqrt();
        let mut forces = vec![0.0; moving.system().mode_count()];
        let normal_force_map = fiber
            .model()
            .displacement_map(fiber.station_fraction(), FiberDirection::Normal)
            .unwrap();
        let mut previous_station = fiber.station_fraction();
        let (mut balance, mut normal_work, mut loaded, mut released) = (initial, 0.0, false, false);
        for sample in 1..=192 {
            let old: Vec<_> = moving.system().components()[first..first + 3]
                .iter()
                .map(|c| c.states().to_vec())
                .collect();
            let speed = if sample <= 96 {
                0.01 * sample as f64 / 96.0
            } else {
                0.01 - 0.02 * (sample - 96) as f64 / 96.0
            };
            // Apply a real normal force at a retained material point of each
            // fiber. Its solved deformation opens the gap; no translated normal
            // frame without the corresponding inertia is invented.
            if sample > 128 {
                for cell in 0..3 {
                    for (mode, shape) in normal_force_map.iter().enumerate() {
                        forces[first_mode + cell * per_cell + mode] = 0.05 * shape;
                    }
                }
            }
            let frame = moving
                .step(&forces, &vec![drive(speed, 0.0); 3])
                .unwrap()
                .clone();
            let station = moving.fiber_support_kinematics().unwrap().station_fraction;
            let normals = moving.system().last_obstacle_frames().unwrap();
            let mut expected_work = 0.0;
            for cell in 0..3 {
                assert_eq!(
                    frame.ports[cell].normal_force_n,
                    normals[cell].normal_force_n
                );
                loaded |= normals[cell].normal_force_n > 0.0;
                if normals[cell].normal_force_n == 0.0 {
                    assert_eq!(frame.ports[cell].traction_n, 0.0);
                    released = true;
                }
                let column = moving.system().obstacle_normal_column(cell).unwrap();
                for mode in 0..2 {
                    let n = (mode + 1) as f64;
                    // Independently evaluate the mass-normalized continuum row.
                    let end = (n * core::f64::consts::PI * station).sin() * scale;
                    let start = (n * core::f64::consts::PI * previous_station).sin() * scale;
                    let index = first_mode + cell * per_cell + 2 + mode;
                    assert!((column[index] + end).abs() < 1e-12 * scale);
                    assert_eq!(column[first_mode + cell * per_cell + mode], 0.0);
                    let next = moving.system().components()[first + cell].states()[2 + mode];
                    expected_work -= normals[cell].normal_force_n
                        * (end - start)
                        * 0.5
                        * (old[cell][2 + mode].displacement_m_sqrt_kg
                            + next.displacement_m_sqrt_kg);
                }
            }
            let support = moving.fiber_support_frame().unwrap();
            assert!((support.normal_attachment_work_j - expected_work).abs() < 1e-17);
            normal_work += support.normal_attachment_work_j.abs();
            // Inertial coordinate work is already in native external work;
            // normal map work is already in native obstacle work. Remove each
            // before substituting physical laboratory support work exactly once.
            balance += frame.external_work_j - support.inertial_coordinate_work_j
                + support.support_work_j
                + frame.obstacle_work_j
                - support.normal_attachment_work_j
                - frame.network_dissipation_j
                - frame.friction_dissipation_j
                - frame.obstacle_dissipation_j
                - frame.obstacle_residual_loss_j;
            assert!((support.laboratory_mechanical_energy_j - balance).abs() < 2e-9);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            previous_station = station;
        }
        assert!(loaded && released);
        assert!(
            normal_work > 1e-16,
            "moving normal fiber geometry did no work"
        );
    });
}

#[test]
fn a_refused_joint_normal_and_support_step_preserves_all_accepted_state() {
    with_cx(|cx, gate| {
        let make = || {
            normal_fiber_bank(cx, gate, 40.0, -1e-6)
                .with_fiber_support_motion(0.0, support_limits(), gate)
                .unwrap()
        };
        let mut actual = make();
        let mut expected = make();
        let zero = vec![0.0; actual.system().mode_count()];
        let drives = vec![drive(0.01, 0.0); 3];
        assert_eq!(
            actual.step(&zero, &drives).unwrap(),
            expected.step(&zero, &drives).unwrap()
        );
        let before = states(&actual);
        let accepted = actual.system().last_frame().cloned();
        let normals = actual.system().last_obstacle_frames().unwrap().to_vec();
        let maps: Vec<_> = (0..3)
            .map(|i| actual.system().obstacle_normal_column(i).unwrap().to_vec())
            .collect();
        let frame = actual.fiber_support_frame().copied();
        let support = actual.fiber_support_kinematics();
        // This refusal occurs after support/map preparation, in the original
        // normal translation owner: that extra motion was never admitted.
        // A subsequent valid retry must match a twin.
        assert!(
            actual
                .step_with_stop_offsets(&zero, &drives, &[0.002; 3])
                .is_err()
        );
        assert_eq!(states(&actual), before);
        assert_eq!(actual.system().last_frame(), accepted.as_ref());
        assert_eq!(actual.system().last_obstacle_frames().unwrap(), normals);
        assert_eq!(actual.fiber_support_kinematics(), support);
        assert_eq!(actual.fiber_support_frame().copied(), frame);
        for (i, map) in maps.iter().enumerate() {
            assert_eq!(actual.system().obstacle_normal_column(i).unwrap(), map);
        }
        assert_eq!(
            actual.step(&zero, &drives).unwrap(),
            expected.step(&zero, &drives).unwrap()
        );
        assert_eq!(states(&actual), states(&expected));
        assert_eq!(actual.fiber_support_frame(), expected.fiber_support_frame());
    });
}

#[test]
fn continuous_bow_travel_solved_normal_contact_and_fiber_reversal_share_the_same_sample_clock() {
    with_cx(|cx, gate| {
        let make = |tension| {
            let bank = normal_fiber_bank(cx, gate, tension, -1e-6);
            let ribbon = bank.source().ribbon().unwrap();
            let length = bank.cards()[0].length_m;
            let stations = |centre| {
                UniformBowRibbon::try_new(
                    StringBowLocation {
                        string: 0,
                        station_fraction: centre,
                    },
                    length,
                    ribbon.width_m(),
                    3,
                )
                .unwrap()
                .stations()
                .to_vec()
            };
            let bow_path = BowStationProgram {
                knots: vec![
                    BowStationKnot {
                        sample: 256,
                        station_fractions: stations(0.201),
                    },
                    BowStationKnot {
                        sample: 384,
                        station_fractions: stations(0.199),
                    },
                ],
                limits: BowStationMotionLimits {
                    maximum_speed_m_s: 1.0,
                    maximum_step_m: 0.001,
                    max_update_terms: 1_048_576,
                },
                max_knots: 4,
                max_coordinate_values: 12,
            };
            bank.into_scheduled_ribbon(fiber_load_schedule(0.0), "bow".into(), performance_config())
                .unwrap()
                .with_bow_station_motion(bow_path, gate)
                .unwrap()
                .with_fiber_support_program(
                    support_program(&[(128, 0.02), (256, -0.02), (384, 0.0)]),
                    gate,
                )
                .unwrap()
        };
        let mut one = make(10.0);
        let mut chunks = make(10.0);
        let mut taut = make(100.0);
        let mut a = vec![0.0; 513];
        let mut b = a.clone();
        let mut c = a.clone();
        for block in a.chunks_mut(1) {
            one.block(block).unwrap();
        }
        for block in b.chunks_mut(37) {
            chunks.block(block).unwrap();
        }
        for block in c.chunks_mut(128) {
            taut.block(block).unwrap();
        }
        assert_eq!(a, b);
        assert_eq!(states(one.bank()), states(chunks.bank()));
        assert_eq!(one.fiber_support_frame(), chunks.fiber_support_frame());
        assert_eq!(
            one.bank().bow_station_fractions(),
            one.bank().stop_station_fractions()
        );
        assert_eq!(one.bank().system().samples_rendered(), 513);
        assert!(a.iter().any(|x| x.abs() > 1e-12));
        assert!(a.iter().zip(c).any(|(x, y)| (x - y).abs() > 1e-12));
    });
}
