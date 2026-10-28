//! Normal load, moving material crossings and support work share one solve.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    performance::bow_stations::{BowStationKnot, BowStationMotionLimits, BowStationProgram},
    ribbon::UniformBowRibbon,
};

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
