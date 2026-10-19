//! G0/G1 work, translation and atomicity on the original contact/friction owner.
use super::*;
use crate::render::schedule::force::coupled::contact::friction::multiple::obstacles::motion::ObstacleMotionLimits;

fn motion_limits() -> ObstacleMotionLimits {
    ObstacleMotionLimits { maximum_abs_offset_m: 0.01, maximum_speed_m_s: 1.0 }
}
fn moving(gap: f64, pressure: f64) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    multi(friction(network(pressure)).with_acceleration_pressure(vec![1.0]).unwrap())
        .with_obstacles(vec![(obstacle(1.0, gap), normal_limits())], limits(), &gate).unwrap()
        .with_obstacle_motion(motion_limits(), &gate).unwrap()
}

#[test]
fn zero_motion_matches_the_original_fixed_support_path() {
    let gate = CancelGate::new_clock_free();
    let build = || multi(friction(coupled_network()))
        .with_obstacles(vec![(obstacle(1.0, 1e-5), normal_limits())], limits(), &gate).unwrap();
    let mut a = build();
    let mut b = build().with_obstacle_motion(motion_limits(), &gate).unwrap();
    for _ in 0..128 {
        let drive = FrictionDrive { speed_m_s: 0.03, normal_force_n: 0.1 };
        let expected = a.step(&[0.0, 0.0], &[drive]).unwrap().clone();
        let actual = b.step_with_obstacle_offsets(&[0.0, 0.0], &[drive], &[0.0]).unwrap().clone();
        assert_eq!(states(&a), states(&b));
        assert_eq!(expected, actual);
        assert_eq!(a.last_obstacle_frames(), b.last_obstacle_frames());
    }
}

#[test]
fn motion_work_equals_stored_energy_change_with_a_stationary_sticking_bow() {
    let mut system = moving(-1e-4, 1000.0);
    let hold = FrictionDrive { speed_m_s: 0.0, normal_force_n: 1.0 };
    let initial = system.total_energy_j().unwrap();
    let mut supply = 0.0;
    let mut previous = 0.0;
    // Press and then withdraw. Holding the surface at rest isolates the
    // manufactured potential U=k*(1e-4-offset)^2/2 from oscillator motion.
    for endpoint in [-1e-5, -2e-5, -3e-5, -2e-5, -1e-5, 0.0] {
        let frame = system.step_with_obstacle_offsets(&[0.0], &[hold], &[endpoint]).unwrap().clone();
        let normal = system.last_obstacle_frames().unwrap()[0];
        let expected = 1000.0*((1e-4_f64-endpoint).powi(2)-(1e-4_f64-previous).powi(2));
        assert!((frame.obstacle_work_j-expected).abs() < 1e-12);
        assert!((normal.actuator_work_j + normal.normal_force_n*(endpoint-previous)).abs() < 1e-18);
        assert_eq!(normal.support_offset_m, endpoint);
        assert_eq!(frame.drive_work_j, 0.0);
        assert_eq!(frame.external_work_j, 0.0);
        assert_eq!(frame.ports[0].regime, ModalFrictionRegime::Sticking);
        supply += frame.obstacle_work_j;
        assert!((system.total_energy_j().unwrap()-initial-supply).abs() < 1e-11);
        assert!(frame.observer_pressure_pa.abs() < 1e-8);
        previous = endpoint;
    }
    assert!(supply.abs() < 1e-12);
    assert_eq!(system.obstacle(0).unwrap().law.gaps(), &[-1e-4]);
}

#[test]
fn normal_motion_loads_both_attachment_bodies_and_holding_does_not_snap_back() {
    let gate = CancelGate::new_clock_free();
    let bodies = CoupledModalSystem::new(vec![model(300.0, 0.0), model(500.0, 0.0)],
        vec![], config(1e6), &gate).unwrap();
    let mut normal = obstacle(1.0, 0.0);
    normal.right = attachment(1, 0.8);
    normal.law = normal.law.with_internal_loss(0.2).unwrap();
    let mut system = multi(friction(bodies)).with_obstacles(vec![(normal, normal_limits())], limits(), &gate).unwrap()
        .with_obstacle_motion(motion_limits(), &gate).unwrap();
    let frame = system.step_with_obstacle_offsets(&[0.0, 0.0], &[released()], &[-1e-5]).unwrap().clone();
    assert!(frame.obstacle_work_j > 0.0);
    assert!(frame.obstacle_dissipation_j > 0.0);
    assert!(system.components()[0].states()[0].velocity_m_sqrt_kg_per_s < 0.0);
    assert!(system.components()[1].states()[0].velocity_m_sqrt_kg_per_s > 0.0);
    let mut energy = system.total_energy_j().unwrap();
    for _ in 0..256 {
        let frame = system.step(&[0.0, 0.0], &[released()]).unwrap().clone();
        assert_eq!(system.obstacle_offsets_m().unwrap(), &[-1e-5]);
        assert_eq!(frame.obstacle_work_j, 0.0);
        let next = system.total_energy_j().unwrap();
        assert!(next <= energy+frame.energy_tolerance_j);
        energy = next;
    }
    assert!(energy > 0.0);
}

#[test]
fn rejected_pressure_and_cancelled_motion_leave_positions_and_reports_untouched() {
    let mut system = moving(0.0, 1e-6);
    system.step(&[0.0], &[released()]).unwrap();
    let before = states(&system);
    let frame = system.last_frame().unwrap().clone();
    let normals = system.last_obstacle_frames().unwrap().to_vec();
    assert!(system.step_with_obstacle_offsets(&[0.0], &[released()], &[-1e-5]).is_err());
    assert_eq!(states(&system), before);
    assert_eq!(system.last_frame(), Some(&frame));
    assert_eq!(system.last_obstacle_frames().unwrap(), normals);
    assert_eq!(system.obstacle_offsets_m().unwrap(), &[0.0]);
    let stop = CancelGate::new_clock_free(); stop.request();
    assert!(system.step_with_obstacle_offsets_under_gate(&[0.0], &[released()], &[-1e-5], &stop).is_err());
    assert_eq!(system.obstacle_offsets_m().unwrap(), &[0.0]);
    // Scratch endpoint from the rejected trial is not a pending command.
    assert_eq!(system.step(&[0.0], &[released()]).unwrap().observer_pressure_pa, 0.0);
    assert_eq!(system.samples_rendered(), 2);
}

#[test]
fn absent_unbound_duplicate_and_out_of_budget_motion_refuse() {
    let gate = CancelGate::new_clock_free();
    assert!(multi(friction(network(1000.0))).with_obstacle_motion(motion_limits(), &gate).is_err());
    let mut unbound = multi(friction(network(1000.0)))
        .with_obstacles(vec![(obstacle(1.0, 0.0), normal_limits())], limits(), &gate).unwrap();
    assert!(unbound.step_with_obstacle_offsets(&[0.0], &[released()], &[0.0]).is_err());
    let mut system = moving(0.0, 1000.0);
    for endpoint in [vec![], vec![0.0, 0.0], vec![f64::NAN], vec![0.02], vec![0.0001]] {
        assert!(system.step_with_obstacle_offsets(&[0.0], &[released()], &endpoint).is_err());
        assert_eq!(system.samples_rendered(), 0);
        assert_eq!(system.obstacle_offsets_m().unwrap(), &[0.0]);
    }
    assert!(system.with_obstacle_motion(motion_limits(), &gate).is_err());
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(ObstacleMotionLimits { maximum_abs_offset_m: value, maximum_speed_m_s: 1.0 }.validate().is_err());
        assert!(ObstacleMotionLimits { maximum_abs_offset_m: 1.0, maximum_speed_m_s: value }.validate().is_err());
    }
}
