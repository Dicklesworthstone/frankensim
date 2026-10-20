//! Independent unilateral storage/product-rule work, observation and rollback.
use super::*;
use super::super::geometry::{ObstacleGeometry, ObstacleGeometryConfig};
use super::super::motion::ObstacleMotionLimits;

fn admission() -> ObstacleGeometryConfig { ObstacleGeometryConfig { max_update_terms: 100_000 } }
fn maps(left: f64, right: f64) -> ObstacleGeometry {
    ObstacleGeometry { left: attachment(0, left), right: attachment(1, right) }
}
fn value(system: &MultiFrictionModalSystem, i: usize) -> f64 {
    system.components()[i].states()[0].displacement_m_sqrt_kg
}
fn stored(x: f64) -> f64 { 1000.0*x.max(0.0).powi(2) }

#[test]
fn constant_geometry_preserves_original_joint_arithmetic_exactly() {
    let gate = CancelGate::new_clock_free();
    let build = || {
        let mut stop = obstacle(1.0, 1e-5); stop.right = attachment(1, 0.8);
        multi(friction(coupled_network()).with_acceleration_pressure(vec![1.0, -0.4]).unwrap())
            .with_obstacles(vec![(stop, normal_limits())], limits(), &gate).unwrap()
    };
    let mut fixed = build();
    let mut moving = build().with_obstacle_geometry_motion(admission(), &gate).unwrap();
    for _ in 0..128 {
        let drive = FrictionDrive { speed_m_s: 0.03, normal_force_n: 0.1 };
        assert_eq!(fixed.step(&[0.0, 0.0], &[drive]).unwrap(),
            moving.step_with_obstacle_geometry(&[0.0, 0.0], &[drive], &[maps(1.0, 0.8)], None).unwrap());
        assert_eq!(states(&fixed), states(&moving));
        assert_eq!(fixed.last_obstacle_frames(), moving.last_obstacle_frames());
        assert_eq!(moving.last_obstacle_frames().unwrap()[0].attachment_work_j, 0.0);
    }
}

#[test]
fn moving_attachment_and_support_close_actual_work_and_endpoint_pressure() {
    let gate = CancelGate::new_clock_free();
    let mut models = vec![model(300.0, 0.02), model(500.0, -0.005)];
    for (m, displacement) in models.iter_mut().zip([3e-5, -1e-5]) {
        let mut state = m.states()[0]; state.displacement_m_sqrt_kg = displacement;
        m.restore_states(&[state]).unwrap();
    }
    let network = CoupledModalSystem::new(models, vec![], config(1e6), &gate).unwrap();
    let core = friction(network).with_acceleration_pressure(vec![1.0, -0.4]).unwrap();
    let per_port = core.config;
    let mut stop = obstacle(1.0, 0.0); stop.right = attachment(1, 0.8);
    stop.law = stop.law.with_internal_loss(0.2).unwrap();
    let mut system = core.with_additional_ports(vec![(ModalFrictionPort {
        left: attachment(1, 1.0), right: None,
        law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 },
    }, per_port)], MultiFrictionConfig { max_ports: 2, max_sweeps: 128, max_setup_terms: 100_000 }, &gate).unwrap()
        .with_obstacles(vec![(stop, normal_limits())], limits(), &gate).unwrap()
        .with_obstacle_motion(ObstacleMotionLimits { maximum_abs_offset_m: 0.001, maximum_speed_m_s: 1.0 }, &gate).unwrap()
        .with_obstacle_geometry_motion(admission(), &gate).unwrap();
    let initial = system.total_energy_j().unwrap();
    let (mut supply, mut loss, mut old_offset) = (0.0, 0.0, 0.0);
    let mut old_map = [1.0, -0.8];
    let (mut separated, mut recontact, mut migration_work) = (false, false, 0.0);
    let dt = system.sample_period_s();
    for sample in 0..192 {
        let old_q = [value(&system, 0), value(&system, 1)];
        let next_map = [1.0+0.001*(sample+1) as f64, -0.8+0.0005*(sample+1) as f64];
        let offset = if sample < 32 { 0.0 } else if sample < 96 { 4e-6*(sample-31) as f64 }
            else if sample < 128 { 256e-6 } else { 4e-6*(191-sample) as f64 };
        let external = [0.03, -0.01];
        let drives = [FrictionDrive { speed_m_s: 0.2, normal_force_n: 0.15 },
            FrictionDrive { speed_m_s: -0.12, normal_force_n: 0.1 }];
        let frame = system.step_with_obstacle_geometry(&external, &drives,
            &[maps(next_map[0], -next_map[1])], Some(&[offset])).unwrap().clone();
        let contact = system.last_obstacle_frames().unwrap()[0];
        let q = [value(&system, 0), value(&system, 1)];
        let delta = [q[0]-old_q[0], q[1]-old_q[1]];
        let old_x = old_map[0]*old_q[0]+old_map[1]*old_q[1]-old_offset;
        let x = next_map[0]*q[0]+next_map[1]*q[1]-offset;
        let dx = x-old_x; let du = stored(x)-stored(old_x);
        // Analytic unilateral quadratic spring integral, including separation.
        let elastic = if x > 0.0 && old_x > 0.0 { 1000.0*(x+old_x) }
            else if dx == 0.0 { 2000.0*x.max(0.0) } else { du/dx };
        let reaction = (elastic*(1.0+0.2*dx/dt)).max(0.0);
        assert!((contact.normal_force_n-reaction).abs() < 1e-9);
        assert!((contact.stored_energy_j-stored(x)).abs() < 1e-15);
        assert!((contact.dissipation_j-(contact.normal_force_n*dx-du)).abs() < 1e-12);
        let attachment_work = contact.normal_force_n * (0..2).map(|i|
            (next_map[i]-old_map[i])*0.5*(q[i]+old_q[i])).sum::<f64>();
        let actuator_work = attachment_work-contact.normal_force_n*(offset-old_offset);
        let mechanical_work = -contact.normal_force_n*(0..2).map(|i|
            0.5*(next_map[i]+old_map[i])*delta[i]).sum::<f64>();
        assert!((contact.attachment_work_j-attachment_work).abs() < 1e-16);
        assert!((frame.obstacle_work_j-actuator_work).abs() < 1e-16);
        assert!((mechanical_work+du+contact.dissipation_j-actuator_work).abs() < 1e-12);
        let external_work = external[0]*delta[0]+external[1]*delta[1];
        let bow_work = (0..2).map(|i| frame.ports[i].traction_n*drives[i].speed_m_s*dt).sum::<f64>();
        assert!((frame.external_work_j-external_work).abs() < 1e-16);
        assert!((frame.drive_work_j-bow_work).abs() < 1e-16);
        let pressure = (0..2).map(|i| {
            let mode = system.components()[i].modes()[0];
            let v = system.components()[i].states()[0].velocity_m_sqrt_kg_per_s;
            [1.0, -0.4][i]*(external[i]+frame.ports[i].traction_n-next_map[i]*contact.normal_force_n
                -2.0*mode.damping_ratio*mode.angular_frequency_rad_s*v-mode.angular_frequency_rad_s.powi(2)*q[i])
        }).sum::<f64>();
        assert!((frame.observer_pressure_pa-pressure).abs() < 1e-10);
        supply += external_work+bow_work+actuator_work;
        loss += frame.network_dissipation_j+frame.friction_dissipation_j+contact.dissipation_j;
        assert!((system.total_energy_j().unwrap()-initial+loss-supply).abs() < 1e-9);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert_eq!(system.obstacle_normal_column(0).unwrap(), &next_map);
        separated |= contact.normal_force_n == 0.0;
        recontact |= sample >= 128 && contact.normal_force_n > 0.0;
        migration_work += attachment_work.abs();
        old_map = next_map; old_offset = offset;
    }
    assert!(separated && recontact && migration_work > 0.0);
    assert_eq!(system.obstacle(0).unwrap().left.shapes, vec![1.0], "immutable source map changed");
    // HOLD must use the accepted endpoint, never the preceding moving midpoint.
    let frame = system.step(&[0.03, -0.01], &[released(), released()]).unwrap();
    assert_eq!(frame.obstacle_work_j, 0.0);
    assert_eq!(system.obstacle_normal_column(0).unwrap(), &old_map);
}

#[test]
fn failed_geometry_pressure_and_cancel_are_atomic_and_retry_exactly() {
    let gate = CancelGate::new_clock_free();
    let build = || multi(friction(network(0.001)).with_acceleration_pressure(vec![1.0]).unwrap())
        .with_obstacles(vec![(obstacle(1.0, 0.0), normal_limits())], limits(), &gate).unwrap()
        .with_obstacle_motion(ObstacleMotionLimits { maximum_abs_offset_m: 0.001, maximum_speed_m_s: 1.0 }, &gate).unwrap()
        .with_obstacle_geometry_motion(admission(), &gate).unwrap();
    let endpoint = |shape| ObstacleGeometry { left: attachment(0, shape), right: attachment(0, 0.0) };
    let mut actual = build(); let mut reference = build();
    actual.step(&[0.0], &[released()]).unwrap(); reference.step(&[0.0], &[released()]).unwrap();
    let state = states(&actual); let report = actual.last_frame().unwrap().clone();
    let normals = actual.last_obstacle_frames().unwrap().to_vec();
    let error = actual.step_with_obstacle_geometry(&[0.0], &[released()], &[endpoint(1.2)], Some(&[-1e-5])).unwrap_err();
    assert!(error.to_string().contains("pressure"), "{error}");
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(actual.step_with_obstacle_geometry_under_gate(&[0.0], &[released()], &[endpoint(1.2)], None, &cancelled).is_err());
    for bad in [endpoint(-1.0), endpoint(f64::NAN), endpoint(f64::INFINITY)] {
        assert!(actual.step_with_obstacle_geometry(&[0.0], &[released()], &[bad], None).is_err());
    }
    assert_eq!(states(&actual), state); assert_eq!(actual.last_frame(), Some(&report));
    assert_eq!(actual.last_obstacle_frames().unwrap(), normals);
    assert_eq!(actual.obstacle_normal_column(0).unwrap(), &[1.0]);
    assert_eq!(actual.obstacle_offsets_m().unwrap(), &[0.0]);
    for step in 0..32 {
        let geometry = [endpoint(1.0+0.001*f64::from(step))];
        assert_eq!(actual.step_with_obstacle_geometry(&[0.0001], &[released()], &geometry, None).unwrap(),
            reference.step_with_obstacle_geometry(&[0.0001], &[released()], &geometry, None).unwrap());
        assert_eq!(states(&actual), states(&reference));
    }
    for _ in 0..16 {
        assert_eq!(actual.step(&[0.0001], &[released()]).unwrap(), reference.step(&[0.0001], &[released()]).unwrap());
        assert_eq!(states(&actual), states(&reference));
    }
    assert!(build().with_obstacle_geometry_motion(admission(), &gate).is_err());
    assert!(multi(friction(network(1.0))).with_obstacle_geometry_motion(admission(), &gate).is_err());
    let unadmitted = || multi(friction(network(1.0)))
        .with_obstacles(vec![(obstacle(1.0, 0.0), normal_limits())], limits(), &gate).unwrap();
    assert!(unadmitted().step_with_obstacle_geometry(&[0.0], &[released()], &[endpoint(1.0)], None).is_err());
    assert!(unadmitted().with_obstacle_geometry_motion(ObstacleGeometryConfig { max_update_terms: 1 }, &gate).is_err());
}
