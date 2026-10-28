//! Source material admission around the original joint contact transaction.
use super::*;
use super::super::geometry::{ObstacleGeometry, ObstacleGeometryConfig};
use super::super::motion::ObstacleMotionLimits;
use fs_contact::normal_patch::source::{PointContactGeometry, SourcePointContact};
#[path = "../../../../../../../../../../fs-contact/tests/source_point/common.rs"]
mod source_fixture;

fn source(young: f64, rate: f64) -> SourcePointContact {
    let a = source_fixture::solid("A", 1000.0, young);
    let b = source_fixture::solid("B", 2000.0, 2e7);
    let mut envelope = source_fixture::envelope();
    envelope.maximum_absolute_rate_m_per_s = rate;
    SourcePointContact::new(source_fixture::bind(&a, &b, Some(0.2), 0.0),
        PointContactGeometry::SpherePlane { effective_radius_m: 0.1 }, envelope).unwrap()
}
fn contact(source: &SourcePointContact, gap: f64) -> (ModalContact, ModalContactConfig) {
    let mut limits = normal_limits();
    limits.maximum_penetration_m = source.envelope().maximum_indentation_m;
    (ModalContact { left: attachment(0, 1.0), right: attachment(0, 0.0),
        law: crate::unilateral_contact::binding::bound_point_obstacle(
            source, vec![-1.0], gap, "synthetic sourced support").unwrap() }, limits)
}
fn endpoint(shape: f64) -> ObstacleGeometry {
    ObstacleGeometry { left: attachment(0, shape), right: attachment(0, 0.0) }
}

#[test]
fn source_modulus_reaches_the_same_hertz_potential_and_joint_normal_reaction() {
    let gate = CancelGate::new_clock_free();
    let mut reactions = Vec::new();
    for young in [1e7, 4e7] {
        let source = source(young, 0.2);
        let k = source.stiffness_n_per_m_3_2();
        let indentation = 1e-5_f64;
        let mut system = multi(friction(network(1e6)))
            .with_obstacles(vec![contact(&source, -indentation)], limits(), &gate).unwrap()
            .with_obstacle_sources(vec![(0, source.clone())], &gate).unwrap();
        assert_eq!(system.obstacle_source(0), Some(&source));
        // A sticking bow balances the real normal reaction, holding q exactly
        // at its original position; no fitted modal frequency enters this check.
        for _ in 0..16 {
            system.step(&[0.0], &[FrictionDrive { speed_m_s: 0.0, normal_force_n: 1.0 }]).unwrap();
            let normal = system.last_obstacle_frames().unwrap()[0];
            assert!((normal.normal_force_n-k*indentation.powf(1.5)).abs() < 1e-9);
            assert!((normal.stored_energy_j-0.4*k*indentation.powf(2.5)).abs() < 1e-14);
            assert!(normal.dissipation_j < 1e-16);
        }
        reactions.push(system.last_obstacle_frames().unwrap()[0].normal_force_n);
    }
    assert!(reactions[1] > reactions[0]);
}

#[test]
fn moving_map_and_support_source_rate_refusal_preserves_state_work_and_retry() {
    let gate = CancelGate::new_clock_free();
    let source = source(1e7, 0.01);
    let build = |bound: bool| {
        let mut component = model(100.0, 0.0);
        component.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: 2e-5,
            velocity_m_sqrt_kg_per_s: 0.0 }]).unwrap();
        let net = CoupledModalSystem::new(vec![component], vec![], config(1e6), &gate).unwrap();
        let system = multi(friction(net))
            .with_obstacles(vec![contact(&source, -1e-5)], limits(), &gate).unwrap();
        let system = if bound { system.with_obstacle_sources(vec![(0, source.clone())], &gate).unwrap() }
            else { system };
        system.with_obstacle_motion(ObstacleMotionLimits {
            maximum_abs_offset_m: 0.001, maximum_speed_m_s: 1.0 }, &gate).unwrap()
            .with_obstacle_geometry_motion(ObstacleGeometryConfig { max_update_terms: 100000 }, &gate).unwrap()
    };
    let mut actual = build(true);
    let mut reference = build(false);
    let before = states(&actual);
    // Outer nonlinear/thermal roots may visit this mathematical candidate. The
    // material guard belongs after root acceptance, not inside stage_constitutive.
    actual.prepare_step(&[released()], None, Some(&[endpoint(1.1)]), false, Some(&gate)).unwrap();
    actual.stage_constitutive(&[0.0], &[released()], None, Some(&gate)).unwrap();
    assert!(actual.validate_obstacle_sources(Some(&gate)).is_err());
    assert_eq!(states(&actual), before);
    for (map, offset) in [(1.1, 0.0), (1.0, -1e-6)] {
        let error = actual.step_with_obstacle_geometry(&[0.0], &[released()],
            &[endpoint(map)], Some(&[offset])).unwrap_err();
        assert!(error.to_string().contains("source point absolute rate"), "{error}");
        assert_eq!(states(&actual), before);
        assert_eq!(actual.samples_rendered(), 0);
        assert!(actual.last_frame().is_none());
        assert_eq!(actual.obstacle_offsets_m().unwrap(), &[0.0]);
        assert_eq!(actual.obstacle_normal_column(0).unwrap(), &[1.0]);
    }
    let (mut old_map, mut old_offset) = (1.0, 0.0);
    for step in 1..=32 {
        let map = 1.0+0.001*f64::from(step);
        let offset = -1e-8*f64::from(step);
        let old_q = actual.components()[0].states()[0].displacement_m_sqrt_kg;
        assert_eq!(actual.step_with_obstacle_geometry(&[0.0], &[released()], &[endpoint(map)], Some(&[offset])).unwrap(),
            reference.step_with_obstacle_geometry(&[0.0], &[released()], &[endpoint(map)], Some(&[offset])).unwrap());
        assert_eq!(states(&actual), states(&reference));
        let normal = actual.last_obstacle_frames().unwrap()[0];
        let q = actual.components()[0].states()[0].displacement_m_sqrt_kg;
        let work = normal.normal_force_n*((map-old_map)*0.5*(q+old_q)-(offset-old_offset));
        assert!((normal.actuator_work_j-work).abs() < 1e-16);
        assert!(actual.last_frame().unwrap().energy_residual_j.abs()
            <= actual.last_frame().unwrap().energy_tolerance_j);
        (old_map, old_offset) = (map, offset);
    }
}

#[test]
fn mismatched_coefficients_point_weight_initial_rate_and_clock_refuse() {
    let gate = CancelGate::new_clock_free();
    let source = source(1e7, 0.01);
    let base = || multi(friction(network(1e6)));
    let (original, limits_normal) = contact(&source, 0.0);
    for (k, weight) in [(source.stiffness_n_per_m_3_2()*2.0, 1.0),
        (source.stiffness_n_per_m_3_2(), 2.0)] {
        let mut changed = original.clone();
        changed.law = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![weight],
            k, 1.5, "synthetic mismatched source".into()).unwrap().with_internal_loss(0.2).unwrap();
        assert!(base().with_obstacles(vec![(changed, limits_normal)], limits(), &gate).unwrap()
            .with_obstacle_sources(vec![(0, source.clone())], &gate).is_err());
    }
    assert!(base().with_obstacles(vec![contact(&source, 0.0)], limits(), &gate).unwrap()
        .with_obstacle_sources(vec![(0, source.clone()), (0, source.clone())], &gate).is_err());
    let net = CoupledModalSystem::new(vec![model(100.0, 0.02)], vec![], config(1e6), &gate).unwrap();
    assert!(multi(friction(net)).with_obstacles(vec![contact(&source, 0.0)], limits(), &gate).unwrap()
        .with_obstacle_sources(vec![(0, source.clone())], &gate).is_err());
    let mut envelope = source.envelope();
    envelope.maximum_step_s = 1e-6;
    let short = SourcePointContact::new(source.model().clone(), source.geometry(), envelope).unwrap();
    assert!(base().with_obstacles(vec![contact(&short, 0.0)], limits(), &gate).unwrap()
        .with_obstacle_sources(vec![(0, short)], &gate).is_err());
}
