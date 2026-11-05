//! Real finite-body motion continues on the same event-timed caloric history.
use super::*;
use fs_couple::friction_heat::adhesive::snap::native::*;
use fs_mbd::{
    MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench,
    tangential_friction::{RigidContactBinding, adhesive::jkr_normal_motion::*},
};
use fs_tribo::partial_slip::TangentFrame;

#[derive(Clone)]
struct NativePair {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
}
impl NativePair {
    fn collision(source: &SourceJkrSphereContact) -> Self {
        let (h, _, v, _) = scales(source);
        let properties = [MassProperties::new(
            2.0 * RELATIVE_MASS,
            Vec3::ZERO,
            Vec3::new(1e-9, 1.2e-9, 1.4e-9),
        )
        .unwrap(); 2];
        let bodies = std::array::from_fn(|side| {
            RigidBodyState::new(
                Pose::new(
                    Vec3::new(0.0, 0.0, side as f64 * 0.25 * h),
                    UnitQuaternion::IDENTITY,
                )
                .unwrap(),
                Vec3::new(0.0, 0.0, if side == 0 { 10.0 * v } else { 0.0 })
                    .scale(properties[side].mass()),
                Vec3::ZERO,
            )
            .unwrap()
        });
        Self { bodies, properties }
    }
    fn prepare<'a>(&self, source: &'a SourceJkrSphereContact) -> PreparedRigidJkrPair<'a> {
        let (h, force, speed, _) = scales(source);
        PreparedRigidJkrPair::new(
            &self.bodies,
            &self.properties,
            RigidJkrContact {
                binding: RigidContactBinding {
                    body_a: 0,
                    body_b: 1,
                    arm_a_body: Vec3::ZERO,
                    arm_b_body: Vec3::ZERO,
                },
                frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
                reference_approach_m: 0.0,
                normal: source,
            },
            RigidJkrControls {
                motion: motion(source).budget(),
                maximum_bodies: 2,
                geometry_tolerance_m: 1e-10 * h,
                displacement_tolerance_m: 1e-9 * h,
                velocity_tolerance_m_s: 1e-9 * speed,
                linear_momentum_tolerance_ns: 1e-12 * RELATIVE_MASS * speed,
                angular_momentum_tolerance_nms: 1e-12 * RELATIVE_MASS * speed * h,
                energy_tolerance_j: 1e-8 * force * h,
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    }
    fn consume_travel(&mut self, after: [RigidBodyState; 2], report: &RigidJkrReport) {
        for side in 0..2 {
            assert_eq!(report.angular_displacement_body_rad[side], Vec3::ZERO);
            let body = after[side];
            self.bodies[side] = RigidBodyState::new(
                Pose::new(
                    body.pose()
                        .position_world()
                        .add(report.body_displacement_world_m[side]),
                    body.pose().orientation(),
                )
                .unwrap(),
                body.linear_momentum_world(),
                body.angular_momentum_body(),
            )
            .unwrap();
        }
    }
}

fn native_controls(source: &SourceJkrSphereContact) -> RigidJkrTimedSnapHeatControls {
    let (h, force, _, _) = scales(source);
    RigidJkrTimedSnapHeatControls {
        heat: controls(source),
        total_energy_tolerance_j: 2e-8 * force * h + controls(source).total_energy_tolerance_j,
    }
}

#[test]
fn loaded_native_collision_continues_one_snap_history_and_real_thermal_clock() {
    let source = source();
    let caloric = curves(&source, [1.0; 2]);
    let mut pair = NativePair::collision(&source);
    let prepared = pair.prepare(&source);
    let (_, force, _, time) = scales(&source);
    let origin = 0.125 * time;
    let thermal = thermal(&source, &caloric, origin, vec![], 0.0);
    let initial = source.detached_state(prepared.approach_m()).unwrap();
    let heat = binding(&source, &caloric, &thermal, initial);
    // Equal body accelerations add real kinetic energy and external work while
    // leaving the relative JKR trajectory unchanged.
    let wrenches = [Wrench {
        force_world: Vec3::new(0.02 * force, 0.0, 0.0),
        torque_body: Vec3::ZERO,
    }; 2];
    let first = heat
        .propose_timed_rigid_motion(
            &thermal,
            heat.stamp(),
            &prepared,
            wrenches,
            0.1 * time,
            origin,
            native_controls(&source),
            &CancelGate::new_clock_free(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(first.report().mechanics.motion.event_count(), 1);
    assert_eq!(first.report().heat.events[0].jump.kind, JkrJumpKind::Attach);
    assert_eq!(
        first.thermal().cells(),
        manual(&source, &thermal, &first.report().mechanics.motion, origin).cells()
    );
    let state = first.motion_state();
    let first_release = first.report().heat.snap.work.dissipated_work_j();
    let (after, heat, thermal, report) = first.into_parts();
    pair.consume_travel(after, &report.mechanics);
    let second = heat
        .propose_timed_rigid_motion(
            &thermal,
            heat.stamp(),
            &pair.prepare(&source),
            wrenches,
            7.9 * time,
            origin,
            native_controls(&source),
            &CancelGate::new_clock_free(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(second.report().mechanics.motion.initial_state(), state);
    assert_eq!(second.motion_state().steps(), state.steps() + 1);
    assert_eq!(second.report().mechanics.motion.event_count(), 1);
    assert_eq!(
        second.report().heat.events[0].jump.kind,
        JkrJumpKind::Detach
    );
    let event = second.report().mechanics.motion.events().next().unwrap();
    assert_eq!(
        second.report().heat.events[0].thermal_time_s,
        origin + (state.time_s() + event.offset_s)
    );
    assert_eq!(
        second.thermal().cells(),
        manual(&source, &thermal, &second.report().mechanics.motion, origin).cells()
    );
    let (after, heat, thermal, report) = second.into_parts();
    assert!(
        after
            .iter()
            .all(|body| body.linear_momentum_world().x > 0.0)
    );
    assert_eq!(
        thermal.time_s(),
        origin + report.mechanics.motion.state().time_s()
    );
    let expected_release = first_release + report.heat.snap.work.dissipated_work_j();
    assert!((heat.ledger().dissipated_work_j() - expected_release).abs() < 1e-22);
    assert!(
        report.combined_energy_residual_j.abs()
            <= native_controls(&source).total_energy_tolerance_j
    );
    assert!(
        report.absolute_stage_energy_residual_j
            <= native_controls(&source).total_energy_tolerance_j
    );
}

#[test]
fn continued_heat_refuses_stale_native_geometry_mass_and_stamp_without_advancing() {
    let source = source();
    let caloric = curves(&source, [1.0; 2]);
    let mut pair = NativePair::collision(&source);
    let stale_geometry = pair.prepare(&source);
    let time = scales(&source).3;
    let thermal = thermal(&source, &caloric, 0.0, vec![], 0.0);
    let initial = source.detached_state(stale_geometry.approach_m()).unwrap();
    let heat = binding(&source, &caloric, &thermal, initial);
    let old_stamp = heat.stamp();
    let first = heat
        .propose_timed_rigid_motion(
            &thermal,
            old_stamp,
            &stale_geometry,
            [Wrench::ZERO; 2],
            0.1 * time,
            0.0,
            native_controls(&source),
            &CancelGate::new_clock_free(),
            || Ok(()),
        )
        .unwrap();
    let (after, heat, thermal, report) = first.into_parts();
    pair.consume_travel(after, &report.mechanics);
    let saved_cells = thermal.cells().clone();
    let saved_stamp = heat.stamp();
    let saved_bodies = pair.bodies;
    let attempt = |native: &PreparedRigidJkrPair<'_>, stamp| {
        heat.propose_timed_rigid_motion(
            &thermal,
            stamp,
            native,
            [Wrench::ZERO; 2],
            0.1 * time,
            0.0,
            native_controls(&source),
            &CancelGate::new_clock_free(),
            || Ok(()),
        )
    };
    assert!(matches!(
        attempt(&stale_geometry, saved_stamp),
        Err(RigidJkrTimedSnapHeatError::Mechanics(_))
    ));
    assert!(matches!(
        attempt(&pair.prepare(&source), old_stamp),
        Err(RigidJkrTimedSnapHeatError::Heat(
            JkrSnapHeatError::HistoryMismatch
        ))
    ));
    let mut changed_mass = pair.clone();
    changed_mass.properties[0] = MassProperties::new(
        4.0 * RELATIVE_MASS,
        Vec3::ZERO,
        Vec3::new(1e-9, 1.2e-9, 1.4e-9),
    )
    .unwrap();
    assert!(matches!(
        attempt(&changed_mass.prepare(&source), saved_stamp),
        Err(RigidJkrTimedSnapHeatError::Mechanics(_))
    ));
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(heat.stamp(), saved_stamp);
    assert_eq!(pair.bodies, saved_bodies);
    assert!(attempt(&pair.prepare(&source), saved_stamp).is_ok());
}

#[test]
fn late_caloric_refusal_and_last_checkpoint_cancel_the_complete_native_candidate() {
    let source = source();
    let caloric = curves(&source, [1.0; 2]);
    let pair = NativePair::collision(&source);
    let prepared = pair.prepare(&source);
    let time = scales(&source).3;
    let thermal = thermal(&source, &caloric, 0.0, vec![], 0.0);
    let initial = source.detached_state(prepared.approach_m()).unwrap();
    let heat = binding(&source, &caloric, &thermal, initial);
    let saved_stamp = heat.stamp();
    let saved_bodies = pair.bodies;
    let saved_cells = thermal.cells().clone();
    let mut checkpoints = 0;
    let reference = heat
        .propose_timed_rigid_motion(
            &thermal,
            saved_stamp,
            &prepared,
            [Wrench::ZERO; 2],
            0.1 * time,
            0.0,
            native_controls(&source),
            &CancelGate::new_clock_free(),
            || {
                checkpoints += 1;
                Ok(())
            },
        )
        .unwrap();
    assert!(checkpoints > 4);
    let mut calls = 0;
    let rejected = heat.propose_timed_rigid_motion(
        &thermal,
        saved_stamp,
        &prepared,
        [Wrench::ZERO; 2],
        0.1 * time,
        0.0,
        native_controls(&source),
        &CancelGate::new_clock_free(),
        || {
            calls += 1;
            if calls == checkpoints {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        rejected,
        Err(RigidJkrTimedSnapHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(heat.stamp(), saved_stamp);
    assert_eq!(pair.bodies, saved_bodies);
    let retry = heat
        .propose_timed_rigid_motion(
            &thermal,
            saved_stamp,
            &prepared,
            [Wrench::ZERO; 2],
            0.1 * time,
            0.0,
            native_controls(&source),
            &CancelGate::new_clock_free(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(retry.report(), reference.report());
    assert_eq!(retry.thermal().cells(), reference.thermal().cells());

    let narrow = curves(&source, [1e-4, 1.0]);
    let narrow_thermal = self::thermal(&source, &narrow, 0.0, vec![], 0.0);
    let narrow_heat = binding(&source, &narrow, &narrow_thermal, initial);
    assert!(matches!(
        narrow_heat.propose_timed_rigid_motion(
            &narrow_thermal,
            narrow_heat.stamp(),
            &prepared,
            [Wrench::ZERO; 2],
            0.1 * time,
            0.0,
            native_controls(&source),
            &CancelGate::new_clock_free(),
            || Ok(()),
        ),
        Err(RigidJkrTimedSnapHeatError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::MaterialDomain
        )))
    ));
    assert_eq!(narrow_thermal.time_s(), 0.0);
    assert_eq!(narrow_heat.ledger().dissipated_work_j(), 0.0);
    assert_eq!(pair.bodies, saved_bodies);
}
