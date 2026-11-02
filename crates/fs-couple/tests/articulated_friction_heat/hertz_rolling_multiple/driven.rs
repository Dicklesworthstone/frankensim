//! Powered shared-body contact islands retain distinct physical heat owners.
use super::*;
use fs_mbd::Wrench;

fn powered<'m>(
    f: &Fixture,
    thermal: &DepositionSimulation<'m>,
    histories: &[RigidRollingHistory; 2],
    external: &[Wrench],
    outside: [f64; 3],
) -> HeatedRigidHertzRollingIslandCandidate<'m, Vec<f64>> {
    let routes = routes(f, thermal);
    propose_rigid_hertz_rolling_island_driven_heat(
        &f.bodies,
        &f.properties,
        &f.heat_contacts(&routes, histories),
        external,
        thermal,
        thermal.time_s(),
        DT,
        controls(),
        &[0.0; 3],
        |state, ids, faces| source(f, state, ids, faces, outside),
        &CancelGate::new(),
    )
    .unwrap()
}

fn drive() -> [Wrench; 3] {
    [
        Wrench {
            force_world: Vec3::new(8.0, 0.0, 0.0),
            torque_body: Vec3::new(0.0, 0.04, 0.0),
        },
        Wrench::ZERO,
        Wrench::ZERO,
    ]
}

// The existing independent reference solves scalar body force/torque balances
// and the three finite capacities. In a midpoint formula, external COM force
// and torque add half their endpoint velocity changes to the free midpoint.
// Both original contacts remain positively sliding; this reference algebra
// does not choose or replace the production owner's friction branch.
fn powered_reference(f: &Fixture, old: [f64; 3], outside: [f64; 3]) -> (Reference, f64) {
    let external = drive()[0];
    let mass = f.properties[0].mass();
    let inertia = 0.4 * mass * R * R;
    let mut midpoint = f.clone();
    midpoint.speed += 0.5 * DT * external.force_world.x / mass;
    midpoint.omega += 0.5 * DT * external.torque_body.y / inertia;
    let mut result = reference(&midpoint, old, outside, DT);
    result.momentum[0] += 0.5 * DT * external.force_world.x;
    result.angular[0] += 0.5 * DT * external.torque_body.y;
    let work = 0.5
        * DT
        * (external.force_world.x * (f.speed + result.momentum[0] / mass)
            + external.torque_body.y * (f.omega + result.angular[0] / inertia));
    (result, work)
}

fn actual_drive_work(f: &Fixture, next: &[RigidBodyState], external: &[Wrench]) -> f64 {
    (0..3)
        .map(|i| {
            let a = f.bodies[i];
            let b = next[i];
            let p = f.properties[i];
            let velocity = a
                .linear_momentum_world()
                .add(b.linear_momentum_world())
                .scale(0.5 / p.mass());
            let omega = p
                .angular_velocity_body_checked(
                    a.angular_momentum_body()
                        .add(b.angular_momentum_body())
                        .scale(0.5),
                )
                .unwrap();
            DT * (external[i].force_world.dot(velocity) + external[i].torque_body.dot(omega))
        })
        .sum()
}

#[test]
fn powered_two_patch_thermal_feedback_matches_independent_shared_sphere_equations() {
    let f = Fixture::new(TAU, 0.2, 0.05);
    let old = [299.0, 301.0, 303.0];
    let outside = [0.01, -0.02, 0.015];
    let t = thermal(&f, old, false);
    let histories = [RigidRollingHistory::new(), RigidRollingHistory::new()];
    let result = powered(&f, &t, &histories, &drive(), outside);
    let (expected, work) = powered_reference(&f, old, outside);
    assert_eq!(result.mechanical.body_indices, vec![0, 1, 2]);
    assert_eq!(result.bulk_heat.len(), 1);
    assert_eq!(result.bulk_heat[0].binding, bulk_binding());
    near(
        result.bulk_heat[0].dissipated_energy_j,
        expected.bulk.iter().sum(),
        8e-11,
    );
    for i in 0..2 {
        let contact = &result.mechanical.report.contacts[i];
        assert_eq!(contact.patch_id, PATCHES[i]);
        near(contact.hertz.force_n[0], expected.force[i], 4e-8);
        near(
            contact.rolling_angular_impulse_world_nms.y,
            DT * expected.moment[i],
            2e-11,
        );
        near(contact.surface_heat_j, expected.surface[i], 2e-10);
        near(contact.bulk_heat_j, expected.bulk[i], 8e-11);
        near(contact.rolling_temperature_k, expected.temperature[0], 2e-7);
        near(
            result.mechanical.rolling_histories[i]
                .history
                .state()
                .cumulative_heat_j(),
            expected.bulk[i],
            8e-11,
        );
        assert_eq!(
            result.mechanical.rolling_histories[i]
                .history
                .state()
                .accepted_steps(),
            1
        );
    }
    for i in 0..3 {
        near(
            result.thermal.temperature_k(IDS[i]).unwrap(),
            expected.temperature[i],
            2e-7,
        );
        near(
            result.mechanical.bodies_after[i].linear_momentum_world().x,
            expected.momentum[i],
            5e-10,
        );
        near(
            result.mechanical.bodies_after[i].angular_momentum_body().y,
            expected.angular[i],
            5e-11,
        );
    }
    for (j, i) in [1, 0].into_iter().enumerate() {
        assert_eq!(result.friction_sources[j].patch_id, PATCHES[i]);
        near(result.source[j], expected.face[i], 2e-7);
        let radius = (f.patches[i].contact_area_m2() / core::f64::consts::PI).sqrt();
        assert_eq!(
            result.friction_sources[j].friction_source,
            f.sources[i]
                .hertz_at(result.source[j], None, f.patches[i].normal_load_n(), radius)
                .unwrap()
        );
    }
    let dk: f64 = (0..3)
        .map(|i| {
            result.mechanical.bodies_after[i]
                .kinetic_energy(f.properties[i])
                .unwrap()
                - f.bodies[i].kinetic_energy(f.properties[i]).unwrap()
        })
        .sum();
    let dh: f64 = t
        .cells()
        .iter()
        .map(|(id, cell)| {
            cell.mass_kg * (result.thermal.cells()[id].specific_h_j_kg - cell.specific_h_j_kg)
        })
        .sum();
    near(result.external_mechanical_work_j, work, 2e-10);
    near(
        actual_drive_work(&f, &result.mechanical.bodies_after, &drive()),
        work,
        2e-10,
    );
    near(
        result.external_heat_j,
        DT * outside.iter().sum::<f64>(),
        1e-14,
    );
    near(dk + dh - result.external_heat_j - work, 0.0, 7e-10);
    near(
        dh - result.external_heat_j
            - expected.surface.iter().sum::<f64>()
            - expected.bulk.iter().sum::<f64>(),
        0.0,
        7e-10,
    );
    near(result.combined_energy_residual_j, 0.0, 7e-10);
    let unforced = reference(&f, old, outside, DT);
    assert!(
        (expected.surface.iter().sum::<f64>() - unforced.surface.iter().sum::<f64>()).abs() > 1e-6
    );
    assert!(expected.bulk.iter().all(|q| *q > 1e-9));
}

#[test]
fn zero_drive_equivalence_and_shared_common_acceleration_preserve_each_heat_owner() {
    let f = Fixture::new(TAU, 0.2, 0.05);
    let t = thermal(&f, [300.0; 3], false);
    let routes = routes(&f, &t);
    let histories = [RigidRollingHistory::new(), RigidRollingHistory::new()];
    let plain = propose_rigid_hertz_rolling_island_heat(
        &f.bodies,
        &f.properties,
        &f.heat_contacts(&routes, &histories),
        &t,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |state, ids, faces| source(&f, state, ids, faces, [0.0; 3]),
        &CancelGate::new(),
    )
    .unwrap();
    let zero = powered(&f, &t, &histories, &[Wrench::ZERO; 3], [0.0; 3]);
    assert_eq!(plain.mechanical, zero.mechanical);
    assert_eq!(plain.thermal.cells(), zero.thermal.cells());
    assert_eq!(plain.source, zero.source);
    assert_eq!(zero.external_mechanical_work_j, 0.0);
    let base = powered(&f, &t, &histories, &drive(), [0.0; 3]);
    let acceleration = 2.0;
    let mut common = drive();
    for i in 0..3 {
        common[i].force_world.x += f.properties[i].mass() * acceleration;
    }
    let boosted = powered(&f, &t, &histories, &common, [0.0; 3]);
    for i in 0..2 {
        near(
            boosted.mechanical.report.contacts[i].surface_heat_j,
            base.mechanical.report.contacts[i].surface_heat_j,
            2e-10,
        );
        near(
            boosted.mechanical.report.contacts[i].bulk_heat_j,
            base.mechanical.report.contacts[i].bulk_heat_j,
            8e-11,
        );
    }
    for i in 0..3 {
        near(
            boosted.thermal.temperature_k(IDS[i]).unwrap(),
            base.thermal.temperature_k(IDS[i]).unwrap(),
            2e-7,
        );
        near(
            boosted.mechanical.bodies_after[i].linear_momentum_world().x,
            base.mechanical.bodies_after[i].linear_momentum_world().x
                + f.properties[i].mass() * acceleration * DT,
            5e-10,
        );
    }
    let mass: f64 = f.properties.iter().map(|p| p.mass()).sum();
    let initial_momentum: f64 = f.bodies.iter().map(|b| b.linear_momentum_world().x).sum();
    let expected_extra = acceleration * DT * initial_momentum
        + acceleration * DT * DT * drive()[0].force_world.x
        + 0.5 * mass * (acceleration * DT).powi(2);
    near(
        boosted.external_mechanical_work_j - base.external_mechanical_work_j,
        expected_extra,
        3e-10,
    );
    near(
        boosted.external_mechanical_work_j,
        actual_drive_work(&f, &boosted.mechanical.bodies_after, &common),
        3e-10,
    );
    assert!(expected_extra > 0.01);
}

#[test]
fn loaded_powered_network_source_failure_and_late_cancellation_preserve_exact_retry() {
    let mut f = Fixture::new(TAU, 0.2, 0.05);
    let t0 = thermal(&f, [300.0; 3], false);
    let fresh = [RigidRollingHistory::new(), RigidRollingHistory::new()];
    let first = powered(&f, &t0, &fresh, &drive(), [0.0; 3]);
    f.bodies.copy_from_slice(&first.mechanical.bodies_after);
    let histories = std::array::from_fn(|i| first.mechanical.rolling_histories[i].history.clone());
    let t = first.thermal;
    let routes = routes(&f, &t);
    let accepted = (f.bodies, histories.clone(), t.cells().clone());
    let baseline = powered(&f, &t, &histories, &drive(), [0.0; 3]);
    for cancel in [false, true] {
        let gate = CancelGate::new();
        let calls = Cell::new(0);
        let refused = propose_rigid_hertz_rolling_island_driven_heat(
            &f.bodies,
            &f.properties,
            &f.heat_contacts(&routes, &histories),
            &drive(),
            &t,
            t.time_s(),
            DT,
            controls(),
            &[0.0; 3],
            |state, ids, faces| {
                calls.set(calls.get() + 1);
                if calls.get() == 3 {
                    if cancel {
                        gate.request();
                    } else {
                        return Err("late powered network source");
                    }
                }
                source(&f, state, ids, faces, [0.0; 3])
            },
            &gate,
        );
        assert!(refused.is_err());
        assert_eq!(calls.get(), 3);
        assert_eq!(f.bodies, accepted.0);
        assert_eq!(histories, accepted.1);
        assert_eq!(t.cells(), &accepted.2);
        assert_eq!(t.time_s(), DT);
    }
    let retry = powered(&f, &t, &histories, &drive(), [0.0; 3]);
    assert_eq!(baseline.mechanical, retry.mechanical);
    assert_eq!(baseline.thermal.cells(), retry.thermal.cells());
    assert_eq!(baseline.source, retry.source);
    assert_eq!(
        baseline.external_mechanical_work_j,
        retry.external_mechanical_work_j
    );
}
