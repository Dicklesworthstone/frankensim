//! Powered contact against independent native work and finite-cell equations.
use super::*;
use fs_mbd::Wrench;

fn drive() -> [Wrench; 2] {
    [
        Wrench {
            force_world: Vec3::new(2.0, 0.0, 0.0),
            torque_body: Vec3::new(0.0, -0.01, 0.0),
        },
        Wrench {
            force_world: Vec3::new(-0.5, 0.0, 0.0),
            torque_body: Vec3::new(0.0, 0.005, 0.0),
        },
    ]
}

fn propose(
    f: &Fixture,
    t: &DepositionSimulation<'static>,
    external: [Wrench; 2],
) -> HeatedRigidHertzRollingCandidate<'static, [f64; 3]> {
    propose_rigid_hertz_rolling_driven_heat(
        &f.bodies,
        &f.properties,
        &RigidRollingHistory::new(),
        f.contact(),
        external,
        t,
        &route(f, t),
        bulk_binding(),
        t.time_s(),
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |t, face| source(f, t, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap()
}

// Independent planar force/torque and two caloric-node equations. This uses
// original material coefficients but no production mobility or work receipt.
fn reference_driven(f: &Fixture, external: [Wrench; 2]) -> Reference {
    let mass = f.properties[0].mass();
    let ia = 0.4 * mass * R * R;
    let ib = 0.13;
    let a = mass.recip() + 1.0 / 7.0 + R * R / ia + D * D / ib;
    let b = -R / ia + D / ib;
    let c = ia.recip() + ib.recip();
    let drive_omega = external[0].torque_body.y / ia - external[1].torque_body.y / ib;
    let drive_slip = external[0].force_world.x / mass
        - external[1].force_world.x / 7.0
        - R * external[0].torque_body.y / ia
        - D * external[1].torque_body.y / ib;
    let damping = f.tau * f.patch.normal_load_n() * R;
    let capacities = [mass * CP[0], 7.0 * CP[1]];
    let at = |face: f64| {
        let mu = super::super::MU[0] - super::super::SLOPE[0] * (face - 300.0);
        let force = -mu * f.patch.normal_load_n();
        let omega =
            (f.roll + 0.5 * DT * (b * force + drive_omega)) / (1.0 + 0.5 * DT * damping * c);
        let moment = -damping * omega;
        let slip = f.slip + 0.5 * DT * (a * force + b * moment + drive_slip);
        let surface_heat = -DT * force * slip;
        let bulk_heat = DT * damping * omega * omega;
        let temperatures = std::array::from_fn(|i| {
            (capacities[i] * 300.0 + DT * SIDE_G[i] * face + if i == 0 { bulk_heat } else { 0.0 })
                / (capacities[i] + DT * SIDE_G[i])
        });
        let residual = DT
            * (0..2)
                .map(|i| SIDE_G[i] * (face - temperatures[i]))
                .sum::<f64>()
            - surface_heat;
        (
            residual,
            Reference {
                face,
                temperatures,
                force,
                moment,
                surface_heat,
                bulk_heat,
            },
        )
    };
    let (mut lo, mut hi) = (280.0, 320.0);
    assert!(at(lo).0 < 0.0 && at(hi).0 > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if at(mid).0 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    at(0.5 * (lo + hi)).1
}

fn independent_external_work(
    f: &Fixture,
    after: &[RigidBodyState; 2],
    external: [Wrench; 2],
) -> f64 {
    (0..2)
        .map(|i| {
            let v = f.bodies[i]
                .center_of_mass_velocity_world(f.properties[i])
                .unwrap()
                .add(
                    after[i]
                        .center_of_mass_velocity_world(f.properties[i])
                        .unwrap(),
                );
            let w = f.properties[i]
                .angular_velocity_body_checked(f.bodies[i].angular_momentum_body())
                .unwrap()
                .add(
                    f.properties[i]
                        .angular_velocity_body_checked(after[i].angular_momentum_body())
                        .unwrap(),
                );
            0.5 * DT * (external[i].force_world.dot(v) + external[i].torque_body.dot(w))
        })
        .sum()
}

fn energy_audit(
    f: &Fixture,
    before: &DepositionSimulation<'_>,
    answer: &HeatedRigidHertzRollingCandidate<'_, [f64; 3]>,
    external: [Wrench; 2],
) {
    let work = independent_external_work(f, &answer.bodies_after, external);
    near(answer.external_mechanical_work_j, work, 2e-14);
    near(answer.mechanical.external_wrench_work_j, work, 2e-14);
    near(
        pair_energy(f, &answer.bodies_after) - pair_energy(f, &f.bodies)
            + energy_change(before, &answer.thermal)
            - answer.external_heat_j
            - work,
        0.0,
        3e-10,
    );
    near(
        energy_change(before, &answer.thermal) - answer.external_heat_j,
        answer.mechanical.surface_heat_j + answer.mechanical.bulk_heat_j,
        2e-10,
    );
}

#[test]
fn powered_surface_and_bulk_heat_match_independent_finite_temperature_solution() {
    let f = Fixture::new(TAU, SLIP, ROLL);
    let t = thermal(&f, [300.0; 2], false);
    let expected = reference_driven(&f, drive());
    let answer = propose(&f, &t, drive());
    energy_audit(&f, &t, &answer, drive());
    assert!(answer.external_mechanical_work_j > 0.0);
    assert!(answer.mechanical.surface_heat_j > 0.0 && answer.mechanical.bulk_heat_j > 0.0);
    near(answer.source[0], expected.face, 2e-7);
    for i in 0..2 {
        near(
            answer.thermal.temperature_k(IDS[i]).unwrap(),
            expected.temperatures[i],
            2e-7,
        );
    }
    near(answer.mechanical.hertz.force_n[0], expected.force, 2e-9);
    near(
        answer.mechanical.surface_heat_j,
        expected.surface_heat,
        2e-10,
    );
    near(answer.mechanical.bulk_heat_j, expected.bulk_heat, 2e-10);
    // Heating changes the original coefficient actually used by mechanics.
    assert!(answer.source[0] > 300.0);
    assert!(
        answer.mechanical.hertz.force_n[0].abs() < super::super::MU[0] * f.patch.normal_load_n()
    );
    assert_eq!(t.time_s(), 0.0);
}

#[test]
fn driven_mixed_sliding_spin_and_rolling_keep_common_acceleration_out_of_heat() {
    let mut f = Fixture::new(TAU, SLIP, ROLL);
    let mass = f.properties[0].mass();
    let inertia = 0.4 * mass * R * R;
    f.bodies[0] = RigidBodyState::new(
        f.bodies[0].pose(),
        Vec3::new(mass * (SLIP + R * ROLL), mass * 0.1, 0.0),
        Vec3::new(inertia * 0.2, inertia * ROLL, inertia * 0.3),
    )
    .unwrap();
    let t = thermal(&f, [300.0; 2], false);
    let mut external = drive();
    external[0].force_world.y = 0.3;
    external[0].torque_body = Vec3::new(0.001, -0.003, 0.001);
    let base = propose(&f, &t, external);
    energy_audit(&f, &t, &base, external);
    assert!(base.mechanical.hertz.force_n.iter().all(|v| v.abs() > 1e-5));
    assert!(base.mechanical.hertz.spin_moment_nm.abs() > 1e-6);
    assert!(base.mechanical.bulk_heat_j > 0.0);
    let common = Vec3::new(0.3, -0.2, 0.0);
    for i in 0..2 {
        external[i].force_world = external[i]
            .force_world
            .add(common.scale(f.properties[i].mass()));
    }
    let accelerated = propose(&f, &t, external);
    energy_audit(&f, &t, &accelerated, external);
    near(
        accelerated.mechanical.surface_heat_j,
        base.mechanical.surface_heat_j,
        2e-10,
    );
    near(
        accelerated.mechanical.bulk_heat_j,
        base.mechanical.bulk_heat_j,
        2e-10,
    );
    for i in 0..2 {
        near(
            accelerated.thermal.temperature_k(IDS[i]).unwrap(),
            base.thermal.temperature_k(IDS[i]).unwrap(),
            2e-7,
        );
    }
    assert!(
        (accelerated.external_mechanical_work_j - base.external_mechanical_work_j).abs() > 1e-5
    );
}

#[test]
fn driven_heat_preserves_exact_zero_drive_and_refuses_unresolved_normal_motion() {
    let f = Fixture::new(TAU, SLIP, ROLL);
    let t = thermal(&f, [300.0; 2], false);
    let history = RigidRollingHistory::new();
    let old = propose_rigid_hertz_rolling_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route(&f, &t),
        bulk_binding(),
        t.time_s(),
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |t, face| source(&f, t, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap();
    let zero = propose(&f, &t, [Wrench::ZERO; 2]);
    assert_eq!(old.bodies_after, zero.bodies_after);
    assert_eq!(old.mechanical, zero.mechanical);
    assert_eq!(old.rolling_history, zero.rolling_history);
    assert_eq!(old.thermal.cells(), zero.thermal.cells());
    assert_eq!(zero.external_mechanical_work_j, 0.0);
    let mut external = drive();
    external[0].force_world.z = 1.0;
    let saved = t.clone();
    let failed = propose_rigid_hertz_rolling_driven_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        external,
        &t,
        &route(&f, &t),
        bulk_binding(),
        t.time_s(),
        DT,
        controls(),
        &[0.0; 2],
        300.0,
        |t, face| source(&f, t, face, [0.0; 2]),
        &CancelGate::new(),
    );
    assert!(matches!(
        failed,
        Err(RigidHertzRollingHeatError::Mechanical(_))
    ));
    assert_eq!(t.time_s(), saved.time_s());
    assert_eq!(t.cells(), saved.cells());
    assert_eq!(t.links(), saved.links());
    assert_eq!(history, RigidRollingHistory::new());
}
