//! Held external work, actual reversal, original sources and finite heat.
use super::*;
use fs_mbd::Wrench;
use fs_mbd::tangential_friction::hertz_rolling::collinear::RigidHertzRollingCollinearMode;

fn driven<'m>(
    f: &Fixture,
    thermal: &DepositionSimulation<'m>,
    history: &RigidHertzRollingCollinearHistory,
    external: [Wrench; 2],
    dt: f64,
) -> HeatedRigidHertzRollingCollinearCandidate<'m, [f64; 3]> {
    let route = route(f, thermal);
    let guess =
        0.5 * (thermal.temperature_k(IDS[0]).unwrap() + thermal.temperature_k(IDS[1]).unwrap());
    propose_rigid_hertz_rolling_collinear_driven_heat(
        &f.bodies,
        &f.properties,
        history,
        f.contact(),
        external,
        thermal,
        &route,
        bulk_binding(),
        thermal.time_s(),
        dt,
        event_controls(),
        &[0.0; 2],
        guess,
        |state, face| source(f, state, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap()
}

fn contact_drive(f: &Fixture, shear: f64, common_acceleration: f64) -> [Wrench; 2] {
    [
        Wrench {
            force_world: Vec3::new(
                shear + f.properties[0].mass() * common_acceleration,
                0.0,
                0.0,
            ),
            torque_body: Vec3::new(0.0, -R * shear, 0.0),
        },
        Wrench {
            force_world: Vec3::new(
                -shear + f.properties[1].mass() * common_acceleration,
                0.0,
                0.0,
            ),
            torque_body: Vec3::new(0.0, -D * shear, 0.0),
        },
    ]
}

#[test]
fn zero_drive_retains_exact_unforced_stopping_history_and_both_heat_destinations() {
    let mut f = Fixture::new(0.0005, 0.2, 0.05);
    f.source = original_source(0.8, 0.4, 0.01);
    let t = thermal(&f, [300.0; 2], false);
    let route = route(&f, &t);
    let history = RigidHertzRollingCollinearHistory::default();
    let plain = propose_rigid_hertz_rolling_collinear_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        0.0,
        5.0,
        event_controls(),
        &[0.0; 2],
        300.0,
        |state, face| source(&f, state, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap();
    let forced = driven(&f, &t, &history, [Wrench::ZERO; 2], 5.0);
    assert_eq!(
        plain.mechanical.bodies_after,
        forced.mechanical.bodies_after
    );
    assert_eq!(plain.mechanical.history, forced.mechanical.history);
    assert_eq!(plain.mechanical.report, forced.mechanical.report);
    assert_eq!(plain.thermal.cells(), forced.thermal.cells());
    assert_eq!(plain.source, forced.source);
    assert_eq!(forced.external_mechanical_work_j, 0.0);
    assert!(forced.mechanical.report.surface_heat_j > 0.0);
    assert!(forced.mechanical.report.bulk_heat_j > 0.0);
}

#[derive(Clone, Copy)]
struct Reversal {
    stop: f64,
    slip: f64,
    surface: f64,
    drive_work: f64,
}

// Independent two-branch constant-force integration. Applied force is at A's
// COM; Hertz shear acts at the two original arms. No native report is an input.
fn reversal(f: &Fixture, mu: f64, drive: f64, dt: f64) -> Reversal {
    let ma = f.properties[0].mass();
    let inertia = 0.4 * ma * R * R;
    let a = ma.recip() + 1.0 / 7.0 + R * R / inertia + D * D / 0.13;
    let traction = mu * f.patch.normal_load_n();
    let first_acceleration = drive / ma - a * traction;
    let second_acceleration = drive / ma + a * traction;
    let stop = -f.slip / first_acceleration;
    assert!(stop > 0.0 && stop < dt && second_acceleration < 0.0);
    let remaining = dt - stop;
    let path = 0.5 * f.slip * stop - 0.5 * second_acceleration * remaining.powi(2);
    let va = f.slip + R * f.roll;
    let aa = (drive - traction) / ma;
    let drive_displacement = va * stop
        + 0.5 * aa * stop * stop
        + (va + aa * stop) * remaining
        + 0.5 * (drive + traction) / ma * remaining.powi(2);
    Reversal {
        stop,
        slip: second_acceleration * remaining,
        surface: traction * path,
        drive_work: drive * drive_displacement,
    }
}

fn reversal_thermal(f: &Fixture, drive: f64, dt: f64, slope: f64) -> (f64, [f64; 2], Reversal) {
    let capacity = [
        f.properties[0].mass() * CP[0],
        f.properties[1].mass() * CP[1],
    ];
    let at = |face: f64| {
        let mechanical = reversal(f, 0.4 - slope * (face - 300.0), drive, dt);
        let temperatures = std::array::from_fn(|i| {
            (capacity[i] * 300.0 + dt * SIDE_G[i] * face) / (capacity[i] + dt * SIDE_G[i])
        });
        let residual = dt
            * (0..2)
                .map(|i| SIDE_G[i] * (face - temperatures[i]))
                .sum::<f64>()
            - mechanical.surface;
        (residual, temperatures, mechanical)
    };
    let (mut lo, mut hi) = (300.0, 320.0);
    assert!(at(lo).0 < 0.0 && at(hi).0 > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if at(mid).0 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let face = 0.5 * (lo + hi);
    let (_, temperatures, motion) = at(face);
    (face, temperatures, motion)
}

#[test]
fn original_thermal_weakening_changes_actual_forced_reversal_and_closes_external_work() {
    let mut values = Vec::new();
    for slope in [0.0, 0.01] {
        let mut f = Fixture::new(0.0, 0.2, 0.05);
        f.source = original_source(0.8, 0.4, slope);
        let t = thermal(&f, [300.0; 2], false);
        let drive = -60.0;
        let dt = 0.04;
        let external = [
            Wrench {
                force_world: Vec3::new(drive, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            },
            Wrench::ZERO,
        ];
        let result = driven(
            &f,
            &t,
            &RigidHertzRollingCollinearHistory::default(),
            external,
            dt,
        );
        let (face, temperatures, motion) = reversal_thermal(&f, drive, dt, slope);
        let phases = &result.mechanical.report.phases;
        assert_eq!(phases.len(), 2);
        assert_eq!(
            phases[0].mode,
            RigidHertzRollingCollinearMode::SlidingPositive
        );
        assert!(phases[0].reached_stop);
        assert_eq!(
            phases[1].mode,
            RigidHertzRollingCollinearMode::SlidingNegative
        );
        near(phases[0].duration_s, motion.stop, 5e-11);
        near(phases[1].shear_velocity_after_m_s, motion.slip, 2e-10);
        near(
            result.mechanical.report.surface_heat_j,
            motion.surface,
            2e-10,
        );
        assert_eq!(result.mechanical.report.bulk_heat_j, 0.0);
        near(result.external_mechanical_work_j, motion.drive_work, 2e-10);
        near(result.source[0], face, 2e-7);
        for i in 0..2 {
            near(
                result.thermal.temperature_k(IDS[i]).unwrap(),
                temperatures[i],
                2e-7,
            );
        }
        near(energy_change(&t, &result.thermal), motion.surface, 3e-10);
        near(
            pair_energy(&f, &result.mechanical.bodies_after) - pair_energy(&f, &f.bodies)
                + energy_change(&t, &result.thermal)
                - motion.drive_work,
            0.0,
            3e-10,
        );
        near(result.combined_energy_residual_j, 0.0, 3e-10);
        near(result.heat_transfer_residual_j, 0.0, 2e-10);
        values.push((motion.stop, motion.surface));
    }
    assert!(values[1].0 > values[0].0);
    assert!((values[1].1 - values[0].1).abs() > 1e-4);
}

#[test]
fn subthreshold_contact_drive_and_common_acceleration_supply_kinetic_energy_without_heat() {
    let mut f = Fixture::new(0.0005, 0.0, 0.0);
    f.source = original_source(0.8, 0.4, 0.0);
    let t = thermal(&f, [300.0; 2], false);
    let shear = 0.6 * f.patch.normal_load_n();
    let acceleration = 1.0;
    let dt = 0.2;
    let result = driven(
        &f,
        &t,
        &RigidHertzRollingCollinearHistory::default(),
        contact_drive(&f, shear, acceleration),
        dt,
    );
    assert!(result.mechanical.history.is_stopped());
    assert!(
        result
            .mechanical
            .report
            .phases
            .iter()
            .all(|phase| phase.mode == RigidHertzRollingCollinearMode::Sticking)
    );
    near(result.mechanical.report.surface_heat_j, 0.0, 1e-28);
    near(result.mechanical.report.bulk_heat_j, 0.0, 1e-28);
    let supplied =
        0.5 * (f.properties[0].mass() + f.properties[1].mass()) * (acceleration * dt).powi(2);
    near(result.external_mechanical_work_j, supplied, 3e-10);
    near(
        pair_energy(&f, &result.mechanical.bodies_after),
        supplied,
        3e-10,
    );
    near(energy_change(&t, &result.thermal), 0.0, 2e-10);
    for i in 0..2 {
        near(
            result.mechanical.bodies_after[i].linear_momentum_world().x,
            f.properties[i].mass() * acceleration * dt,
            2e-10,
        );
    }
    near(result.combined_energy_residual_j, 0.0, 3e-10);
}

#[test]
fn driven_late_source_and_cancellation_leave_original_inputs_and_exact_retry() {
    let mut f = Fixture::new(0.0, 0.2, 0.05);
    f.source = original_source(0.8, 0.4, 0.01);
    let t = thermal(&f, [300.0; 2], false);
    let route = route(&f, &t);
    let history = RigidHertzRollingCollinearHistory::default();
    let original = (f.bodies, t.cells().clone(), history.clone());
    let external = [
        Wrench {
            force_world: Vec3::new(-60.0, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ];
    let baseline = driven(&f, &t, &history, external, 0.04);
    for cancel in [false, true] {
        let gate = CancelGate::new();
        let calls = Cell::new(0);
        let failed = propose_rigid_hertz_rolling_collinear_driven_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(),
            external,
            &t,
            &route,
            bulk_binding(),
            0.0,
            0.04,
            event_controls(),
            &[0.0; 2],
            300.0,
            |state, face| {
                calls.set(calls.get() + 1);
                if calls.get() == 3 {
                    if cancel {
                        gate.request();
                    } else {
                        return Err("late driven source");
                    }
                }
                source(&f, state, face, [0.0; 2])
            },
            &gate,
        );
        assert!(failed.is_err());
        assert_eq!(calls.get(), 3);
        assert_eq!(f.bodies, original.0);
        assert_eq!(t.cells(), &original.1);
        assert_eq!(history, original.2);
        assert_eq!(t.time_s(), 0.0);
    }
    let retry = driven(&f, &t, &history, external, 0.04);
    assert_eq!(
        baseline.mechanical.bodies_after,
        retry.mechanical.bodies_after
    );
    assert_eq!(baseline.mechanical.history, retry.mechanical.history);
    assert_eq!(baseline.thermal.cells(), retry.thermal.cells());
    assert_eq!(
        baseline.external_mechanical_work_j,
        retry.external_mechanical_work_j
    );
}
