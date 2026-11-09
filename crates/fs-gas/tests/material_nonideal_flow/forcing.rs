//! Mechanical impulse/work checks independent of the real-gas caloric solver.

use super::*;
use fs_gas::body_force::{BodyForce, ForceControl};
use fs_gas::dynamics::nonideal::{RealGasWallContact, RealGasWallOperator};
use fs_gas::multicomponent::wall_exchange::{WallContact, WallControl, WallDrag, WallThermal};

fn force_controls() -> ForceControl {
    ForceControl { maximum_mach: 4.0 }
}

#[test]
fn exact_dense_gas_kick_preserves_thermal_state_and_original_energy_gauge() {
    let original = gas();
    let shifted = bind(
        common::model([1e22, -3e21], [2e18, -7e17], [200.0, 2000.0]),
        domain(),
    );
    for source in [original, shifted] {
        let initial = [
            state(&source, 5.0, -30.0, 350.0),
            state(&source, 20.0, 60.0, 500.0),
        ];
        let mut flow = uniform_mesh_field(source, initial.to_vec());
        let before = flow.primitives().to_vec();
        let acceleration = [800.0, -1200.0];
        let dt = 0.02;
        let exchange = flow
            .apply_body_force(
                Time::new(dt),
                BodyForce::Acceleration(&acceleration),
                force_controls(),
                || false,
            )
            .unwrap();
        let mut impulse = 0.0;
        let mut work = 0.0;
        for i in 0..2 {
            let expected_velocity = before[i].velocity + acceleration[i] * dt;
            let mass = initial[i].density * flow.mesh().volumes()[i];
            impulse += mass * acceleration[i] * dt;
            work += 0.5 * mass * (expected_velocity.powi(2) - before[i].velocity.powi(2));
            let after = flow.primitives()[i];
            near(after.velocity, expected_velocity, 2e-13);
            near(after.temperature, before[i].temperature, 2e-9);
            near(after.pressure, before[i].pressure, 2e-5);
            near(after.sound_speed, before[i].sound_speed, 2e-9);
            assert_eq!(flow.states()[i].density, initial[i].density);
        }
        near(exchange.impulse, impulse, 2e-13);
        near(exchange.work, work, 2e-11);
        near(exchange.energy_residual, 0.0, 2e-10);
        assert_eq!(flow.time().value(), 0.0);
    }
}

#[test]
fn nonuniform_volume_force_density_uses_volume_not_mass() {
    let source = gas();
    let initial = vec![
        state(&source, 5.0, -10.0, 360.0),
        state(&source, 20.0, 40.0, 500.0),
    ];
    let mesh = Mesh1d::new(vec![0.0, 0.2, 1.0], vec![0.07; 3], vec![0.014, 0.056]).unwrap();
    let mut flow = RealGasEuler1d::new_on_mesh(
        source,
        mesh,
        initial.clone(),
        RealGasBoundaries::Periodic,
        || false,
    )
    .unwrap();
    let force_density = 9000.0;
    let dt = 0.004;
    let report = flow
        .apply_body_force(
            Time::new(dt),
            BodyForce::UniformForceDensity(force_density),
            force_controls(),
            || false,
        )
        .unwrap();
    near(report.impulse, force_density * 0.07 * dt, 1e-13);
    for (i, q) in initial.iter().enumerate() {
        near(
            flow.states()[i].momentum - q.momentum,
            force_density * dt,
            1e-12,
        );
        near(
            flow.primitives()[i].velocity,
            q.momentum / q.density + force_density * dt / q.density,
            1e-12,
        );
    }
}

#[test]
fn forced_periodic_transport_closes_impulse_work_and_keeps_dense_gas_cold() {
    for acceleration in [-200_000.0, 200_000.0] {
        let source = gas();
        let q = state(&source, 20.0, -11.0, 430.0);
        let mut flow = uniform_mesh_field(source, vec![q; 5]);
        let mut impulse = 0.0;
        let mut work = 0.0;
        for _ in 0..6 {
            let report = flow
                .step_with_body_force(
                    Time::new(0.00005),
                    controls(),
                    BodyForce::UniformAcceleration(acceleration),
                    force_controls(),
                    || false,
                )
                .unwrap();
            impulse += report.source_change.impulse;
            work += report.source_change.work;
            near(
                report.source_change.impulse,
                report.before.mass * acceleration * report.elapsed.value(),
                2e-13,
            );
            near(report.conservation_residual.momentum, 0.0, 5e-12);
            near(report.conservation_residual.energy, 0.0, 5e-10);
            assert_eq!(report.boundary_change, Default::default());
            assert!(report.maximum_cfl <= controls().cfl * (1.0 + 64.0 * f64::EPSILON));
        }
        let expected = -11.0 + acceleration * flow.time().value();
        for p in flow.primitives() {
            near(p.velocity, expected, 5e-12);
            near(p.temperature, 430.0, 2e-9);
        }
        let total = flow.totals(|| false).unwrap();
        near(impulse, total.mass * (expected + 11.0), 2e-12);
        near(work, 0.5 * total.mass * (expected * expected - 121.0), 1e-9);
    }
}

#[test]
fn force_endpoint_cfl_retries_restart_both_kicks_from_accepted_state() {
    let source = gas();
    let q = state(&source, 20.0, 0.0, 430.0);
    let mut flow = uniform_mesh_field(source, vec![q; 3]);
    let acceleration = 2e7;
    let report = flow
        .step_with_body_force(
            Time::new(1.0),
            controls(),
            BodyForce::UniformAcceleration(acceleration),
            ForceControl {
                maximum_mach: 100.0,
            },
            || false,
        )
        .unwrap();
    assert!(
        report.halvings > 0,
        "unforced CFL cannot admit the accelerated endpoint"
    );
    near(
        report.source_change.impulse,
        report.before.mass * acceleration * report.elapsed.value(),
        1e-10,
    );
    for p in flow.primitives() {
        near(p.velocity, acceleration * report.elapsed.value(), 1e-10);
        near(p.temperature, 430.0, 3e-8);
    }
    assert!(report.maximum_cfl <= controls().cfl * (1.0 + 64.0 * f64::EPSILON));
    near(report.conservation_residual.energy, 0.0, 1e-8);
}

#[test]
fn zero_forcing_retains_the_original_wave_path_exactly() {
    let source = gas();
    let states = vec![
        state(&source, 10.0, 15.0, 400.0),
        state(&source, 15.0, -20.0, 450.0),
    ];
    let mut original = uniform_mesh_field(source.clone(), states.clone());
    let mut forced = uniform_mesh_field(source, states);
    let a = original
        .step(Time::new(1e-5), controls(), || false)
        .unwrap();
    let b = forced
        .step_with_body_force(
            Time::new(1e-5),
            controls(),
            BodyForce::UniformAcceleration(0.0),
            force_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(original.states(), forced.states());
    assert_eq!(original.primitives(), forced.primitives());
}

fn wall_case() -> (RealGasEuler1d, RealGasWallOperator<'static>) {
    let source = gas();
    let q = state(&source, 20.0, 50.0, 430.0);
    let flow = uniform_mesh_field(source, vec![q; 3]);
    let wall = RealGasWallContact {
        contact: WallContact {
            drag: WallDrag::MovingLinear {
                rate_per_s: 2000.0,
                wall_velocity_m_s: 20.0,
            },
            thermal: WallThermal::Lumped {
                temperature_k: 500.0,
                heat_capacity_j_per_k: 80.0,
                conductance_w_per_k: 6000.0,
            },
        },
        heat_storage: None,
    };
    let walls = RealGasWallOperator::new(&flow, vec![wall; 3], || false).unwrap();
    (flow, walls)
}

fn wall_controls() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-7,
        relative_energy_tolerance: 1e-10,
        maximum_iterations: 128,
    }
}

#[test]
fn real_gas_wall_drag_heat_and_external_body_work_have_distinct_balances() {
    let (mut flow, mut walls) = wall_case();
    let acceleration = 200_000.0;
    let report = flow
        .step_with_walls_and_body_force(
            Time::new(0.00004),
            controls(),
            &mut walls,
            wall_controls(),
            BodyForce::UniformAcceleration(acceleration),
            force_controls(),
            || false,
        )
        .unwrap();
    let thermal = walls.last_report().unwrap();
    let forced_velocity = 50.0 + acceleration * report.elapsed.value();
    let final_velocity = 20.0 + (forced_velocity - 20.0) * (-2000.0 * report.elapsed.value()).exp();
    for p in flow.primitives() {
        near(p.velocity, final_velocity, 2e-12);
    }
    near(
        report.source_change.work,
        0.5 * report.before.mass * (forced_velocity.powi(2) - 2500.0),
        1e-10,
    );
    near(
        report.wall_momentum_impulse,
        report.before.mass * (final_velocity - forced_velocity),
        1e-10,
    );
    near(report.wall_work, 20.0 * report.wall_momentum_impulse, 1e-9);
    near(report.conservation_residual.momentum, 0.0, 1e-10);
    near(report.conservation_residual.energy, 0.0, 3e-7);
    near(thermal.gas_wall_storage_residual_j, 0.0, 3e-7);
    near(
        report.after.energy - report.before.energy + thermal.wall_storage_change_j,
        report.boundary_change.energy
            + report.source_change.work
            + thermal.external_heat_to_gas_j
            + report.wall_work,
        3e-7,
    );
}

#[test]
fn invalid_mach_and_late_cancellation_leave_gas_cache_wall_inventory_and_clock_unchanged() {
    let (mut flow, mut walls) = wall_case();
    let states = flow.states().to_vec();
    let primitive = flow.primitives().to_vec();
    assert!(
        flow.apply_body_force(
            Time::new(1.0),
            BodyForce::UniformAcceleration(1e8),
            force_controls(),
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.states(), states);
    assert_eq!(flow.primitives(), primitive);
    let contact = walls.contacts().to_vec();
    let mut polls = 0;
    {
        let (mut probe, mut probe_walls) = wall_case();
        probe
            .step_with_walls_and_body_force(
                Time::new(0.00004),
                controls(),
                &mut probe_walls,
                wall_controls(),
                BodyForce::UniformAcceleration(200_000.0),
                force_controls(),
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
    }
    let mut seen = 0;
    assert_eq!(
        flow.step_with_walls_and_body_force(
            Time::new(0.00004),
            controls(),
            &mut walls,
            wall_controls(),
            BodyForce::UniformAcceleration(200_000.0),
            force_controls(),
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(RealGasFlowError::Cancelled)
    );
    assert_eq!(flow.states(), states);
    assert_eq!(flow.primitives(), primitive);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(walls.contacts(), contact);
    assert!(walls.last_report().is_none());
}

#[test]
fn dense_gas_ramp_uses_absolute_time_and_admits_the_hidden_interior_velocity() {
    let source = gas();
    let q = state(&source, 20.0, -11.0, 430.0);
    let mut flow = uniform_mesh_field(source, vec![q; 5]);
    let initial = flow.totals(|| false).unwrap();
    let (reference_time, acceleration, jerk) = (-0.001, -200_000.0, 400_000_000.0);
    let force = BodyForce::UniformAccelerationRamp {
        reference_time,
        acceleration,
        jerk,
    };
    let (mut impulse, mut work) = (0.0, 0.0);
    for _ in 0..12 {
        let start = flow.time().value();
        let report = flow
            .step_with_body_force(
                Time::new(0.000025), controls(), force, force_controls(), || false,
            )
            .unwrap();
        assert_eq!(report.elapsed, Time::new(0.000025));
        let h = report.elapsed.value();
        let delta = h * (acceleration + jerk * (start - reference_time) + 0.5 * jerk * h);
        near(report.source_change.impulse, initial.mass * delta, 2e-11);
        near(report.conservation_residual.momentum, 0.0, 2e-11);
        near(report.conservation_residual.energy, 0.0, 5e-9);
        impulse += report.source_change.impulse;
        work += report.source_change.work;
    }
    let time = flow.time().value();
    let velocity = -11.0 + (acceleration - jerk * reference_time) * time
        + 0.5 * jerk * time * time;
    near(impulse, initial.mass * (velocity + 11.0), 5e-10);
    near(work, 0.5 * initial.mass * (velocity * velocity - 121.0), 5e-8);
    near(flow.totals(|| false).unwrap().energy - initial.energy, work, 5e-8);
    for p in flow.primitives() {
        near(p.velocity, velocity, 5e-10);
        near(p.temperature, 430.0, 5e-8);
    }

    let source = gas();
    let q = state(&source, 20.0, 0.0, 430.0);
    let mut reversal = uniform_mesh_field(source, vec![q; 3]);
    let before = reversal.states().to_vec();
    let primitives = reversal.primitives().to_vec();
    let limit = ForceControl {
        maximum_mach: 25.0 / primitives[0].sound_speed,
    };
    // This source returns to zero velocity, but its interior 50 m/s state
    // must pass the original PR inversion and then the authored Mach gate.
    assert!(
        reversal.apply_body_force(
            Time::new(0.02),
            BodyForce::UniformAccelerationRamp {
                reference_time: 0.0,
                acceleration: 10_000.0,
                jerk: -1_000_000.0,
            },
            limit,
            || false,
        )
        .is_err()
    );
    assert_eq!(reversal.states(), before);
    assert_eq!(reversal.primitives(), primitives);
    assert_eq!(reversal.time(), Time::new(0.0));
}
