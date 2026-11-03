//! G0/G1/G3/G4: body-force work, closed-wall impulse, and atomic forced steps.

use fs_duct::superfluid::{
    BodyForce, ChannelError, ComponentForceDensity, HeatStepError, LinearizationLimits,
    MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState};
use std::cell::Cell;

fn material(alpha: f64) -> TwoFluidState {
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: alpha,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap()
}

fn channel(n: usize, alpha: f64) -> TwoFluidChannel {
    TwoFluidChannel::try_new(
        material(alpha),
        1.0,
        0.01,
        n,
        LinearizationLimits {
            max_pressure_perturbation_pa: 1.0,
            max_temperature_perturbation_k: 0.02,
            max_component_speed_m_s: 0.1,
        },
    )
    .unwrap()
}

fn force(normal: f64, superfluid: f64) -> ComponentForceDensity {
    ComponentForceDensity {
        normal_n_m3: normal,
        superfluid_n_m3: superfluid,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn momentum(duct: &TwoFluidChannel) -> f64 {
    // Full-domain trapezoidal quadrature, independently reconstructed from
    // component velocities. Closed endpoints have zero half-cell momentum.
    let volume = duct.area_m2() * duct.length_m() / duct.cells() as f64;
    (0..=duct.cells())
        .map(|i| {
            let face = duct.face(i).unwrap();
            let weight = if i == 0 || i == duct.cells() {
                0.5
            } else {
                1.0
            };
            weight
                * volume
                * (duct.material().normal_density_kg_m3() * face.normal_velocity_m_s
                    + duct.material().superfluid_density_kg_m3() * face.superfluid_velocity_m_s)
        })
        .sum()
}

#[test]
fn finite_expansion_hydrostatic_balance_includes_both_wall_half_cells() {
    for n in [2, 23] {
        for alpha in [-0.3, 0.3] {
            // The force half-kicks are algebraic auxiliaries: a balanced
            // stationary flow must satisfy even a tight physical speed bound.
            let mut duct = TwoFluidChannel::try_new(
                material(alpha),
                1.0,
                0.01,
                n,
                LinearizationLimits {
                    max_pressure_perturbation_pa: 1.0,
                    max_temperature_perturbation_k: 0.02,
                    max_component_speed_m_s: 1e-6,
                },
            )
            .unwrap();
            let acceleration = 0.02;
            let left_pressure = 0.01;
            let gradient = duct.material().density_kg_m3() * acceleration;
            let dx = duct.length_m() / n as f64;
            let pressure: Vec<_> = (0..n)
                .map(|i| left_pressure + gradient * (i as f64 + 0.5) * dx)
                .collect();
            duct.initialize(
                &pressure,
                &vec![0.0; n],
                &vec![0.0; n + 1],
                &vec![0.0; n + 1],
                || false,
            )
            .unwrap();
            let initial_energy = duct.wave_energy_j().unwrap();
            let dt = 0.4 * duct.maximum_step_s();
            for _ in 0..3 {
                let report = duct
                    .step_forced_closed(
                        dt,
                        BodyForce::Acceleration(acceleration),
                        MutualFriction::Linear { gamma_kg_m3_s: 3.0 },
                        || false,
                    )
                    .unwrap();
                let body = gradient * duct.area_m2() * duct.length_m() * dt;
                let left = duct.area_m2() * left_pressure * dt;
                let right = -duct.area_m2() * (left_pressure + gradient * duct.length_m()) * dt;
                close(report.body_impulse_n_s, body, body.abs());
                close(report.wall_impulse_n_s[0], left, body.abs());
                close(report.wall_impulse_n_s[1], right, body.abs());
                close(report.momentum_after_kg_m_s, 0.0, body.abs());
                close(momentum(&duct), 0.0, body.abs());
                close(report.momentum_defect_n_s, 0.0, body.abs());
                close(report.channel.body_force_work_j, 0.0, initial_energy);
                close(report.channel.mutual_friction_heat_j, 0.0, initial_energy);
                close(
                    duct.wave_energy_j().unwrap(),
                    initial_energy,
                    initial_energy,
                );
                for (i, expected) in pressure.iter().enumerate() {
                    let cell = duct.cell(i).unwrap();
                    close(cell.pressure_perturbation_pa, *expected, 0.1);
                    close(cell.temperature_perturbation_k, 0.0, 0.001);
                }
                for i in 0..=n {
                    let face = duct.face(i).unwrap();
                    close(face.normal_velocity_m_s, 0.0, acceleration * dt);
                    close(face.superfluid_velocity_m_s, 0.0, acceleration * dt);
                }
            }
        }
    }
}

#[test]
fn endpoint_body_forces_are_absorbed_by_stationary_walls_without_work() {
    let n = 8;
    let mut duct = channel(n, 0.3);
    let mut forces: Vec<_> = (0..=n).map(|_| force(0.0, 0.0)).collect();
    forces[0] = force(0.4, -0.1);
    forces[n] = force(-0.2, -0.6);
    let dt = 0.5 * duct.maximum_step_s();
    let half_volume = duct.area_m2() * duct.length_m() / (2.0 * n as f64);
    let left_body = half_volume * dt * (forces[0].normal_n_m3 + forces[0].superfluid_n_m3);
    let right_body = half_volume * dt * (forces[n].normal_n_m3 + forces[n].superfluid_n_m3);
    let report = duct
        .step_forced_closed(
            dt,
            BodyForce::Faces(&forces),
            MutualFriction::Linear {
                gamma_kg_m3_s: 10.0,
            },
            || false,
        )
        .unwrap();
    let scale = left_body.abs() + right_body.abs();
    close(report.body_impulse_n_s, left_body + right_body, scale);
    close(report.wall_impulse_n_s[0], -left_body, scale);
    close(report.wall_impulse_n_s[1], -right_body, scale);
    close(report.momentum_defect_n_s, 0.0, scale);
    assert_eq!(report.channel.boundary_pressure_impulse_n_s, [0.0; 2]);
    assert_eq!(report.momentum_before_kg_m_s, 0.0);
    assert_eq!(report.momentum_after_kg_m_s, 0.0);
    assert_eq!(report.channel.body_force_work_j, 0.0);
    assert_eq!(report.channel.mutual_friction_heat_j, 0.0);
    assert_eq!(duct.wave_energy_j().unwrap(), 0.0);
    for i in 0..=n {
        let face = duct.face(i).unwrap();
        assert_eq!(face.normal_velocity_m_s, 0.0);
        assert_eq!(face.superfluid_velocity_m_s, 0.0);
    }
}

#[test]
fn opposite_component_forces_drive_heat_without_creating_total_momentum() {
    let mut duct = channel(40, 0.0);
    let dt = 0.6 * duct.maximum_step_s();
    let mut supplied_work = 0.0;
    let mut heat = 0.0;
    for _ in 0..80 {
        let report = duct
            .step_forced_closed(
                dt,
                BodyForce::Uniform(force(0.004, -0.004)),
                MutualFriction::Linear {
                    gamma_kg_m3_s: 12.0,
                },
                || false,
            )
            .unwrap();
        supplied_work += report.channel.body_force_work_j;
        heat += report.channel.mutual_friction_heat_j;
        assert_eq!(report.body_impulse_n_s, 0.0);
        close(report.momentum_before_kg_m_s, 0.0, 1e-8);
        close(report.momentum_after_kg_m_s, 0.0, 1e-8);
        close(report.wall_impulse_n_s.iter().sum(), 0.0, 1e-8);
        let cell_heat: f64 = duct.last_mutual_friction_heat_j().iter().sum();
        close(
            cell_heat,
            report.channel.mutual_friction_heat_j,
            supplied_work,
        );
    }
    assert!(supplied_work > 0.0);
    assert!(heat > 0.01 * supplied_work);
    close(
        duct.wave_energy_j().unwrap() + heat,
        supplied_work,
        supplied_work,
    );
    close(
        duct.cumulative_mutual_friction_heat_j(),
        heat,
        supplied_work,
    );
    let middle = duct.face(20).unwrap();
    assert!(middle.normal_velocity_m_s > 0.0);
    assert!(middle.superfluid_velocity_m_s < 0.0);
    close(momentum(&duct), 0.0, 1e-8);
}

fn moving_channel(n: usize) -> TwoFluidChannel {
    let mut duct = channel(n, 0.23);
    let pressure: Vec<_> = (0..n)
        .map(|i| 0.01 * (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
        .collect();
    let temperature: Vec<_> = (0..n)
        .map(|i| 0.0003 * (std::f64::consts::TAU * (i as f64 + 0.5) / n as f64).sin())
        .collect();
    let mut normal = vec![0.0; n + 1];
    let mut superfluid = vec![0.0; n + 1];
    for i in 1..n {
        let shape = (std::f64::consts::PI * i as f64 / n as f64).sin();
        normal[i] = 0.02 * shape;
        superfluid[i] = -0.01 * shape;
    }
    duct.initialize(&pressure, &temperature, &normal, &superfluid, || false)
        .unwrap();
    duct
}

#[test]
fn nonuniform_forcing_has_independent_impulse_and_negative_work_balances() {
    let n = 13;
    let mut duct = moving_channel(n);
    let before_momentum = momentum(&duct);
    let before_energy = duct.wave_energy_j().unwrap();
    let old: Vec<_> = (0..=n).map(|i| duct.face(i).unwrap()).collect();
    let old_left = duct.cell(0).unwrap().pressure_perturbation_pa;
    let old_right = duct.cell(n - 1).unwrap().pressure_perturbation_pa;
    // Prescribed opposing forces, fixed during this interval. Their work must
    // be negative; it is not heat or a forcibly nonnegative source budget.
    let mut forces: Vec<_> = old
        .iter()
        .map(|f| {
            force(
                -0.2 * duct.material().normal_density_kg_m3() * f.normal_velocity_m_s,
                -0.2 * duct.material().superfluid_density_kg_m3() * f.superfluid_velocity_m_s,
            )
        })
        .collect();
    forces[0] = force(0.004, 0.005);
    forces[n] = force(-0.006, 0.002);
    let dt = 0.7 * duct.maximum_step_s();
    let volume = duct.area_m2() * duct.length_m() / n as f64;
    let report = duct
        .step_forced_closed(dt, BodyForce::Faces(&forces), MutualFriction::None, || {
            false
        })
        .unwrap();
    let expected_body: f64 = forces
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let weight = if i == 0 || i == n { 0.5 } else { 1.0 };
            weight * volume * dt * (f.normal_n_m3 + f.superfluid_n_m3)
        })
        .sum();
    let expected_work: f64 = (1..n)
        .map(|i| {
            let new = duct.face(i).unwrap();
            dt * volume
                * (forces[i].normal_n_m3
                    * 0.5
                    * (old[i].normal_velocity_m_s + new.normal_velocity_m_s)
                    + forces[i].superfluid_n_m3
                        * 0.5
                        * (old[i].superfluid_velocity_m_s + new.superfluid_velocity_m_s))
        })
        .sum();
    let pressure_impulses = [
        duct.area_m2() * dt * 0.5 * (old_left + duct.cell(0).unwrap().pressure_perturbation_pa),
        -duct.area_m2()
            * dt
            * 0.5
            * (old_right + duct.cell(n - 1).unwrap().pressure_perturbation_pa),
    ];
    let expected_walls = [
        pressure_impulses[0]
            - 0.5 * volume * dt * (forces[0].normal_n_m3 + forces[0].superfluid_n_m3),
        pressure_impulses[1]
            - 0.5 * volume * dt * (forces[n].normal_n_m3 + forces[n].superfluid_n_m3),
    ];
    let impulse_scale = expected_body.abs() + expected_walls.iter().map(|j| j.abs()).sum::<f64>();
    close(
        report.momentum_before_kg_m_s,
        before_momentum,
        before_momentum.abs(),
    );
    close(
        report.momentum_after_kg_m_s,
        momentum(&duct),
        before_momentum.abs(),
    );
    close(report.body_impulse_n_s, expected_body, impulse_scale);
    for end in 0..2 {
        close(
            report.channel.boundary_pressure_impulse_n_s[end],
            pressure_impulses[end],
            impulse_scale,
        );
        close(
            report.wall_impulse_n_s[end],
            expected_walls[end],
            impulse_scale,
        );
    }
    close(
        momentum(&duct) - before_momentum,
        expected_body + expected_walls.iter().sum::<f64>(),
        impulse_scale,
    );
    close(report.momentum_defect_n_s, 0.0, impulse_scale);
    assert!(expected_work < 0.0);
    close(
        report.channel.body_force_work_j,
        expected_work,
        expected_work.abs(),
    );
    close(
        duct.wave_energy_j().unwrap() - before_energy,
        expected_work,
        before_energy,
    );
    assert!(duct.wave_energy_j().unwrap() < before_energy);
}

fn snapshot(duct: &TwoFluidChannel) -> Vec<f64> {
    let mut values = vec![
        duct.time_s(),
        duct.cumulative_mutual_friction_heat_j(),
        duct.wave_energy_j().unwrap(),
    ];
    values.extend(duct.last_mutual_friction_heat_j());
    for i in 0..duct.cells() {
        let c = duct.cell(i).unwrap();
        values.extend([
            c.pressure_perturbation_pa,
            c.temperature_perturbation_k,
            c.density_perturbation_kg_m3,
            c.specific_entropy_perturbation_j_kg_k,
        ]);
    }
    for i in 0..=duct.cells() {
        let f = duct.face(i).unwrap();
        values.extend([
            f.normal_velocity_m_s,
            f.superfluid_velocity_m_s,
            f.mass_flux_kg_m2_s,
            f.entropy_flux_w_m2_k,
            f.counterflow_heat_flux_w_m2,
        ]);
    }
    values
}

#[test]
fn forced_refusal_and_heat_transaction_never_publish_partial_impulse_or_work() {
    let n = 32;
    let mut duct = moving_channel(n);
    let before = snapshot(&duct);
    let dt = 0.5 * duct.maximum_step_s();
    let friction = MutualFriction::Linear { gamma_kg_m3_s: 3.0 };
    let short: Vec<_> = (0..n).map(|_| force(0.0, 0.0)).collect();
    let mut bad_endpoint: Vec<_> = (0..=n).map(|_| force(0.0, 0.0)).collect();
    bad_endpoint[n].normal_n_m3 = f64::NAN;
    for invalid in [
        BodyForce::Faces(&short),
        BodyForce::Faces(&bad_endpoint),
        BodyForce::Acceleration(f64::NAN),
        BodyForce::Uniform(force(f64::MAX, f64::MAX)),
    ] {
        assert!(
            duct.step_forced_closed(dt, invalid, friction, || false)
                .is_err()
        );
        assert_eq!(snapshot(&duct), before);
    }
    assert_eq!(
        duct.step_forced_closed(dt, BodyForce::Acceleration(0.01), friction, || true),
        Err(ChannelError::Cancelled)
    );
    assert_eq!(snapshot(&duct), before);
    let refused = duct.step_forced_closed_with_heat(
        dt,
        BodyForce::Acceleration(0.01),
        friction,
        || false,
        |report, heat| {
            assert!(report.body_impulse_n_s > 0.0);
            assert_eq!(heat.len(), n);
            Err::<(), _>("thermal refusal")
        },
    );
    assert_eq!(refused, Err(HeatStepError::Heat("thermal refusal")));
    assert_eq!(snapshot(&duct), before);
    let prepared = Cell::new(false);
    let cancelled = duct.step_forced_closed_with_heat(
        dt,
        BodyForce::Acceleration(0.01),
        friction,
        || prepared.get(),
        |_, heat| {
            prepared.set(true);
            Ok::<_, &'static str>(heat.iter().sum::<f64>())
        },
    );
    assert!(prepared.get());
    assert_eq!(
        cancelled,
        Err(HeatStepError::Channel(ChannelError::Cancelled))
    );
    assert_eq!(snapshot(&duct), before);

    let mut baseline = moving_channel(n);
    let expected = baseline
        .step_forced_closed(dt, BodyForce::Acceleration(0.01), friction, || false)
        .unwrap();
    let (actual, heat) = duct
        .step_forced_closed_with_heat(
            dt,
            BodyForce::Acceleration(0.01),
            friction,
            || false,
            |_, heat| Ok::<_, &'static str>(heat.iter().sum::<f64>()),
        )
        .unwrap();
    assert_eq!(actual, expected);
    close(
        heat,
        actual.channel.mutual_friction_heat_j,
        actual.channel.wave_energy_before_j,
    );
    assert_eq!(snapshot(&duct), snapshot(&baseline));
}

#[test]
fn forcing_cannot_silently_close_an_already_moving_endpoint() {
    let n = 8;
    let mut duct = channel(n, 0.3);
    let mut normal = vec![0.0; n + 1];
    normal[0] = 0.001;
    duct.initialize(
        &vec![0.0; n],
        &vec![0.0; n],
        &normal,
        &vec![0.0; n + 1],
        || false,
    )
    .unwrap();
    let before = snapshot(&duct);
    assert!(
        duct.step_forced_closed(
            0.5 * duct.maximum_step_s(),
            BodyForce::None,
            MutualFriction::None,
            || false
        )
        .is_err()
    );
    assert_eq!(snapshot(&duct), before);
}
