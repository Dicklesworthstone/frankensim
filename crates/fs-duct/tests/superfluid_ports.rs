//! Reciprocal moving ports, authored momentum, and nonpublishing response queries.

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ChannelError, ComponentForceDensity, EndVelocities,
    HeatStepError, LinearizationLimits, MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState};
use std::cell::Cell;

fn channel(n: usize, tight: bool) -> TwoFluidChannel {
    let material = TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: -0.2,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap();
    TwoFluidChannel::try_new(
        material,
        1.0,
        0.01,
        n,
        LinearizationLimits {
            max_pressure_perturbation_pa: if tight { 1e-10 } else { 1.0 },
            max_temperature_perturbation_k: if tight { 1e-12 } else { 0.02 },
            max_component_speed_m_s: 0.1,
        },
    )
    .unwrap()
}

fn ports(outward: [f64; 2]) -> BoundaryVelocities {
    BoundaryVelocities {
        left: EndVelocities {
            normal_m_s: -outward[0],
            superfluid_m_s: -outward[0],
        },
        right: EndVelocities {
            normal_m_s: outward[1],
            superfluid_m_s: outward[1],
        },
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn interior_momentum(channel: &TwoFluidChannel) -> f64 {
    let volume = channel.area_m2() * channel.length_m() / channel.cells() as f64;
    (1..channel.cells())
        .map(|i| {
            let face = channel.face(i).unwrap();
            volume
                * (channel.material().normal_density_kg_m3() * face.normal_velocity_m_s
                    + channel.material().superfluid_density_kg_m3() * face.superfluid_velocity_m_s)
        })
        .sum()
}

#[test]
fn reciprocal_two_end_response_predicts_actual_forced_pressure_and_momentum() {
    for n in [2, 9, 64] {
        for friction in [
            MutualFriction::None,
            MutualFriction::Linear { gamma_kg_m3_s: 3.0 },
            MutualFriction::GorterMellink { a_m_s_kg: 20.0 },
            MutualFriction::IsotropicVortexTangle {
                circulation_quantum_m2_s: 1e-7,
                line_density_per_m2: 3e8,
                dissipative_b: 1.2,
            },
        ] {
            let mut duct = channel(n, false);
            let pressure: Vec<_> = (0..n)
                .map(|i| 0.005 * (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
                .collect();
            let temperature: Vec<_> = (0..n)
                .map(|i| 0.0002 * (std::f64::consts::TAU * (i as f64 + 0.5) / n as f64).sin())
                .collect();
            let normal: Vec<_> = (0..=n)
                .map(|i| 0.001 * (std::f64::consts::PI * i as f64 / n as f64).sin())
                .collect();
            let superfluid: Vec<_> = normal.iter().map(|v| -0.7 * v).collect();
            duct.initialize(&pressure, &temperature, &normal, &superfluid, || false)
                .unwrap();
            let forces: Vec<_> = (0..=n)
                .map(|i| ComponentForceDensity {
                    normal_n_m3: 0.001 * (i as f64).cos(),
                    superfluid_n_m3: -0.0004 * (i as f64).sin(),
                })
                .collect();
            let force = BodyForce::Faces(&forces);
            let dt = 0.6 * duct.maximum_step_s();
            let before_energy = duct.wave_energy_j().unwrap();
            let before_momentum = interior_momentum(&duct);
            // Querying another duration first exercises factor/cache invalidation.
            duct.impermeable_ends_response(0.5 * dt, [0.0; 2], force, friction, || false)
                .unwrap();
            let reference = [-0.002, 0.003];
            let response = duct
                .impermeable_ends_response(dt, reference, force, friction, || false)
                .unwrap();
            let r = response.pressure_velocity_impedance_pa_s_m;
            assert_eq!(r[0][1], r[1][0]);
            assert!(r[0][0] > r[0][1] && r[1][1] > r[1][0]);
            if n == 2 {
                assert!(r[0][1] > 0.0);
            }
            assert_eq!(duct.time_s(), 0.0);
            assert_eq!(duct.wave_energy_j().unwrap(), before_energy);
            assert_eq!(duct.cumulative_mutual_friction_heat_j(), 0.0);
            assert!(duct.last_mutual_friction_heat_j().iter().all(|q| *q == 0.0));
            let actual_u = [-0.0015, 0.002];
            let report = duct
                .step_forced(dt, ports(actual_u), force, friction, || false)
                .unwrap();
            for (end, cell) in [0, n - 1].into_iter().enumerate() {
                let expected = response.reference_pressure_pa[end]
                    - r[end][0] * (actual_u[0] - reference[0])
                    - r[end][1] * (actual_u[1] - reference[1]);
                let actual =
                    0.5 * (pressure[cell] + duct.cell(cell).unwrap().pressure_perturbation_pa);
                close(actual, expected, 0.1);
                let sign = if end == 0 { 1.0 } else { -1.0 };
                close(
                    report.boundary_pressure_impulse_n_s[end],
                    sign * duct.area_m2() * dt * actual,
                    1e-5,
                );
            }
            let volume = duct.area_m2() * duct.length_m() / n as f64;
            let impulse: f64 = forces[1..n]
                .iter()
                .map(|f| dt * volume * (f.normal_n_m3 + f.superfluid_n_m3))
                .sum();
            let after_momentum = interior_momentum(&duct);
            let scale = before_momentum.abs()
                + after_momentum.abs()
                + report
                    .boundary_pressure_impulse_n_s
                    .iter()
                    .map(|j| j.abs())
                    .sum::<f64>();
            close(
                report.resolved_momentum_before_kg_m_s,
                before_momentum,
                scale,
            );
            close(report.resolved_momentum_after_kg_m_s, after_momentum, scale);
            close(report.interior_body_impulse_n_s, impulse, scale);
            close(
                after_momentum - before_momentum,
                impulse + report.boundary_pressure_impulse_n_s.iter().sum::<f64>(),
                scale,
            );
        }
    }
}

#[test]
fn common_acceleration_and_translation_do_not_create_pressure_or_counterflow() {
    let n = 17;
    let mut duct = channel(n, true);
    let initial_u = 0.01;
    duct.initialize(
        &vec![0.0; n],
        &vec![0.0; n],
        &vec![initial_u; n + 1],
        &vec![initial_u; n + 1],
        || false,
    )
    .unwrap();
    let dt = 0.5 * duct.maximum_step_s();
    let acceleration = 0.02;
    let u_mid = initial_u + 0.5 * acceleration * dt;
    let reference = [-u_mid, u_mid];
    let response = duct
        .impermeable_ends_response(
            dt,
            reference,
            BodyForce::Acceleration(acceleration),
            MutualFriction::None,
            || false,
        )
        .unwrap();
    for p in response.reference_pressure_pa {
        close(p, 0.0, 1e-4);
    }
    let report = duct
        .step_forced(
            dt,
            ports(reference),
            BodyForce::Acceleration(acceleration),
            MutualFriction::None,
            || false,
        )
        .unwrap();
    let mass = duct.material().density_kg_m3() * duct.area_m2() * duct.length_m() * (n - 1) as f64
        / n as f64;
    close(
        report.interior_body_impulse_n_s,
        mass * acceleration * dt,
        1e-5,
    );
    close(
        report.body_force_work_j,
        mass * acceleration * dt * u_mid,
        1e-7,
    );
    assert_eq!(report.mutual_friction_heat_j, 0.0);
    for i in 0..n {
        close(duct.cell(i).unwrap().pressure_perturbation_pa, 0.0, 1e-4);
        close(duct.cell(i).unwrap().temperature_perturbation_k, 0.0, 1e-6);
    }
    for i in 1..n {
        let face = duct.face(i).unwrap();
        close(
            face.normal_velocity_m_s,
            initial_u + acceleration * dt,
            initial_u,
        );
        close(
            face.superfluid_velocity_m_s,
            initial_u + acceleration * dt,
            initial_u,
        );
    }
}

#[test]
fn auxiliary_response_never_admits_an_invalid_real_step_or_publishes_on_cancellation() {
    let mut duct = channel(8, true);
    let dt = 0.5 * duct.maximum_step_s();
    let artificial = [0.001, 0.002];
    let response = duct
        .impermeable_ends_response(
            dt,
            artificial,
            BodyForce::None,
            MutualFriction::None,
            || false,
        )
        .unwrap();
    assert!(
        response
            .reference_pressure_pa
            .iter()
            .any(|p| p.abs() > 1e-10)
    );
    assert_eq!(duct.time_s(), 0.0);
    assert_eq!(duct.wave_energy_j().unwrap(), 0.0);
    assert!(matches!(
        duct.step_forced(
            dt,
            ports(artificial),
            BodyForce::None,
            MutualFriction::None,
            || false
        ),
        Err(ChannelError::Domain { .. })
    ));
    assert_eq!(duct.time_s(), 0.0);
    assert_eq!(duct.wave_energy_j().unwrap(), 0.0);

    let mut duct = channel(8, false);
    let cancelled = Cell::new(false);
    let result = duct.step_forced_with_heat(
        dt,
        ports([0.0001, 0.0002]),
        BodyForce::Acceleration(0.01),
        MutualFriction::Linear { gamma_kg_m3_s: 3.0 },
        || cancelled.get(),
        |_, _| {
            cancelled.set(true);
            Ok::<_, ()>(())
        },
    );
    assert!(matches!(
        result,
        Err(HeatStepError::Channel(ChannelError::Cancelled))
    ));
    assert_eq!(duct.time_s(), 0.0);
    assert_eq!(duct.wave_energy_j().unwrap(), 0.0);
    assert_eq!(duct.cumulative_mutual_friction_heat_j(), 0.0);
    assert!(duct.last_mutual_friction_heat_j().iter().all(|q| *q == 0.0));
}

#[test]
fn counterflow_ports_transport_physical_heat_and_entropy_without_net_mass() {
    let mut duct = channel(32, false);
    let q_out = [0.001, -0.002];
    let boundary = duct.counterflow_boundary_velocities(q_out).unwrap();
    let dt = 0.6 * duct.maximum_step_s();
    let report = duct
        .step(dt, boundary, MutualFriction::None, || false)
        .unwrap();
    let t0 = duct.material().temperature_k();
    let physical_heat = -dt * (q_out[0] + q_out[1]);
    close(
        report.boundary_counterflow_heat_j,
        physical_heat,
        physical_heat,
    );
    close(
        report.boundary_entropy_j_k,
        physical_heat / t0,
        physical_heat / t0,
    );
    assert!(report.boundary_mass_kg.abs() < 1e-20);
    let mut wave_work = 0.0;
    for (end, cell, face, sign) in [(0, 0, 0, -1.0), (1, 31, 32, 1.0)] {
        let velocity = duct.face(face).unwrap();
        let theta_mid = 0.5 * duct.cell(cell).unwrap().temperature_perturbation_k;
        close(
            sign * duct.area_m2() * velocity.counterflow_heat_flux_w_m2,
            q_out[end],
            q_out[end].abs(),
        );
        close(
            sign * duct.area_m2() * velocity.entropy_flux_w_m2_k,
            q_out[end] / t0,
            q_out[end].abs() / t0,
        );
        assert!(velocity.mass_flux_kg_m2_s.abs() < 1e-17);
        wave_work -= dt * theta_mid * q_out[end] / t0;
    }
    close(report.boundary_wave_work_j, wave_work, wave_work.abs());
    assert!(wave_work.abs() < physical_heat.abs() / 100.0);
    let volume = duct.area_m2() * duct.length_m() / duct.cells() as f64;
    let stored_entropy: f64 = (0..duct.cells())
        .map(|i| {
            let sample = duct.cell(i).unwrap();
            volume
                * (duct.material().density_kg_m3() * sample.specific_entropy_perturbation_j_kg_k
                    + duct.material().entropy_j_kg_k() * sample.density_perturbation_kg_m3)
        })
        .sum();
    close(t0 * stored_entropy, physical_heat, physical_heat);
}

#[test]
fn reciprocal_thermal_response_predicts_forced_two_sound_endpoint_temperatures() {
    for n in [2, 9, 64] {
        for friction in [
            MutualFriction::None,
            MutualFriction::Linear { gamma_kg_m3_s: 3.0 },
            MutualFriction::GorterMellink { a_m_s_kg: 20.0 },
            MutualFriction::IsotropicVortexTangle {
                circulation_quantum_m2_s: 1e-7,
                line_density_per_m2: 3e8,
                dissipative_b: 1.2,
            },
        ] {
            let mut duct = channel(n, false);
            let pressure: Vec<_> = (0..n)
                .map(|i| 0.005 * (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
                .collect();
            let temperature: Vec<_> = (0..n)
                .map(|i| 0.0002 * (std::f64::consts::TAU * (i as f64 + 0.5) / n as f64).sin())
                .collect();
            let normal: Vec<_> = (0..=n)
                .map(|i| 0.001 * (std::f64::consts::PI * i as f64 / n as f64).sin())
                .collect();
            let superfluid: Vec<_> = normal.iter().map(|v| -0.7 * v).collect();
            duct.initialize(&pressure, &temperature, &normal, &superfluid, || false)
                .unwrap();
            let force = BodyForce::Uniform(ComponentForceDensity {
                normal_n_m3: 0.002,
                superfluid_n_m3: -0.001,
            });
            let dt = 0.6 * duct.maximum_step_s();
            let reference = [0.0001, -0.0002];
            let q_out = [-0.0003, 0.0004];
            let original = accepted_state(&duct);
            let response = duct
                .counterflow_ends_response(dt, reference, force, friction, || false)
                .unwrap();
            assert_eq!(accepted_state(&duct), original);
            let resistance = response.temperature_heat_impedance_k_w;
            assert_eq!(resistance[0][1], resistance[1][0]);
            assert!(resistance[0][0] > resistance[0][1]);
            assert!(resistance[1][1] > resistance[0][1]);
            let boundary = duct.counterflow_boundary_velocities(q_out).unwrap();
            let report = duct
                .step_forced(dt, boundary, force, friction, || false)
                .unwrap();
            let mut work = 0.0;
            for (end, cell) in [0, n - 1].into_iter().enumerate() {
                assert_eq!(
                    response.reference_temperature_k[end],
                    duct.material().temperature_k()
                        + response.reference_temperature_perturbation_k[end]
                );
                let expected = response.reference_temperature_perturbation_k[end]
                    - resistance[end][0] * (q_out[0] - reference[0])
                    - resistance[end][1] * (q_out[1] - reference[1]);
                let actual =
                    0.5 * (temperature[cell] + duct.cell(cell).unwrap().temperature_perturbation_k);
                close(actual, expected, actual.abs().max(1e-6));
                work -= dt * actual * q_out[end] / duct.material().temperature_k();
            }
            close(report.boundary_wave_work_j, work, work.abs().max(1e-12));
        }
    }
}

#[test]
fn two_cell_thermal_response_matches_independent_caloric_midpoint_equations() {
    let mut properties = channel(2, false).material().properties();
    properties.thermal_expansion_per_k = 0.0;
    let material = TwoFluidState::try_new(properties).unwrap();
    let mut duct = TwoFluidChannel::try_new(
        material,
        1.0,
        0.01,
        2,
        LinearizationLimits {
            max_pressure_perturbation_pa: 1.0,
            max_temperature_perturbation_k: 0.02,
            max_component_speed_m_s: 0.1,
        },
    )
    .unwrap();
    let dt = 0.7 * duct.maximum_step_s();
    let response = duct
        .counterflow_ends_response(dt, [0.0; 2], BodyForce::None, MutualFriction::None, || {
            false
        })
        .unwrap();
    // At alpha=0 the thermal equations are rho*cp*theta_t=-q_x and
    // q_t=-rho*cp*c2^2*theta_x. Eliminating the single interior heat current
    // gives [[1+r,-r],[-r,1+r]] with r=(c2*dt/(2*dx))^2. Boundary Q is in W.
    let dx = 0.5;
    let c2_squared = 3.0 * 2.0 * 16.0 / 16.0;
    let r = c2_squared * (dt / dx).powi(2) / 4.0;
    let heat_capacity =
        material.density_kg_m3() * material.specific_heat_cp_j_kg_k() * duct.area_m2() * dx;
    let factor = dt / (2.0 * heat_capacity * (1.0 + 2.0 * r));
    for i in 0..2 {
        for j in 0..2 {
            let expected = factor * if i == j { 1.0 + r } else { r };
            close(
                response.temperature_heat_impedance_k_w[i][j],
                expected,
                expected,
            );
        }
    }
    assert_eq!(response.reference_temperature_perturbation_k, [0.0; 2]);
    assert_eq!(
        response.reference_temperature_k,
        [material.temperature_k(); 2]
    );
}

#[test]
fn thermal_response_refusals_preserve_state_and_do_not_admit_invalid_real_heat() {
    let mut duct = channel(16, false);
    let original = accepted_state(&duct);
    let dt = 0.5 * duct.maximum_step_s();
    let polls = Cell::new(0);
    assert_eq!(
        duct.counterflow_ends_response(
            dt,
            [0.0001; 2],
            BodyForce::None,
            MutualFriction::None,
            || {
                polls.set(polls.get() + 1);
                polls.get() >= 6
            }
        ),
        Err(ChannelError::Cancelled)
    );
    for q in [[f64::NAN, 0.0], [0.0, f64::INFINITY]] {
        assert!(duct.counterflow_boundary_velocities(q).is_err());
        assert!(
            duct.counterflow_ends_response(dt, q, BodyForce::None, MutualFriction::None, || false)
                .is_err()
        );
    }
    assert_eq!(accepted_state(&duct), original);
    let q = [0.0001, -0.0002];
    let mut fresh = channel(16, false);
    for owner in [&mut duct, &mut fresh] {
        let boundary = owner.counterflow_boundary_velocities(q).unwrap();
        owner
            .step(dt, boundary, MutualFriction::None, || false)
            .unwrap();
    }
    assert_eq!(accepted_state(&duct), accepted_state(&fresh));

    let mut tight = channel(16, true);
    let original = accepted_state(&tight);
    let response = tight
        .counterflow_ends_response(dt, q, BodyForce::None, MutualFriction::None, || false)
        .unwrap();
    assert!(
        response
            .reference_temperature_perturbation_k
            .iter()
            .any(|t| t.abs() > 1e-12)
    );
    assert_eq!(accepted_state(&tight), original);
    let boundary = tight.counterflow_boundary_velocities(q).unwrap();
    assert!(matches!(
        tight.step(dt, boundary, MutualFriction::None, || false),
        Err(ChannelError::Domain { .. })
    ));
    assert_eq!(accepted_state(&tight), original);
}

fn accepted_state(duct: &TwoFluidChannel) -> Vec<f64> {
    let mut values = vec![
        duct.time_s(),
        duct.cumulative_mutual_friction_heat_j(),
        duct.wave_energy_j().unwrap(),
    ];
    values.extend(duct.last_mutual_friction_heat_j());
    for i in 0..duct.cells() {
        let c = duct.cell(i).unwrap();
        values.extend([c.pressure_perturbation_pa, c.temperature_perturbation_k]);
    }
    for i in 0..=duct.cells() {
        let f = duct.face(i).unwrap();
        values.extend([f.normal_velocity_m_s, f.superfluid_velocity_m_s]);
    }
    values
}
