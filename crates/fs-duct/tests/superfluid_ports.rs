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
