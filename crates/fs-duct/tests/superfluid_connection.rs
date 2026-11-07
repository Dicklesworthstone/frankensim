//! G0/G1/G3/G4: native split-channel equivalence, superleaks and atomic coupling.
#![cfg(feature = "superfluid")]

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ChannelError, ComponentForceDensity, ConnectionDrive,
    ConnectionKind, EndVelocities, HeatStepError, LinearizationLimits, MutualFriction,
    TwoFluidChannel, TwoFluidConnection,
};
use fs_material::superfluid::{NormalViscousSpectrum, TwoFluidProperties, TwoFluidState};
use std::cell::Cell;

fn material() -> TwoFluidState {
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: -0.2,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap()
}

fn limits() -> LinearizationLimits {
    LinearizationLimits {
        max_pressure_perturbation_pa: 1.0,
        max_temperature_perturbation_k: 0.02,
        max_component_speed_m_s: 0.1,
    }
}

fn channel(length: f64, area: f64, cells: usize) -> TwoFluidChannel {
    TwoFluidChannel::try_new(material(), length, area, cells, limits()).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:e}"
    );
}

fn snapshot(channel: &TwoFluidChannel) -> Vec<f64> {
    let mut result = vec![
        channel.time_s(),
        channel.represented_energy_j().unwrap(),
        channel.cumulative_mutual_friction_heat_j(),
        channel.cumulative_normal_viscous_heat_j(),
    ];
    for i in 0..channel.cells() {
        let c = channel.cell(i).unwrap();
        result.extend([c.pressure_perturbation_pa, c.temperature_perturbation_k]);
    }
    for i in 0..=channel.cells() {
        let f = channel.face(i).unwrap();
        result.extend([f.normal_velocity_m_s, f.superfluid_velocity_m_s]);
        if let Some(spectrum) = channel.normal_viscous_spectrum() {
            for mode in 0..spectrum.modes().len() {
                if let Some(value) = channel.normal_profile_mode_velocity_m_s(i, mode) {
                    result.push(value);
                }
            }
        }
    }
    result.extend(channel.last_dissipated_heat_j());
    result
}

#[test]
fn lossless_split_matches_original_unsplit_two_fluid_solver() {
    let n = 8;
    let mut whole = channel(2.0, 0.01, 2 * n);
    let mut first = channel(1.0, 0.01, n);
    let mut second = channel(1.0, 0.01, n);
    let pressure: Vec<_> = (0..2 * n)
        .map(|i| 0.01 * (std::f64::consts::PI * (i as f64 + 0.5) / (2 * n) as f64).cos())
        .collect();
    let temperature: Vec<_> = (0..2 * n)
        .map(|i| 0.0003 * (std::f64::consts::TAU * (i as f64 + 0.5) / (2 * n) as f64).cos())
        .collect();
    let normal: Vec<_> = (0..=2 * n)
        .map(|i| {
            if i == 0 || i == 2 * n {
                0.0
            } else {
                0.0004 * (std::f64::consts::PI * i as f64 / (2 * n) as f64).sin()
            }
        })
        .collect();
    let superfluid: Vec<_> = normal.iter().map(|v| -0.7 * v).collect();
    whole
        .initialize(&pressure, &temperature, &normal, &superfluid, || false)
        .unwrap();
    first
        .initialize(
            &pressure[..n],
            &temperature[..n],
            &normal[..=n],
            &superfluid[..=n],
            || false,
        )
        .unwrap();
    second
        .initialize(
            &pressure[n..],
            &temperature[n..],
            &normal[n..],
            &superfluid[n..],
            || false,
        )
        .unwrap();
    let mut connection =
        TwoFluidConnection::try_new(&first, &second, ConnectionKind::Open).unwrap();
    let energy = whole.wave_energy_j().unwrap();
    near(
        first.wave_energy_j().unwrap()
            + second.wave_energy_j().unwrap()
            + connection.kinetic_energy_j().unwrap(),
        energy,
        1e-20,
    );
    let dt = 0.7 * whole.maximum_step_s();
    for _ in 0..60 {
        whole
            .step(dt, BoundaryVelocities::CLOSED, MutualFriction::None, || {
                false
            })
            .unwrap();
        let report = connection
            .step(
                &mut first,
                &mut second,
                dt,
                ConnectionDrive::default(),
                || false,
            )
            .unwrap();
        for i in 0..2 * n {
            let actual = if i < n {
                first.cell(i)
            } else {
                second.cell(i - n)
            }
            .unwrap();
            let expected = whole.cell(i).unwrap();
            near(
                actual.pressure_perturbation_pa,
                expected.pressure_perturbation_pa,
                2e-12,
            );
            near(
                actual.temperature_perturbation_k,
                expected.temperature_perturbation_k,
                2e-14,
            );
        }
        for i in 1..2 * n {
            let expected = whole.face(i).unwrap();
            let (vn, vs) = if i == n {
                let [u, w] = connection.flow_m3_s().map(|q| q / 0.01);
                (u + 0.75 * w, u - 0.25 * w)
            } else {
                let actual = if i < n {
                    first.face(i)
                } else {
                    second.face(i - n)
                }
                .unwrap();
                (actual.normal_velocity_m_s, actual.superfluid_velocity_m_s)
            };
            near(vn, expected.normal_velocity_m_s, 2e-13);
            near(vs, expected.superfluid_velocity_m_s, 2e-13);
        }
        near(report.combined_energy_defect_j, 0.0, 1e-18);
        near(report.mass_balance_defect_kg, 0.0, 1e-17);
        near(report.entropy_balance_defect_j_k, 0.0, 1e-17);
    }
    near(
        first.wave_energy_j().unwrap()
            + second.wave_energy_j().unwrap()
            + connection.kinetic_energy_j().unwrap(),
        energy,
        5e-12 * energy,
    );
}

#[test]
fn unequal_areas_conserve_each_component_current_and_joint_energy() {
    let mut first = channel(0.7, 0.01, 7);
    let mut second = channel(1.1, 0.004, 5);
    first
        .initialize(&[0.01; 7], &[0.0004; 7], &[0.0; 8], &[0.0; 8], || false)
        .unwrap();
    let mut connection =
        TwoFluidConnection::try_new(&first, &second, ConnectionKind::Open).unwrap();
    let energy = first.wave_energy_j().unwrap();
    let dt = 0.6 * first.maximum_step_s().min(second.maximum_step_s());
    for _ in 0..40 {
        let report = connection
            .step(
                &mut first,
                &mut second,
                dt,
                ConnectionDrive::default(),
                || false,
            )
            .unwrap();
        let a = first.face(first.cells()).unwrap();
        let b = second.face(0).unwrap();
        near(
            first.area_m2() * a.normal_velocity_m_s,
            second.area_m2() * b.normal_velocity_m_s,
            1e-18,
        );
        near(
            first.area_m2() * a.superfluid_velocity_m_s,
            second.area_m2() * b.superfluid_velocity_m_s,
            1e-18,
        );
        near(
            report.channels[0].boundary_mass_kg,
            -report.channels[1].boundary_mass_kg,
            1e-18,
        );
        near(
            report.channels[0].boundary_entropy_j_k,
            -report.channels[1].boundary_entropy_j_k,
            1e-18,
        );
        near(report.combined_energy_defect_j, 0.0, 1e-18);
    }
    near(
        first.wave_energy_j().unwrap()
            + second.wave_energy_j().unwrap()
            + connection.kinetic_energy_j().unwrap(),
        energy,
        1e-11 * energy,
    );
}

#[test]
fn superleak_preserves_fountain_equilibrium_and_drives_toward_hot_side() {
    for pressure_factor in [1.0, 0.0, 2.0] {
        let mut first = channel(1.0, 0.01, 6);
        let mut second = channel(0.8, 0.007, 5);
        let delta_t = 0.0004;
        // Independent thermodynamic prediction: delta p = rho*s*delta T.
        let fountain_pressure = 4.0 * 4.0 * delta_t;
        second
            .initialize(
                &[pressure_factor * fountain_pressure; 5],
                &[delta_t; 5],
                &[0.0; 6],
                &[0.0; 6],
                || false,
            )
            .unwrap();
        let mut connection =
            TwoFluidConnection::try_new(&first, &second, ConnectionKind::Superleak).unwrap();
        let dt = 0.4 * first.maximum_step_s().min(second.maximum_step_s());
        let report = connection
            .step(
                &mut first,
                &mut second,
                dt,
                ConnectionDrive::default(),
                || false,
            )
            .unwrap();
        let a = first.face(first.cells()).unwrap();
        let b = second.face(0).unwrap();
        near(a.normal_velocity_m_s, 0.0, 1e-17);
        near(b.normal_velocity_m_s, 0.0, 1e-17);
        near(report.channels[0].boundary_entropy_j_k, 0.0, 1e-18);
        near(report.channels[1].boundary_entropy_j_k, 0.0, 1e-18);
        let superfluid_flow = -connection.flow_m3_s()[1];
        if pressure_factor == 1.0 {
            near(superfluid_flow, 0.0, 1e-17);
            near(first.cell(0).unwrap().pressure_perturbation_pa, 0.0, 1e-14);
            near(
                second.cell(4).unwrap().pressure_perturbation_pa,
                fountain_pressure,
                1e-14,
            );
        } else {
            assert!(superfluid_flow * (1.0 - pressure_factor) > 0.0);
        }
        let geometric = 0.5 * (1.0 / 6.0) / 0.01 + 0.5 * (0.8 / 5.0) / 0.007;
        near(
            connection.kinetic_energy_j().unwrap(),
            0.5 * 3.0 * geometric * superfluid_flow * superfluid_flow,
            1e-22,
        );
    }
}

#[test]
fn superleak_refuses_to_remove_existing_normal_current() {
    let mut first = channel(1.0, 0.01, 4);
    let mut second = channel(1.0, 0.01, 4);
    for c in [&mut first, &mut second] {
        c.initialize(&[0.0; 4], &[0.0; 4], &[0.001; 5], &[0.0; 5], || false)
            .unwrap();
    }
    assert!(TwoFluidConnection::try_new(&first, &second, ConnectionKind::Superleak).is_err());
    assert!(TwoFluidConnection::try_new(&first, &second, ConnectionKind::Open).is_ok());
}

#[test]
fn body_force_includes_junction_half_volumes_and_preserves_uniform_acceleration() {
    let mut first = channel(1.0, 0.01, 8);
    let mut second = channel(0.5, 0.01, 4);
    let initial = 0.001;
    for c in [&mut first, &mut second] {
        let n = c.cells();
        c.initialize(
            &vec![0.0; n],
            &vec![0.0; n],
            &vec![initial; n + 1],
            &vec![initial; n + 1],
            || false,
        )
        .unwrap();
    }
    let mut connection =
        TwoFluidConnection::try_new(&first, &second, ConnectionKind::Open).unwrap();
    let dt = 0.5 * first.maximum_step_s();
    let acceleration = 0.02;
    let mean = initial + 0.5 * dt * acceleration;
    let end = EndVelocities {
        normal_m_s: mean,
        superfluid_m_s: mean,
    };
    let drive = ConnectionDrive {
        outer: BoundaryVelocities {
            left: end,
            right: end,
        },
        forces: [BodyForce::Acceleration(acceleration); 2],
        ..ConnectionDrive::default()
    };
    let report = connection
        .step(&mut first, &mut second, dt, drive, || false)
        .unwrap();
    near(
        connection.flow_m3_s()[0],
        0.01 * (initial + dt * acceleration),
        1e-18,
    );
    near(connection.flow_m3_s()[1], 0.0, 1e-18);
    near(
        report.junction_body_work_j,
        4.0 * 0.01 * 0.125 * acceleration * mean * dt,
        1e-20,
    );
    for c in [&first, &second] {
        for i in 0..c.cells() {
            near(c.cell(i).unwrap().pressure_perturbation_pa, 0.0, 1e-13);
            near(c.cell(i).unwrap().temperature_perturbation_k, 0.0, 1e-15);
        }
    }
}

fn profile_pair() -> (TwoFluidChannel, TwoFluidChannel, TwoFluidConnection) {
    let make = || {
        let m = material();
        let mut c = TwoFluidChannel::try_new_with_normal_profile(
            m,
            1.0,
            0.01,
            4,
            limits(),
            NormalViscousSpectrum::parallel_plates(m, 0.1, 4).unwrap(),
        )
        .unwrap();
        c.initialize(
            &[0.001, 0.002, -0.002, -0.001],
            &[0.0001; 4],
            &[0.0, 0.002, -0.001, 0.001, 0.0],
            &[0.0; 5],
            || false,
        )
        .unwrap();
        c.step(
            0.002,
            BoundaryVelocities::CLOSED,
            MutualFriction::None,
            || false,
        )
        .unwrap();
        c
    };
    let a = make();
    let b = make();
    let connection = TwoFluidConnection::try_new(&a, &b, ConnectionKind::Open).unwrap();
    (a, b, connection)
}

#[test]
fn rejected_heat_and_final_cancellation_preserve_both_native_histories() {
    let (mut first, mut second, mut connection) = profile_pair();
    let snapshots = [snapshot(&first), snapshot(&second)];
    let before = (
        connection.time_s(),
        connection.flow_m3_s(),
        connection.kinetic_energy_j().unwrap(),
    );
    let drive = ConnectionDrive {
        friction: [MutualFriction::Linear { gamma_kg_m3_s: 3.0 }; 2],
        ..ConnectionDrive::default()
    };
    let refusal = connection.step_with_heat(
        &mut first,
        &mut second,
        0.002,
        drive,
        || false,
        |_, _, _| Err::<(), _>("thermal refusal"),
    );
    assert!(matches!(
        refusal,
        Err(HeatStepError::Heat("thermal refusal"))
    ));
    let cancelled = Cell::new(false);
    let refusal = connection.step_with_heat(
        &mut first,
        &mut second,
        0.002,
        drive,
        || cancelled.get(),
        |_, _, _| {
            cancelled.set(true);
            Ok::<_, ()>(())
        },
    );
    assert!(matches!(
        refusal,
        Err(HeatStepError::Channel(ChannelError::Cancelled))
    ));
    assert_eq!(snapshot(&first), snapshots[0]);
    assert_eq!(snapshot(&second), snapshots[1]);
    assert_eq!(
        (
            connection.time_s(),
            connection.flow_m3_s(),
            connection.kinetic_energy_j().unwrap()
        ),
        before
    );
    let report = connection
        .step(&mut first, &mut second, 0.002, drive, || false)
        .unwrap();
    assert!(report.channels.iter().all(|r| r.dissipated_heat_j() > 0.0));
    near(report.combined_energy_defect_j, 0.0, 1e-18);
}

#[test]
fn second_channel_domain_refusal_does_not_publish_first_or_connection() {
    let mut first = channel(1.0, 0.01, 4);
    let mut tight = limits();
    tight.max_pressure_perturbation_pa = 1e-10;
    tight.max_temperature_perturbation_k = 1e-12;
    let mut second = TwoFluidChannel::try_new(material(), 1.0, 0.01, 4, tight).unwrap();
    first
        .initialize(&[0.01; 4], &[0.001; 4], &[0.0; 5], &[0.0; 5], || false)
        .unwrap();
    let mut connection =
        TwoFluidConnection::try_new(&first, &second, ConnectionKind::Open).unwrap();
    let snapshots = [snapshot(&first), snapshot(&second)];
    assert!(
        connection
            .step(
                &mut first,
                &mut second,
                0.005,
                ConnectionDrive::default(),
                || false
            )
            .is_err()
    );
    assert_eq!(snapshot(&first), snapshots[0]);
    assert_eq!(snapshot(&second), snapshots[1]);
    assert_eq!(connection.time_s(), 0.0);
    assert_eq!(connection.flow_m3_s(), [0.0; 2]);
}

#[test]
fn instantaneous_junction_speed_is_admitted_independently_of_midpoint_speed() {
    let mut tight = limits();
    tight.max_component_speed_m_s = 0.001;
    let mut first = TwoFluidChannel::try_new(material(), 1.0, 0.01, 4, tight).unwrap();
    let mut second = TwoFluidChannel::try_new(material(), 1.0, 0.01, 4, tight).unwrap();
    let mut connection =
        TwoFluidConnection::try_new(&first, &second, ConnectionKind::Open).unwrap();
    let before = [snapshot(&first), snapshot(&second)];
    let zero = ComponentForceDensity {
        normal_n_m3: 0.0,
        superfluid_n_m3: 0.0,
    };
    let drive_force = ComponentForceDensity {
        normal_n_m3: 0.3,
        superfluid_n_m3: 0.9,
    };
    let mut left_force = [zero; 5];
    let mut right_force = [zero; 5];
    left_force[4] = drive_force;
    right_force[0] = drive_force;
    // Only the junction half volumes accelerate. The midpoint speed is below
    // 0.001 m/s, while its instantaneous endpoint would exceed that limit.
    let drive = ConnectionDrive {
        forces: [
            BodyForce::Faces(&left_force),
            BodyForce::Faces(&right_force),
        ],
        ..ConnectionDrive::default()
    };
    let result = connection.step(&mut first, &mut second, 0.005, drive, || false);
    assert!(
        matches!(
            result,
            Err(ChannelError::Domain {
                field: "junction endpoint velocity",
                ..
            })
        ),
        "{result:?}"
    );
    assert_eq!(snapshot(&first), before[0]);
    assert_eq!(snapshot(&second), before[1]);
    assert_eq!(connection.flow_m3_s(), [0.0; 2]);
    assert_eq!(connection.time_s(), 0.0);
}
