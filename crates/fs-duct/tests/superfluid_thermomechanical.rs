#![cfg(feature = "superfluid")]

//! G0/G1/G4: coupled port efforts, physical currents and atomic response probes.

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ChannelError, ComponentForceDensity, LinearizationLimits,
    MutualFriction, ThermomechanicalEndsResponse, TwoFluidChannel,
};
use fs_material::superfluid::{
    NormalViscousSpectrum, NormalViscousWall, TwoFluidProperties, TwoFluidState,
};
use std::cell::Cell;

#[derive(Clone, Copy, Debug)]
enum Wall {
    None,
    Darcy,
    Profile,
}

fn channel(wall: Wall, n: usize, length: f64, area: f64) -> TwoFluidChannel {
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
    let limits = LinearizationLimits {
        max_pressure_perturbation_pa: 1.0,
        max_temperature_perturbation_k: 0.02,
        max_component_speed_m_s: 0.1,
    };
    let radius = (area / core::f64::consts::PI).sqrt();
    let mut duct = match wall {
        Wall::None => TwoFluidChannel::try_new(material, length, area, n, limits),
        Wall::Darcy => TwoFluidChannel::try_new_with_normal_wall(
            material,
            length,
            area,
            n,
            limits,
            NormalViscousWall::circular_capillary(radius).unwrap(),
        ),
        Wall::Profile => TwoFluidChannel::try_new_with_normal_profile(
            material,
            length,
            area,
            n,
            limits,
            NormalViscousSpectrum::circular_capillary(material, radius, 4).unwrap(),
        ),
    }
    .unwrap();
    let pressure: Vec<_> = (0..n)
        .map(|i| 0.005 * (core::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
        .collect();
    let temperature: Vec<_> = (0..n)
        .map(|i| 0.0002 * (core::f64::consts::TAU * (i as f64 + 0.5) / n as f64).sin())
        .collect();
    let normal: Vec<_> = (0..=n)
        .map(|i| 0.001 * (core::f64::consts::PI * i as f64 / n as f64).sin())
        .collect();
    let superfluid: Vec<_> = normal.iter().map(|v| -0.7 * v).collect();
    duct.initialize(&pressure, &temperature, &normal, &superfluid, || false)
        .unwrap();
    if matches!(wall, Wall::Profile) {
        duct.initialize_normal_profile_steady(|| false).unwrap();
    }
    duct
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * scale,
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:e}"
    );
}

fn snapshot(duct: &TwoFluidChannel) -> Vec<f64> {
    let mut values = vec![
        duct.time_s(),
        duct.wave_energy_j().unwrap(),
        duct.normal_profile_excess_energy_j().unwrap(),
        duct.cumulative_mutual_friction_heat_j(),
        duct.cumulative_normal_viscous_heat_j(),
    ];
    for i in 0..duct.cells() {
        let cell = duct.cell(i).unwrap();
        values.extend([
            cell.pressure_perturbation_pa,
            cell.temperature_perturbation_k,
        ]);
    }
    for i in 0..=duct.cells() {
        let face = duct.face(i).unwrap();
        values.extend([face.normal_velocity_m_s, face.superfluid_velocity_m_s]);
        if i > 0 && i < duct.cells() {
            if let Some(spectrum) = duct.normal_viscous_spectrum() {
                for mode in 0..spectrum.modes().len() {
                    values.push(duct.normal_profile_mode_velocity_m_s(i, mode).unwrap());
                }
            }
        }
    }
    values.extend(duct.last_mutual_friction_heat_j());
    values.extend(duct.last_normal_viscous_heat_j());
    values
}

fn check_reciprocal_passive(response: ThermomechanicalEndsResponse) {
    let r = response.impedance_pa_s_m3;
    for i in 0..4 {
        assert!(r[i][i] > 0.0);
        for j in 0..4 {
            assert_eq!(r[i][j], r[j][i]);
        }
    }
    // Finite thermal expansion must not collapse to separate mechanical and
    // thermal ports. Raw temperature is not the second conjugate effort.
    assert!(r[0][1].abs() > 1e-6 * (r[0][0] * r[1][1]).sqrt());
    for direction in [
        [1.0, 1.0, 1.0, 1.0],
        [1.0, -1.0, 1.0, -1.0],
        [1.0, 1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0, 1.0],
        [-0.3, 0.8, 0.5, -0.2],
    ] {
        // Diagonal scaling exercises both coordinates without the much larger
        // pressure impedance obscuring negative thermal or mixed work.
        let q = core::array::from_fn::<_, 4, _>(|i| direction[i] / r[i][i].sqrt());
        let mut work = 0.0;
        let mut scale = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                let term = q[i] * r[i][j] * q[j];
                work += term;
                scale += term.abs();
            }
        }
        assert!(
            work >= -128.0 * f64::EPSILON * scale,
            "negative port work: {work:e}"
        );
    }
}

#[test]
fn joint_response_predicts_every_forced_port_coordinate_with_native_losses() {
    // Different dx and area prevent unit-area or velocity/volume-flow aliases.
    for (wall, n, length, area) in [
        (Wall::None, 2, 1.0, 0.01),
        (Wall::Darcy, 7, 0.8, 0.003),
        (Wall::Profile, 11, 1.3, 0.0008),
    ] {
        for friction in [
            MutualFriction::None,
            MutualFriction::Linear { gamma_kg_m3_s: 3.0 },
            MutualFriction::IsotropicVortexTangle {
                circulation_quantum_m2_s: 1e-7,
                line_density_per_m2: 3e8,
                dissipative_b: 1.2,
            },
            MutualFriction::GorterMellink { a_m_s_kg: 20.0 },
        ] {
            let mut query = channel(wall, n, length, area);
            let forces: Vec<_> = (0..=n)
                .map(|i| ComponentForceDensity {
                    normal_n_m3: 0.001 * (i as f64).cos(),
                    superfluid_n_m3: -0.0004 * (i as f64).sin(),
                })
                .collect();
            let force = BodyForce::Faces(&forces);
            let dt = 0.6 * query.maximum_step_s();
            let reference = [
                [-0.001 * area, 0.002 * area],
                [0.0015 * area, -0.001 * area],
            ];
            let before = snapshot(&query);
            // Exercise duration-dependent factor/cache invalidation first.
            query
                .thermomechanical_ends_response(0.5 * dt, reference, force, friction, || false)
                .unwrap();
            let response = query
                .thermomechanical_ends_response(dt, reference, force, friction, || false)
                .unwrap();
            assert_eq!(snapshot(&query), before);
            assert_eq!(response.reference_flow_outward_m3_s, reference);
            check_reciprocal_passive(response);

            // Independent real steps isolate every matrix column; the fifth
            // perturbation checks simultaneous mixed currents at both ends.
            for direction in [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, -1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, -1.0],
                [0.7, -0.4, -0.6, 0.9],
            ] {
                let delta = direction.map(|value| 0.0005 * area * value);
                let flows = core::array::from_fn(|end| {
                    core::array::from_fn(|component| {
                        reference[end][component] + delta[2 * end + component]
                    })
                });
                let mut actual = channel(wall, n, length, area);
                let old = [actual.cell(0).unwrap(), actual.cell(n - 1).unwrap()];
                let boundary = actual.thermomechanical_boundary_velocities(flows).unwrap();
                let report = actual
                    .step_forced(dt, boundary, force, friction, || false)
                    .unwrap();
                let entropy_coefficient = actual.material().superfluid_density_kg_m3()
                    * actual.material().entropy_j_kg_k();
                let mut work = 0.0;
                for (end, cell) in [0, n - 1].into_iter().enumerate() {
                    let new = actual.cell(cell).unwrap();
                    let observed = [
                        0.5 * (old[end].pressure_perturbation_pa + new.pressure_perturbation_pa),
                        0.5 * entropy_coefficient
                            * (old[end].temperature_perturbation_k
                                + new.temperature_perturbation_k),
                    ];
                    for component in 0..2 {
                        let predicted = response.reference_effort_pa[end][component]
                            - response.impedance_pa_s_m3[2 * end + component]
                                .iter()
                                .zip(delta)
                                .map(|(r, dq)| r * dq)
                                .sum::<f64>();
                        near(observed[component], predicted, 0.1);
                        work -= dt * observed[component] * flows[end][component];
                    }
                }
                near(report.boundary_wave_work_j, work, work.abs().max(1e-12));
            }
        }
    }
}

#[test]
fn joint_flow_mapping_preserves_mass_advected_entropy_and_counterflow_heat() {
    let mut duct = channel(Wall::None, 9, 1.2, 0.003);
    let material = *duct.material();
    let rho = material.density_kg_m3();
    let s = material.entropy_j_kg_k();
    let rs_s = material.superfluid_density_kg_m3() * s;
    let temperature = material.temperature_k();
    let inventory = |duct: &TwoFluidChannel| {
        let volume = duct.area_m2() * duct.length_m() / duct.cells() as f64;
        let mut total = [0.0; 2];
        for i in 0..duct.cells() {
            let cell = duct.cell(i).unwrap();
            total[0] += volume * cell.density_perturbation_kg_m3;
            total[1] += volume
                * (rho * cell.specific_entropy_perturbation_j_kg_k
                    + s * cell.density_perturbation_kg_m3);
        }
        total
    };
    let before = inventory(&duct);
    let flows = [[-5e-6, 12e-6], [8e-6, -3e-6]];
    let boundary = duct.thermomechanical_boundary_velocities(flows).unwrap();
    let dt = 0.6 * duct.maximum_step_s();
    let report = duct
        .step(dt, boundary, MutualFriction::None, || false)
        .unwrap();
    for (end, face, sign) in [(0, 0, -1.0), (1, duct.cells(), 1.0)] {
        let sample = duct.face(face).unwrap();
        near(
            sign * duct.area_m2() * sample.mass_flux_kg_m2_s,
            rho * flows[end][0],
            1e-4,
        );
        near(
            sign * duct.area_m2() * sample.entropy_flux_w_m2_k,
            rho * s * flows[end][0] + rs_s * flows[end][1],
            1e-3,
        );
        near(
            sign * duct.area_m2() * sample.counterflow_heat_flux_w_m2,
            temperature * rs_s * flows[end][1],
            1e-3,
        );
    }
    let mass = -dt * rho * (flows[0][0] + flows[1][0]);
    let heat = -dt * temperature * rs_s * (flows[0][1] + flows[1][1]);
    let entropy = s * mass + heat / temperature;
    near(report.boundary_mass_kg, mass, mass.abs());
    near(report.boundary_counterflow_heat_j, heat, heat.abs());
    near(report.boundary_entropy_j_k, entropy, entropy.abs());
    assert!((entropy - heat / temperature).abs() > 0.1 * entropy.abs());
    let after = inventory(&duct);
    near(after[0] - before[0], mass, mass.abs());
    near(after[1] - before[1], entropy, entropy.abs());
}

#[test]
fn cancelled_joint_queries_preserve_profile_fields_heat_and_deterministic_retry() {
    let friction = MutualFriction::Linear { gamma_kg_m3_s: 3.0 };
    let force = BodyForce::Uniform(ComponentForceDensity {
        normal_n_m3: 0.001,
        superfluid_n_m3: -0.0002,
    });
    let make = || {
        let mut duct = channel(Wall::Profile, 9, 0.8, 0.002);
        duct.step_forced(
            0.4 * duct.maximum_step_s(),
            BoundaryVelocities::CLOSED,
            force,
            friction,
            || false,
        )
        .unwrap();
        assert!(duct.cumulative_mutual_friction_heat_j() > 0.0);
        assert!(duct.cumulative_normal_viscous_heat_j() > 0.0);
        assert!(duct.normal_profile_excess_energy_j().unwrap() > 0.0);
        duct
    };
    let mut query = make();
    let before = snapshot(&query);
    let dt = 0.6 * query.maximum_step_s();
    let flows = [[-1e-6, 2e-6], [1.5e-6, -1e-6]];
    let calls = Cell::new(0usize);
    let expected_response = query
        .thermomechanical_ends_response(dt, flows, force, friction, || {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    assert_eq!(snapshot(&query), before);
    assert!(calls.get() > 2);
    let mut control = make();
    let boundary = control.thermomechanical_boundary_velocities(flows).unwrap();
    let expected_step = control
        .step_forced(dt, boundary, force, friction, || false)
        .unwrap();
    let expected_state = snapshot(&control);
    for stop in [1, calls.get() / 2, calls.get()] {
        let mut duct = make();
        let attempts = Cell::new(0usize);
        assert!(matches!(
            duct.thermomechanical_ends_response(dt, flows, force, friction, || {
                attempts.set(attempts.get() + 1);
                attempts.get() >= stop
            }),
            Err(ChannelError::Cancelled)
        ));
        assert_eq!(snapshot(&duct), before);
        let response = duct
            .thermomechanical_ends_response(dt, flows, force, friction, || false)
            .unwrap();
        assert_eq!(response, expected_response);
        assert_eq!(snapshot(&duct), before);
        let actual = duct
            .step_forced(dt, boundary, force, friction, || false)
            .unwrap();
        assert_eq!(actual, expected_step);
        assert_eq!(snapshot(&duct), expected_state);
    }
}
