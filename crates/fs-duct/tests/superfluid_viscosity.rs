//! Linear dissipation: exact sources, stiff waves and reciprocal ports.

use std::cell::Cell;

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ComponentForceDensity, EndVelocities, HeatStepError,
    LinearizationLimits, MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::{NormalViscousWall, TwoFluidProperties, TwoFluidState};

fn material(alpha: f64, eta: f64) -> TwoFluidState {
    // Manufactured coefficients make the independent limits below analytic.
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: alpha,
        normal_viscosity_pa_s: eta,
    })
    .unwrap()
}

fn channel(n: usize, material: TwoFluidState, permeability: Option<f64>) -> TwoFluidChannel {
    let limits = LinearizationLimits {
        max_pressure_perturbation_pa: 1.0,
        max_temperature_perturbation_k: 0.02,
        max_component_speed_m_s: 0.1,
    };
    match permeability {
        Some(k) => TwoFluidChannel::try_new_with_normal_wall(
            material,
            1.0,
            0.01,
            n,
            limits,
            NormalViscousWall::try_new(k).unwrap(),
        ),
        None => TwoFluidChannel::try_new(material, 1.0, 0.01, n, limits),
    }
    .unwrap()
}

fn uniform(duct: &mut TwoFluidChannel, normal: f64, superfluid: f64) {
    let n = duct.cells();
    duct.initialize(
        &vec![0.0; n],
        &vec![0.0; n],
        &vec![normal; n + 1],
        &vec![superfluid; n + 1],
        || false,
    )
    .unwrap();
}

fn equal_ports(normal: f64, superfluid: f64) -> BoundaryVelocities {
    let end = EndVelocities {
        normal_m_s: normal,
        superfluid_m_s: superfluid,
    };
    BoundaryVelocities {
        left: end,
        right: end,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:e}"
    );
}

fn snapshot(duct: &TwoFluidChannel) -> Vec<f64> {
    let mut values = vec![
        duct.time_s(),
        duct.wave_energy_j().unwrap(),
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
    }
    values.extend(duct.last_mutual_friction_heat_j());
    values.extend(duct.last_normal_viscous_heat_j());
    values
}

#[test]
fn uniform_normal_relaxation_matches_exact_exponential_including_stiff_limit() {
    let n = 8;
    let dt = 0.005;
    let normal = 0.002;
    let superfluid = -0.0007;
    let resolved_volume = 0.01 * (n - 1) as f64 / n as f64;
    for a in [1e-6_f64, 0.1, 1.0, 1000.0] {
        let mut duct = channel(n, material(0.0, 0.03), Some(0.03 * dt / a));
        uniform(&mut duct, normal, superfluid);
        let mean_normal = normal * (-(-a).exp_m1()) / a;
        let report = duct
            .step(
                dt,
                equal_ports(mean_normal, superfluid),
                MutualFriction::None,
                || false,
            )
            .unwrap();
        for i in 0..n {
            near(duct.cell(i).unwrap().pressure_perturbation_pa, 0.0, 2e-14);
            near(duct.cell(i).unwrap().temperature_perturbation_k, 0.0, 2e-16);
        }
        for i in 1..n {
            let face = duct.face(i).unwrap();
            near(face.normal_velocity_m_s, normal * (-a).exp(), 2e-16);
            near(face.superfluid_velocity_m_s, superfluid, 2e-16);
        }
        let heat = 0.5 * resolved_volume * normal * normal * (-(-2.0 * a).exp_m1());
        let impulse = resolved_volume * normal * (-a).exp_m1();
        near(report.normal_viscous_heat_j, heat, 2e-11 * heat + 1e-25);
        near(
            report.interior_wall_impulse_n_s,
            impulse,
            2e-11 * impulse.abs() + 1e-23,
        );
        assert_eq!(report.mutual_friction_heat_j, 0.0);
        near(
            report.wave_energy_before_j - report.wave_energy_after_j,
            heat,
            2e-20,
        );
        near(duct.last_normal_viscous_heat_j().iter().sum(), heat, 2e-20);
        near(duct.cumulative_dissipated_heat_j(), heat, 2e-20);
        assert!(duct.last_normal_viscous_heat_j().iter().all(|q| *q >= 0.0));
    }
}

#[test]
fn steady_normal_force_balances_wall_drag_without_stiff_kick_error() {
    let n = 8;
    let rate = 1e5;
    let normal = 0.001;
    let dt = 0.005;
    let mut duct = channel(n, material(0.0, 0.03), Some(0.03 / rate));
    uniform(&mut duct, normal, 0.0);
    let force = BodyForce::Uniform(ComponentForceDensity {
        normal_n_m3: rate * normal,
        superfluid_n_m3: 0.0,
    });
    let report = duct
        .step_forced(
            dt,
            equal_ports(normal, 0.0),
            force,
            MutualFriction::None,
            || false,
        )
        .unwrap();
    let work = (n - 1) as f64 / n as f64 * 0.01 * dt * rate * normal * normal;
    for i in 1..n {
        near(duct.face(i).unwrap().normal_velocity_m_s, normal, 3e-15);
        near(duct.face(i).unwrap().superfluid_velocity_m_s, 0.0, 3e-15);
    }
    for i in 0..n {
        near(duct.cell(i).unwrap().pressure_perturbation_pa, 0.0, 3e-13);
        near(duct.cell(i).unwrap().temperature_perturbation_k, 0.0, 3e-15);
    }
    near(report.body_force_work_j, work, 2e-12 * work);
    near(report.normal_viscous_heat_j, work, 2e-12 * work);
    near(
        report.interior_wall_impulse_n_s,
        -work / normal,
        2e-12 * work / normal,
    );
    near(
        report.wave_energy_after_j,
        report.wave_energy_before_j,
        1e-19,
    );
}

#[test]
fn linear_mutual_friction_preserves_driven_equilibrium_at_arbitrary_stiffness() {
    let n = 8;
    let dt = 0.005;
    let bulk = 0.0002;
    let relative = 0.001;
    let resolved_volume = 0.01 * (n - 1) as f64 / n as f64;
    for alpha in [0.0, -0.2] {
        let source = material(alpha, 0.03);
        let normal = bulk + source.superfluid_fraction() * relative;
        let superfluid = bulk - source.normal_fraction() * relative;
        for a in [1e-6_f64, 0.1, 1.0, 1e4] {
            let gamma = source.reduced_density_kg_m3() * a / dt;
            let mut duct = channel(n, source, None);
            uniform(&mut duct, normal, superfluid);
            let report = duct
                .step_forced(
                    dt,
                    equal_ports(normal, superfluid),
                    BodyForce::Uniform(ComponentForceDensity {
                        normal_n_m3: gamma * relative,
                        superfluid_n_m3: -gamma * relative,
                    }),
                    MutualFriction::Linear {
                        gamma_kg_m3_s: gamma,
                    },
                    || false,
                )
                .unwrap();
            for i in 1..n {
                let face = duct.face(i).unwrap();
                near(face.normal_velocity_m_s, normal, 5e-15);
                near(face.superfluid_velocity_m_s, superfluid, 5e-15);
                near(face.barycentric_velocity_m_s, bulk, 5e-15);
            }
            for i in 0..n {
                near(duct.cell(i).unwrap().pressure_perturbation_pa, 0.0, 3e-13);
                near(duct.cell(i).unwrap().temperature_perturbation_k, 0.0, 3e-15);
            }
            let heat = resolved_volume * dt * gamma * relative * relative;
            near(report.mutual_friction_heat_j, heat, 2e-11 * heat + 1e-24);
            near(report.body_force_work_j, heat, 2e-11 * heat + 1e-24);
            near(
                duct.last_mutual_friction_heat_j().iter().sum(),
                heat,
                2e-11 * heat + 1e-24,
            );
            assert_eq!(report.normal_viscous_heat_j, 0.0);
            assert_eq!(report.interior_wall_impulse_n_s, 0.0);
            near(
                report.resolved_momentum_after_kg_m_s,
                report.resolved_momentum_before_kg_m_s,
                1e-18,
            );
            near(
                report.wave_energy_after_j,
                report.wave_energy_before_j,
                1e-18,
            );
        }
    }
}

#[test]
fn forced_linear_counterflow_matches_independent_transient_flux_work_and_heat() {
    let n = 8;
    let dt = 0.005;
    let source = material(-0.2, 0.03);
    let mu = source.reduced_density_kg_m3();
    let rho = source.density_kg_m3();
    let bulk = 0.0001;
    let bulk_acceleration = 0.003;
    let relative = 0.002;
    let equilibrium = -0.0007;
    let difference = relative - equilibrium;
    let resolved_volume = 0.01 * (n - 1) as f64 / n as f64;
    for a in [1e-6_f64, 0.1, 1.0, 1000.0] {
        let rate = a / dt;
        let gamma = mu * rate;
        let relative_acceleration = rate * equilibrium;
        let mean_relative = equilibrium + difference * (-(-a).exp_m1()) / a;
        let final_relative = equilibrium + difference * (-a).exp();
        let mean_bulk = bulk + 0.5 * bulk_acceleration * dt;
        let final_bulk = bulk + bulk_acceleration * dt;
        let mut duct = channel(n, source, None);
        uniform(
            &mut duct,
            bulk + source.superfluid_fraction() * relative,
            bulk - source.normal_fraction() * relative,
        );
        let report = duct
            .step_forced(
                dt,
                equal_ports(
                    mean_bulk + source.superfluid_fraction() * mean_relative,
                    mean_bulk - source.normal_fraction() * mean_relative,
                ),
                BodyForce::Uniform(ComponentForceDensity {
                    normal_n_m3: source.normal_density_kg_m3()
                        * (bulk_acceleration
                            + source.superfluid_fraction() * relative_acceleration),
                    superfluid_n_m3: source.superfluid_density_kg_m3()
                        * (bulk_acceleration - source.normal_fraction() * relative_acceleration),
                }),
                MutualFriction::Linear {
                    gamma_kg_m3_s: gamma,
                },
                || false,
            )
            .unwrap();
        for i in 1..n {
            let face = duct.face(i).unwrap();
            near(face.barycentric_velocity_m_s, final_bulk, 3e-15);
            near(face.counterflow_velocity_m_s, final_relative, 3e-15);
        }
        for i in 0..n {
            near(duct.cell(i).unwrap().pressure_perturbation_pa, 0.0, 3e-13);
            near(duct.cell(i).unwrap().temperature_perturbation_k, 0.0, 3e-15);
        }
        // Independent integral of gamma * (w_eq + delta * exp(-rate*t))^2.
        // exp_m1 retains the weak-drag limit without kinetic-store subtraction.
        let heat = resolved_volume
            * mu
            * (a * equilibrium * equilibrium
                + 2.0 * equilibrium * difference * (-(-a).exp_m1())
                + 0.5 * difference * difference * (-(-2.0 * a).exp_m1()));
        let work = resolved_volume
            * dt
            * (rho * bulk_acceleration * mean_bulk + gamma * equilibrium * mean_relative);
        near(report.mutual_friction_heat_j, heat, 2e-10 * heat + 1e-24);
        near(report.body_force_work_j, work, 2e-10 * work.abs() + 1e-24);
        near(
            report.wave_energy_after_j - report.wave_energy_before_j,
            work - heat,
            1e-18,
        );
        assert_eq!(report.normal_viscous_heat_j, 0.0);
        assert_eq!(report.interior_wall_impulse_n_s, 0.0);
        assert!(
            duct.last_mutual_friction_heat_j()
                .iter()
                .all(|heat| *heat >= 0.0)
        );
    }
}

#[test]
fn fitted_mutual_friction_has_second_order_damped_second_sound_dynamics() {
    let n = 16;
    let rate = 5.0;
    let duration = 0.17;
    let amplitude = 0.003;
    let source = material(0.0, 0.03);
    let shape: Vec<_> = (0..n)
        .map(|i| (core::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
        .collect();
    // Use the exact eigenvalue of the staggered spatial operator, isolating
    // temporal accuracy from spatial dispersion. At alpha=0, c2^2 is 6 m2/s2.
    let wave_number = 2.0 * n as f64 * (core::f64::consts::PI / (2.0 * n as f64)).sin();
    let frequency = (6.0 * wave_number * wave_number - 0.25 * rate * rate).sqrt();
    let exact_amplitude = amplitude
        * (-0.5 * rate * duration).exp()
        * ((frequency * duration).cos() + 0.5 * rate / frequency * (frequency * duration).sin());
    let mut errors = Vec::new();
    for steps in [32, 64, 128] {
        let mut duct = channel(n, source, None);
        duct.initialize(
            &vec![0.0; n],
            &shape.iter().map(|x| amplitude * x).collect::<Vec<_>>(),
            &vec![0.0; n + 1],
            &vec![0.0; n + 1],
            || false,
        )
        .unwrap();
        let initial_energy = duct.wave_energy_j().unwrap();
        for _ in 0..steps {
            let report = duct
                .step(
                    duration / steps as f64,
                    BoundaryVelocities::CLOSED,
                    MutualFriction::Linear {
                        gamma_kg_m3_s: source.reduced_density_kg_m3() * rate,
                    },
                    || false,
                )
                .unwrap();
            assert!(report.mutual_friction_heat_j >= 0.0);
        }
        errors.push(
            (shape
                .iter()
                .enumerate()
                .map(|(i, x)| {
                    (duct.cell(i).unwrap().temperature_perturbation_k - exact_amplitude * x).powi(2)
                })
                .sum::<f64>()
                / n as f64)
                .sqrt(),
        );
        near(
            duct.wave_energy_j().unwrap() + duct.cumulative_mutual_friction_heat_j(),
            initial_energy,
            2e-11 * initial_energy,
        );
    }
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.27 * errors[1], "{errors:?}");
    assert!(errors[2] < 6e-8, "{errors:?}");
}

#[test]
fn changing_friction_at_fixed_dt_refactors_wave_and_endpoint_response() {
    let n = 17;
    let dt = 0.005;
    let source = material(-0.2, 0.03);
    let force = BodyForce::Acceleration(0.001);
    let reference = [0.0003, -0.0001];
    let actual = [-0.0001, 0.0002];
    let mut cached = channel(n, source, None);
    uniform(&mut cached, 0.001, -0.0004);
    let before = snapshot(&cached);
    for friction in [
        MutualFriction::Linear { gamma_kg_m3_s: 1e5 },
        MutualFriction::Linear {
            gamma_kg_m3_s: 0.25,
        },
        MutualFriction::None,
        MutualFriction::IsotropicVortexTangle {
            circulation_quantum_m2_s: 0.03,
            line_density_per_m2: 1e6,
            dissipative_b: 1.0,
        },
        MutualFriction::Linear {
            gamma_kg_m3_s: 1000.0,
        },
    ] {
        let response = cached
            .counterflow_ends_response(dt, reference, force, friction, || false)
            .unwrap();
        let mut fresh = channel(n, source, None);
        uniform(&mut fresh, 0.001, -0.0004);
        assert_eq!(
            response,
            fresh
                .counterflow_ends_response(dt, reference, force, friction, || false)
                .unwrap()
        );
        assert_eq!(snapshot(&cached), before);
        let boundary = fresh.counterflow_boundary_velocities(actual).unwrap();
        fresh
            .step_forced(dt, boundary, force, friction, || false)
            .unwrap();
        let impedance = response.temperature_heat_impedance_k_w;
        assert_eq!(impedance[0][1], impedance[1][0]);
        for (end, cell) in [0, n - 1].into_iter().enumerate() {
            let prediction = response.reference_temperature_perturbation_k[end]
                - impedance[end][0] * (actual[0] - reference[0])
                - impedance[end][1] * (actual[1] - reference[1]);
            near(
                0.5 * fresh.cell(cell).unwrap().temperature_perturbation_k,
                prediction,
                2e-13,
            );
        }
    }
}

#[test]
fn zero_viscosity_retains_original_two_sound_dynamics_exactly() {
    let n = 12;
    let source = material(-0.2, 0.0);
    let mut original = channel(n, source, None);
    let mut configured = channel(n, source, Some(1e-9));
    for duct in [&mut original, &mut configured] {
        uniform(duct, 0.001, -0.002);
    }
    for _ in 0..8 {
        let dt = 0.7 * original.maximum_step_s();
        let run = |duct: &mut TwoFluidChannel| {
            duct.step(
                dt,
                BoundaryVelocities::CLOSED,
                MutualFriction::GorterMellink { a_m_s_kg: 30.0 },
                || false,
            )
            .unwrap()
        };
        assert_eq!(run(&mut original), run(&mut configured));
        assert_eq!(snapshot(&original), snapshot(&configured));
    }
}

#[test]
fn stiff_normal_clamping_propagates_fourth_sound_at_the_physical_speed() {
    // With alpha=0 and vn=0, entropy per volume is frozen and
    // c4²=fs*c1²+fn*c2²=0.75*100+0.25*6=76.5, not c1² or c2².
    let speed = 76.5_f64.sqrt();
    let time = 0.37 / speed;
    let pressure = 0.01;
    let temperature = -0.00125 * pressure;
    let mut errors = Vec::new();
    for n in [12, 24, 48] {
        let mut duct = channel(n, material(0.0, 0.03), Some(3e-10));
        let shape: Vec<_> = (0..n)
            .map(|i| (core::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
            .collect();
        duct.initialize(
            &shape.iter().map(|x| pressure * x).collect::<Vec<_>>(),
            &shape.iter().map(|x| temperature * x).collect::<Vec<_>>(),
            &vec![0.0; n + 1],
            &vec![0.0; n + 1],
            || false,
        )
        .unwrap();
        let initial_energy = duct.wave_energy_j().unwrap();
        let steps = (time / (0.7 * duct.maximum_step_s())).ceil() as usize;
        let dt = time / steps as f64;
        assert!(dt * 1e8 > 1e4, "the test must exercise stiff relaxation");
        for _ in 0..steps {
            let report = duct
                .step(dt, BoundaryVelocities::CLOSED, MutualFriction::None, || {
                    false
                })
                .unwrap();
            assert!(report.normal_viscous_heat_j >= 0.0);
        }
        let phase = (core::f64::consts::PI * speed * time).cos();
        let error = (shape
            .iter()
            .enumerate()
            .map(|(i, x)| {
                (duct.cell(i).unwrap().pressure_perturbation_pa / pressure - x * phase).powi(2)
            })
            .sum::<f64>()
            / n as f64)
            .sqrt();
        errors.push(error);
        near(
            duct.wave_energy_j().unwrap() + duct.cumulative_dissipated_heat_j(),
            initial_energy,
            2e-10 * initial_energy,
        );
        assert!(duct.face(n / 2).unwrap().normal_velocity_m_s.abs() < 2e-9);
    }
    assert!(errors[0] / errors[1] > 3.2, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.2, "{errors:?}");
    assert!(errors[2] < 0.0008, "{errors:?}");
}

#[test]
fn clamped_normal_limit_retains_thermal_expansion_and_entropy_polarization() {
    let n = 64;
    for alpha in [-0.3, 0.3] {
        let source = material(alpha, 0.03);
        // Independently impose delta(rho*s)=0 and eliminate chemical potential
        // in rho'_t=-rho_s*vs_x, vs_t=-(p'/rho-s*T')_x.
        let theta_per_pressure = (alpha
            - source.entropy_j_kg_k()
                * source.density_kg_m3()
                * source.isothermal_compressibility_per_pa())
            / (source.density_kg_m3()
                * (source.specific_heat_cp_j_kg_k() / source.temperature_k()
                    - source.entropy_j_kg_k() * alpha));
        let correction = 1.0 - alpha * 4.0 * 2.0 / 16.0;
        let speed = (0.75 * (100.0 * correction * correction + 2.0)).sqrt();
        let time = 0.37 / speed;
        let mut duct = channel(n, source, Some(3e-10));
        let shape: Vec<_> = (0..n)
            .map(|i| (core::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
            .collect();
        let pressure: Vec<_> = shape.iter().map(|x| 0.01 * x).collect();
        let temperature: Vec<_> = pressure.iter().map(|p| theta_per_pressure * p).collect();
        duct.initialize(
            &pressure,
            &temperature,
            &vec![0.0; n + 1],
            &vec![0.0; n + 1],
            || false,
        )
        .unwrap();
        let steps = (time / (0.7 * duct.maximum_step_s())).ceil() as usize;
        for _ in 0..steps {
            duct.step(
                time / steps as f64,
                BoundaryVelocities::CLOSED,
                MutualFriction::None,
                || false,
            )
            .unwrap();
        }
        let phase = (core::f64::consts::PI * speed * time).cos();
        for i in 0..n {
            let cell = duct.cell(i).unwrap();
            near(cell.pressure_perturbation_pa, pressure[i] * phase, 6e-6);
            near(
                cell.temperature_perturbation_k,
                temperature[i] * phase,
                8e-8,
            );
            let entropy_per_volume = source.density_kg_m3()
                * cell.specific_entropy_perturbation_j_kg_k
                + source.entropy_j_kg_k() * cell.density_perturbation_kg_m3;
            near(entropy_per_volume, 0.0, 2e-7);
        }
    }
}

#[test]
fn viscous_thermal_and_mechanical_responses_predict_the_same_final_operator() {
    for n in [2, 17] {
        for permeability in [0.03, 3e-8] {
            let mut duct = channel(n, material(-0.2, 0.03), Some(permeability));
            uniform(&mut duct, 0.001, -0.0004);
            let friction = MutualFriction::Linear { gamma_kg_m3_s: 2.0 };
            let force = BodyForce::Acceleration(0.001);
            let dt = 0.7 * duct.maximum_step_s();
            let before = snapshot(&duct);
            let reference = [0.0003, -0.0001];
            // Different durations exercise all operator and endpoint caches.
            duct.counterflow_ends_response(0.5 * dt, reference, force, friction, || false)
                .unwrap();
            let response = duct
                .counterflow_ends_response(dt, reference, force, friction, || false)
                .unwrap();
            assert_eq!(before, snapshot(&duct));
            let r = response.temperature_heat_impedance_k_w;
            assert_eq!(r[0][1], r[1][0]);
            assert!(r[0][0] > r[0][1]);
            let q = [-0.0001, 0.0002];
            let boundary = duct.counterflow_boundary_velocities(q).unwrap();
            duct.step_forced(dt, boundary, force, friction, || false)
                .unwrap();
            for (end, cell) in [0, n - 1].into_iter().enumerate() {
                let prediction = response.reference_temperature_perturbation_k[end]
                    - r[end][0] * (q[0] - reference[0])
                    - r[end][1] * (q[1] - reference[1]);
                near(
                    0.5 * duct.cell(cell).unwrap().temperature_perturbation_k,
                    prediction,
                    2e-13,
                );
            }
            let old_pressure = [
                duct.cell(0).unwrap().pressure_perturbation_pa,
                duct.cell(n - 1).unwrap().pressure_perturbation_pa,
            ];
            let uref = [0.001, -0.002];
            let response = duct
                .impermeable_ends_response(0.6 * dt, uref, force, friction, || false)
                .unwrap();
            let u = [0.0007, -0.0013];
            let r = response.pressure_velocity_impedance_pa_s_m;
            assert_eq!(r[0][1], r[1][0]);
            let boundary = BoundaryVelocities {
                left: EndVelocities {
                    normal_m_s: -u[0],
                    superfluid_m_s: -u[0],
                },
                right: EndVelocities {
                    normal_m_s: u[1],
                    superfluid_m_s: u[1],
                },
            };
            duct.step_forced(0.6 * dt, boundary, force, friction, || false)
                .unwrap();
            for (end, cell) in [0, n - 1].into_iter().enumerate() {
                let prediction = response.reference_pressure_pa[end]
                    - r[end][0] * (u[0] - uref[0])
                    - r[end][1] * (u[1] - uref[1]);
                near(
                    0.5 * (old_pressure[end] + duct.cell(cell).unwrap().pressure_perturbation_pa),
                    prediction,
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn viscous_heat_and_momentum_remain_private_on_refusal_and_late_cancellation() {
    let mut duct = channel(8, material(0.1, 0.03), Some(3e-6));
    uniform(&mut duct, 0.001, -0.002);
    let before = snapshot(&duct);
    let dt = 0.6 * duct.maximum_step_s();
    let friction = MutualFriction::Linear { gamma_kg_m3_s: 4.0 };
    let rejected = duct.step_with_heat(
        dt,
        BoundaryVelocities::CLOSED,
        friction,
        || false,
        |report, heat| {
            assert!(report.normal_viscous_heat_j > 0.0);
            assert!(report.mutual_friction_heat_j > 0.0);
            near(heat.iter().sum(), report.dissipated_heat_j(), 1e-20);
            Err::<(), _>("caloric refusal")
        },
    );
    assert!(matches!(
        rejected,
        Err(HeatStepError::Heat("caloric refusal"))
    ));
    assert_eq!(before, snapshot(&duct));
    let cancel = Cell::new(false);
    let rejected = duct.step_with_heat(
        dt,
        BoundaryVelocities::CLOSED,
        friction,
        || cancel.get(),
        |_, _| {
            cancel.set(true);
            Ok::<_, &str>(())
        },
    );
    assert!(matches!(rejected, Err(HeatStepError::Channel(_))));
    assert_eq!(before, snapshot(&duct));
    let report = duct
        .step(dt, BoundaryVelocities::CLOSED, friction, || false)
        .unwrap();
    near(
        duct.cumulative_dissipated_heat_j(),
        report.dissipated_heat_j(),
        1e-20,
    );
    near(
        duct.last_dissipated_heat_j().sum(),
        report.dissipated_heat_j(),
        1e-20,
    );
}
