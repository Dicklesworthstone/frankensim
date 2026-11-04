//! Unsteady normal Stokes profiles: continuum limits, memory and native ports.

use std::cell::Cell;

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ChannelError, ChannelStep, ComponentForceDensity, EndVelocities,
    HeatStepError, LinearizationLimits, MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::{NormalViscousSpectrum, TwoFluidProperties, TwoFluidState};

const SIZE: f64 = 0.01;

#[derive(Clone, Copy, Debug)]
enum Shape {
    Circle,
    Slit,
}

impl Shape {
    fn spectrum(self, material: TwoFluidState, modes: usize) -> NormalViscousSpectrum {
        match self {
            Self::Circle => NormalViscousSpectrum::circular_capillary(material, SIZE, modes),
            Self::Slit => NormalViscousSpectrum::parallel_plates(material, SIZE, modes),
        }
        .unwrap()
    }

    fn area(self) -> f64 {
        match self {
            Self::Circle => core::f64::consts::PI * SIZE * SIZE,
            Self::Slit => 0.1 * SIZE,
        }
    }

    // nu = eta/rho_n = SIZE^2 in these manufactured material states.
    fn dc_mobility(self) -> f64 {
        match self {
            Self::Circle => 1.0 / 8.0,
            Self::Slit => 1.0 / 12.0,
        }
    }

    fn continuum_energy_factor(self) -> f64 {
        match self {
            Self::Circle => 4.0 / 3.0,
            Self::Slit => 6.0 / 5.0,
        }
    }
}

fn material(alpha: f64) -> TwoFluidState {
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: alpha,
        normal_viscosity_pa_s: SIZE * SIZE,
    })
    .unwrap()
}

fn channel(shape: Shape, modes: usize, cells: usize, alpha: f64) -> TwoFluidChannel {
    let material = material(alpha);
    TwoFluidChannel::try_new_with_normal_profile(
        material,
        1.0,
        shape.area(),
        cells,
        LinearizationLimits {
            max_pressure_perturbation_pa: 1.0,
            max_temperature_perturbation_k: 0.02,
            max_component_speed_m_s: 0.1,
        },
        shape.spectrum(material, modes),
    )
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

fn modes(duct: &TwoFluidChannel) -> Vec<(f64, f64)> {
    duct.normal_viscous_spectrum()
        .unwrap()
        .modes()
        .iter()
        .map(|mode| (mode.weight, mode.decay_rate_s_inv))
        .collect()
}

fn volume(duct: &TwoFluidChannel) -> f64 {
    duct.area_m2() * duct.length_m() * (duct.cells() - 1) as f64 / duct.cells() as f64
}

// Solve the independent scalar Stokes initial-value equations analytically.
// Their interval-mean velocity supplies equal end fluxes, isolating transverse
// dynamics from axial compression. The production prepared factors are unused.
fn uniform_source_step(
    duct: &mut TwoFluidChannel,
    states: &mut [f64],
    acceleration: f64,
    dt: f64,
    superfluid: f64,
) -> (ChannelStep, f64) {
    let mut integrated_normal = 0.0;
    for ((weight, rate), state) in modes(duct).into_iter().zip(states.iter_mut()) {
        let loss = -(-rate * dt).exp_m1();
        let duration = loss / rate;
        integrated_normal += weight * (*state * duration + acceleration / rate * (dt - duration));
        *state = *state * (-rate * dt).exp() + acceleration / rate * loss;
    }
    let report = duct
        .step_forced(
            dt,
            equal_ports(integrated_normal / dt, superfluid),
            BodyForce::Uniform(ComponentForceDensity {
                normal_n_m3: duct.material().normal_density_kg_m3() * acceleration,
                superfluid_n_m3: 0.0,
            }),
            MutualFriction::None,
            || false,
        )
        .unwrap();
    for face in 1..duct.cells() {
        for (mode, &expected) in states.iter().enumerate() {
            near(
                duct.normal_profile_mode_velocity_m_s(face, mode).unwrap(),
                expected,
                4e-12 * expected.abs() + 2e-14,
            );
        }
        near(
            duct.face(face).unwrap().superfluid_velocity_m_s,
            superfluid,
            2e-14,
        );
    }
    (report, integrated_normal)
}

fn snapshot(duct: &TwoFluidChannel) -> Vec<f64> {
    let mut values = vec![
        duct.time_s(),
        duct.wave_energy_j().unwrap(),
        duct.normal_profile_excess_energy_j().unwrap(),
        duct.resolved_momentum_kg_m_s().unwrap(),
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
    for face in 0..=duct.cells() {
        let value = duct.face(face).unwrap();
        values.extend([value.normal_velocity_m_s, value.superfluid_velocity_m_s]);
        if face > 0 && face < duct.cells() {
            for mode in 0..modes(duct).len() {
                values.push(duct.normal_profile_mode_velocity_m_s(face, mode).unwrap());
            }
        }
    }
    values.extend(duct.last_mutual_friction_heat_j());
    values.extend(duct.last_normal_viscous_heat_j());
    values
}

fn developed(shape: Shape, cells: usize) -> TwoFluidChannel {
    let mut duct = channel(shape, 8, cells, 0.1);
    let normal: Vec<_> = (0..=cells)
        .map(|i| 0.003 * (core::f64::consts::PI * i as f64 / cells as f64).sin())
        .collect();
    let superfluid: Vec<_> = normal.iter().map(|u| -u / 3.0).collect();
    duct.initialize(
        &vec![0.01; cells],
        &vec![0.0001; cells],
        &normal,
        &superfluid,
        || false,
    )
    .unwrap();
    let dt = 0.4 * duct.maximum_step_s();
    duct.step(
        dt,
        BoundaryVelocities::CLOSED,
        MutualFriction::Linear { gamma_kg_m3_s: 0.5 },
        || false,
    )
    .unwrap();
    assert!(duct.normal_profile_excess_energy_j().unwrap() > 0.0);
    duct
}

#[test]
fn uniform_forcing_startup_refines_to_circular_and_slit_stokes_flow() {
    // Independent continuum values at nu*t/L^2 = [.01,.1,1], L=R or full gap.
    // Circle: 1/8 - sum 4*exp(-j_0m^2*t)/j_0m^4.
    // Slit: 1/12 - sum 8*exp(-(2m+1)^2*pi^2*t)/((2m+1)^4*pi^4).
    // Constants were evaluated with 120-digit Bessel/Fourier series, without
    // the production zero table, tail closure or prepared-source coefficients.
    for (shape, reference) in [
        (
            Shape::Circle,
            [
                0.008546292106486597,
                0.05771930723672449,
                0.12463176935846083,
            ],
        ),
        (
            Shape::Slit,
            [
                0.008495494443872692,
                0.05272349614541028,
                0.08332908541883954,
            ],
        ),
    ] {
        let acceleration = 0.02;
        let mut previous_error = f64::INFINITY;
        for retained in [1, 4, 32] {
            let mut duct = channel(shape, retained, 2, 0.0);
            uniform(&mut duct, 0.0, 0.0);
            let mut states = vec![0.0; modes(&duct).len()];
            let mut worst: f64 = 0.0;
            let mut external_work = 0.0;
            for step in 1..=100 {
                let (report, integrated) =
                    uniform_source_step(&mut duct, &mut states, acceleration, 0.01, 0.0);
                external_work += volume(&duct) * acceleration * integrated;
                near(
                    report.body_force_work_j,
                    volume(&duct) * acceleration * integrated,
                    2e-20,
                );
                near(
                    duct.wave_energy_j().unwrap() + duct.cumulative_normal_viscous_heat_j(),
                    external_work,
                    2e-18,
                );
                if let Some(index) = [1, 10, 100].iter().position(|value| *value == step) {
                    let observed = duct.face(1).unwrap().normal_velocity_m_s / acceleration;
                    let error = observed - reference[index];
                    // The positive moment-matched tail gives the Jensen upper
                    // approximation for this constant-force startup only.
                    assert!(
                        error >= -2e-12,
                        "{shape:?}, modes={retained}, error={error:e}"
                    );
                    worst = worst.max(error.abs());
                }
            }
            assert!(
                worst <= previous_error + 2e-12,
                "{shape:?}: {worst:e} > {previous_error:e}"
            );
            previous_error = worst;
        }
        assert!(previous_error < 2e-10, "{shape:?}: {previous_error:e}");
    }
}

#[test]
fn steady_profiles_recover_poiseuille_resistance_and_true_kinetic_storage() {
    let normal = 0.003;
    let superfluid = -0.001;
    for shape in [Shape::Circle, Shape::Slit] {
        let mut previous_factor = 1.0;
        for retained in [1, 4, 32] {
            let mut duct = channel(shape, retained, 2, 0.0);
            uniform(&mut duct, normal, superfluid);
            assert_eq!(duct.normal_profile_excess_energy_j().unwrap(), 0.0);
            duct.initialize_normal_profile_steady(|| false).unwrap();
            let storage = duct.wave_energy_j().unwrap();
            let superfluid_storage = 0.5 * 3.0 * volume(&duct) * superfluid * superfluid;
            let factor = 2.0 * (storage - superfluid_storage) / (volume(&duct) * normal * normal);
            assert!(factor >= previous_factor - 2e-13);
            assert!(factor <= shape.continuum_energy_factor() + 2e-13);
            previous_factor = factor;
            for (index, (_, rate)) in modes(&duct).into_iter().enumerate() {
                near(
                    duct.normal_profile_mode_velocity_m_s(1, index).unwrap(),
                    normal / shape.dc_mobility() / rate,
                    3e-15,
                );
            }
            let report = duct
                .step_forced(
                    0.04,
                    equal_ports(normal, superfluid),
                    BodyForce::Uniform(ComponentForceDensity {
                        normal_n_m3: normal / shape.dc_mobility(),
                        superfluid_n_m3: 0.0,
                    }),
                    MutualFriction::None,
                    || false,
                )
                .unwrap();
            let work = volume(&duct) * normal * normal / shape.dc_mobility() * 0.04;
            near(duct.face(1).unwrap().normal_velocity_m_s, normal, 3e-15);
            near(duct.wave_energy_j().unwrap(), storage, 2e-20);
            near(report.normal_viscous_heat_j, work, 2e-12 * work);
            near(report.body_force_work_j, work, 2e-12 * work);
            near(
                report.interior_wall_impulse_n_s,
                -work / normal,
                2e-12 * work / normal,
            );
        }
        near(previous_factor, shape.continuum_energy_factor(), 1e-8);
    }
}

#[test]
fn unforced_profile_energy_and_wall_reaction_close_without_heating_the_superfluid() {
    for shape in [Shape::Circle, Shape::Slit] {
        let mut duct = channel(shape, 16, 4, 0.0);
        let normal = 0.02;
        let superfluid = -0.007;
        uniform(&mut duct, normal, superfluid);
        let initial = duct.wave_energy_j().unwrap();
        let mut states = vec![normal; modes(&duct).len()];
        let (report, _) = uniform_source_step(&mut duct, &mut states, 0.0, 0.02, superfluid);
        let weights = modes(&duct);
        let mean: f64 = weights.iter().zip(&states).map(|((w, _), u)| w * u).sum();
        let kinetic: f64 = weights
            .iter()
            .zip(&states)
            .map(|((w, _), u)| w * u * u)
            .sum();
        let variance: f64 = weights
            .iter()
            .zip(&states)
            .map(|((w, _), u)| w * (u - mean).powi(2))
            .sum();
        let full_energy = 0.5 * volume(&duct) * (kinetic + 3.0 * superfluid * superfluid);
        near(duct.wave_energy_j().unwrap(), full_energy, 2e-19);
        near(
            duct.normal_profile_excess_energy_j().unwrap(),
            0.5 * volume(&duct) * variance,
            2e-20,
        );
        assert!(variance > 0.0);
        near(report.normal_viscous_heat_j, initial - full_energy, 2e-19);
        near(
            report.interior_wall_impulse_n_s,
            volume(&duct) * (mean - normal),
            2e-17,
        );
        near(
            duct.last_normal_viscous_heat_j().iter().sum(),
            report.normal_viscous_heat_j,
            2e-20,
        );
        assert!(duct.last_normal_viscous_heat_j().iter().all(|q| *q >= 0.0));
        assert_eq!(report.mutual_friction_heat_j, 0.0);
    }
}

#[test]
fn a_zero_mean_reversed_flow_retains_profile_memory_and_can_rebound() {
    for shape in [Shape::Circle, Shape::Slit] {
        let mut duct = channel(shape, 16, 2, 0.0);
        uniform(&mut duct, 0.0, 0.0);
        let rates = modes(&duct);
        let dt = 0.04;
        let a = 0.05;
        let response: f64 = rates.iter().map(|(w, r)| w * -(-r * dt).exp_m1() / r).sum();
        let carried: f64 = rates
            .iter()
            .map(|(w, r)| w * -(-r * dt).exp_m1() / r * (-r * dt).exp())
            .sum();
        let reverse = -a * carried / response;
        let mut states = vec![0.0; rates.len()];
        uniform_source_step(&mut duct, &mut states, a, dt, 0.0);
        uniform_source_step(&mut duct, &mut states, reverse, dt, 0.0);
        near(duct.face(1).unwrap().normal_velocity_m_s, 0.0, 2e-14);
        let stored = duct.wave_energy_j().unwrap();
        assert!(stored > 1e-14 * volume(&duct));
        near(
            stored,
            duct.normal_profile_excess_energy_j().unwrap(),
            2e-20,
        );
        let (report, _) = uniform_source_step(&mut duct, &mut states, 0.0, dt, 0.0);
        assert!(duct.face(1).unwrap().normal_velocity_m_s > 1e-6);
        assert!(report.normal_viscous_heat_j > 0.0);
        near(
            duct.wave_energy_j().unwrap() + report.normal_viscous_heat_j,
            stored,
            2e-19,
        );
    }
}

#[test]
fn oscillatory_forcing_recovers_womersley_amplitude_and_phase() {
    // Independent circular Bessel and slit tanh solutions of the transverse
    // Stokes PDE, with exp(+i*omega*t), Omega=omega*L^2/nu=10.
    // Circle H=(1-2*J1(z)/(z*J0(z)))/(i*Omega/8), z^2=-i*Omega.
    // Slit H=(1-tanh(z)/z)/(i*Omega/12), z^2=i*Omega/4.
    for (shape, expected) in [
        (Shape::Circle, [0.2795497501957999, -0.4257706707249414]),
        (Shape::Slit, [0.5006104282665160, -0.4941524395021566]),
    ] {
        let omega = 10.0;
        let samples = 256;
        let dt = 2.0 * core::f64::consts::PI / (omega * samples as f64);
        let acceleration = 0.02;
        let mut duct = channel(shape, 32, 2, 0.0);
        uniform(&mut duct, 0.0, 0.0);
        let mut states = vec![0.0; modes(&duct).len()];
        let mut response = [0.0; 2];
        for step in 0..6 * samples {
            let forcing = acceleration * (omega * (step as f64 + 0.5) * dt).cos();
            uniform_source_step(&mut duct, &mut states, forcing, dt, 0.0);
            if step >= 5 * samples {
                let phase = omega * (step + 1) as f64 * dt;
                let velocity = duct.face(1).unwrap().normal_velocity_m_s;
                response[0] += velocity * phase.cos();
                response[1] -= velocity * phase.sin();
            }
        }
        for component in 0..2 {
            let observed =
                2.0 * response[component] / (samples as f64 * acceleration * shape.dc_mobility());
            // Midpoint-held forcing has O((omega*dt)^2) time error; the
            // 32-mode geometry approximation is much smaller at this Omega.
            near(observed, expected[component], 2e-4);
        }
    }
}

#[test]
fn reciprocal_end_queries_preserve_and_advance_the_same_profile_history() {
    let friction = MutualFriction::Linear { gamma_kg_m3_s: 0.5 };
    let force = BodyForce::Uniform(ComponentForceDensity {
        normal_n_m3: 0.001,
        superfluid_n_m3: -0.0002,
    });
    for shape in [Shape::Circle, Shape::Slit] {
        for n in [2, 17] {
            let mut duct = developed(shape, n);
            let before = snapshot(&duct);
            let old = [duct.cell(0).unwrap(), duct.cell(n - 1).unwrap()];
            let dt = 0.5 * duct.maximum_step_s();
            let reference = [0.001, -0.002];
            let delta = [-0.0004, 0.0003];
            let response = duct
                .impermeable_ends_response(dt, reference, force, friction, || false)
                .unwrap();
            assert_eq!(snapshot(&duct), before);
            near(
                response.pressure_velocity_impedance_pa_s_m[0][1],
                response.pressure_velocity_impedance_pa_s_m[1][0],
                0.0,
            );
            let left = -(reference[0] + delta[0]);
            let right = reference[1] + delta[1];
            duct.step_forced(
                dt,
                BoundaryVelocities {
                    left: EndVelocities {
                        normal_m_s: left,
                        superfluid_m_s: left,
                    },
                    right: EndVelocities {
                        normal_m_s: right,
                        superfluid_m_s: right,
                    },
                },
                force,
                friction,
                || false,
            )
            .unwrap();
            for i in 0..2 {
                let predicted = response.reference_pressure_pa[i]
                    - response.pressure_velocity_impedance_pa_s_m[i][0] * delta[0]
                    - response.pressure_velocity_impedance_pa_s_m[i][1] * delta[1];
                let observed = 0.5
                    * (old[i].pressure_perturbation_pa
                        + duct.cell(i * (n - 1)).unwrap().pressure_perturbation_pa);
                near(observed, predicted, 2e-12);
            }

            let mut duct = developed(shape, n);
            let before = snapshot(&duct);
            let old = [duct.cell(0).unwrap(), duct.cell(n - 1).unwrap()];
            let reference = [-1e-6, 0.7e-6];
            let delta = [0.3e-6, -0.2e-6];
            let response = duct
                .counterflow_ends_response(dt, reference, force, friction, || false)
                .unwrap();
            assert_eq!(snapshot(&duct), before);
            near(
                response.temperature_heat_impedance_k_w[0][1],
                response.temperature_heat_impedance_k_w[1][0],
                0.0,
            );
            let boundary = duct
                .counterflow_boundary_velocities([reference[0] + delta[0], reference[1] + delta[1]])
                .unwrap();
            duct.step_forced(dt, boundary, force, friction, || false)
                .unwrap();
            for i in 0..2 {
                let predicted = response.reference_temperature_perturbation_k[i]
                    - response.temperature_heat_impedance_k_w[i][0] * delta[0]
                    - response.temperature_heat_impedance_k_w[i][1] * delta[1];
                let observed = 0.5
                    * (old[i].temperature_perturbation_k
                        + duct.cell(i * (n - 1)).unwrap().temperature_perturbation_k);
                near(observed, predicted, 2e-14);
            }
        }
    }
}

#[test]
fn profile_refusals_and_late_cancellation_preserve_every_accepted_mode() {
    let make = || developed(Shape::Circle, 8);
    let friction = MutualFriction::Linear { gamma_kg_m3_s: 0.5 };
    let mut baseline = make();
    let dt = 0.5 * baseline.maximum_step_s();
    let calls = Cell::new(0usize);
    let expected = baseline
        .step_with_heat(
            dt,
            BoundaryVelocities::CLOSED,
            friction,
            || {
                calls.set(calls.get() + 1);
                false
            },
            |_, _: &[f64]| Ok::<(), ()>(()),
        )
        .unwrap();
    let accepted = snapshot(&baseline);
    for stop in [1, calls.get() / 2, calls.get()] {
        let mut duct = make();
        let before = snapshot(&duct);
        let count = Cell::new(0usize);
        let result = duct.step_with_heat(
            dt,
            BoundaryVelocities::CLOSED,
            friction,
            || {
                count.set(count.get() + 1);
                count.get() >= stop
            },
            |_, _: &[f64]| Ok::<(), ()>(()),
        );
        assert!(
            matches!(result, Err(HeatStepError::Channel(ChannelError::Cancelled))),
            "stop={stop}: {result:?}"
        );
        assert_eq!(snapshot(&duct), before);
        let retry = duct
            .step_with_heat(
                dt,
                BoundaryVelocities::CLOSED,
                friction,
                || false,
                |_, _: &[f64]| Ok::<(), ()>(()),
            )
            .unwrap();
        assert_eq!(retry, expected);
        assert_eq!(snapshot(&duct), accepted);
    }
    let mut duct = make();
    let before = snapshot(&duct);
    let result = duct.step_with_heat(
        dt,
        BoundaryVelocities::CLOSED,
        friction,
        || false,
        |report, _: &[f64]| {
            assert!(report.normal_viscous_heat_j > 0.0);
            Err::<(), _>("original thermal owner refused")
        },
    );
    assert!(matches!(result, Err(HeatStepError::Heat(_))));
    assert_eq!(snapshot(&duct), before);
    assert!(
        duct.step(dt, equal_ports(1.0, 0.0), friction, || false)
            .is_err()
    );
    assert_eq!(snapshot(&duct), before);
    assert!(duct.initialize_normal_profile_steady(|| false).is_err());
    assert_eq!(snapshot(&duct), before);

    let mut fresh = channel(Shape::Slit, 8, 2, 0.0);
    uniform(&mut fresh, 0.099, 0.0);
    let before = snapshot(&fresh);
    // Mean speed is admitted, but its developed Poiseuille profile exceeds
    // the original normal RMS limit; preparation must remain atomic.
    assert!(fresh.initialize_normal_profile_steady(|| false).is_err());
    assert_eq!(snapshot(&fresh), before);
}
