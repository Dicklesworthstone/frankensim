//! G1/G3/G4: two-fluid bulk attenuation, physical heat, availability and rollback.

use std::cell::Cell;

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ChannelEnd, ChannelError, ClosedTransportControls,
    ComponentForceDensity, EndVelocities, HeatStepError, LinearizationLimits, MutualFriction,
    TwoFluidChannel,
};
use fs_material::superfluid::transport::{LongitudinalTransport, SecondViscosities};
use fs_material::superfluid::{NormalViscousWall, TwoFluidProperties, TwoFluidState};

fn material(eta: f64, alpha: f64) -> TwoFluidState {
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

fn channel(state: TwoFluidState, n: usize, wall: bool) -> TwoFluidChannel {
    let limits = LinearizationLimits {
        max_pressure_perturbation_pa: 1.0,
        max_temperature_perturbation_k: 0.02,
        max_component_speed_m_s: 0.1,
    };
    if wall {
        TwoFluidChannel::try_new_with_normal_wall(
            state,
            1.0,
            0.01,
            n,
            limits,
            NormalViscousWall::try_new(0.02).unwrap(),
        )
        .unwrap()
    } else {
        TwoFluidChannel::try_new(state, 1.0, 0.01, n, limits).unwrap()
    }
}

fn initialize(duct: &mut TwoFluidChannel, amplitudes: [f64; 4]) {
    let n = duct.cells();
    let fs = duct.material().superfluid_fraction();
    let fn_ = duct.material().normal_fraction();
    let cells: Vec<_> = (0..n)
        .map(|i| (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
        .collect();
    let faces: Vec<_> = (0..=n)
        .map(|i| {
            if i == 0 || i == n {
                0.0
            } else {
                (std::f64::consts::PI * i as f64 / n as f64).sin()
            }
        })
        .collect();
    duct.initialize(
        &cells.iter().map(|v| v * amplitudes[0]).collect::<Vec<_>>(),
        &cells.iter().map(|v| v * amplitudes[1]).collect::<Vec<_>>(),
        &faces
            .iter()
            .map(|v| v * (amplitudes[2] + fs * amplitudes[3]))
            .collect::<Vec<_>>(),
        &faces
            .iter()
            .map(|v| v * (amplitudes[2] - fn_ * amplitudes[3]))
            .collect::<Vec<_>>(),
        || false,
    )
    .unwrap();
}

fn snapshot(duct: &TwoFluidChannel) -> Vec<u64> {
    let mut values = vec![
        duct.time_s(),
        duct.represented_energy_j().unwrap(),
        duct.cumulative_dissipated_heat_j(),
        duct.cumulative_conductive_availability_loss_j(),
    ];
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
        values.extend([f.normal_velocity_m_s, f.superfluid_velocity_m_s]);
    }
    values.extend(duct.last_dissipated_heat_j());
    values.into_iter().map(f64::to_bits).collect()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 2e-11 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

// Independent semidiscrete Landau--Khalatnikov ODE in physical (p,T,v,w),
// without calling the transport matrix, channel modal transform or solver.
fn derivative(
    y: [f64; 4],
    state: TwoFluidState,
    second: [f64; 3],
    conductivity: f64,
    k: f64,
) -> [f64; 4] {
    let rho = state.density_kg_m3();
    let fs = state.superfluid_fraction();
    let rho_s = state.superfluid_density_kg_m3();
    let rho_n = state.normal_density_kg_m3();
    let s = state.entropy_j_kg_k();
    let a = state.isothermal_compressibility_per_pa();
    let alpha = state.thermal_expansion_per_k();
    let d = rho * state.specific_heat_cp_j_kg_k() / state.temperature_k();
    let det = a * d - alpha * alpha;
    let mass = -k * y[2];
    let entropy = -k * rho_s * s * y[3] - conductivity / state.temperature_k() * k * k * y[1];
    let vn = y[2] + fs * y[3];
    let q0 = -rho_s * y[3];
    let stress = (4.0 * state.normal_viscosity_pa_s() / 3.0 + second[1]) * vn + second[0] * q0;
    let chemical = second[0] * vn + second[2] * q0;
    let vdot = k * y[0] / rho - k * k * stress / rho;
    let wdot = k * rho * s / rho_n * y[1] - k * k * (stress - rho * chemical) / rho_n;
    [
        (d * mass + alpha * entropy) / det,
        (alpha * mass + a * entropy) / det,
        vdot,
        wdot,
    ]
}

fn reference(
    initial: [f64; 4],
    state: TwoFluidState,
    second: [f64; 3],
    conductivity: f64,
    k: f64,
    duration: f64,
) -> [f64; 4] {
    let mut y = initial;
    let h = duration / 4000.0;
    for _ in 0..4000 {
        let a = derivative(y, state, second, conductivity, k);
        let b = derivative(
            std::array::from_fn(|i| y[i] + 0.5 * h * a[i]),
            state,
            second,
            conductivity,
            k,
        );
        let c = derivative(
            std::array::from_fn(|i| y[i] + 0.5 * h * b[i]),
            state,
            second,
            conductivity,
            k,
        );
        let d = derivative(
            std::array::from_fn(|i| y[i] + h * c[i]),
            state,
            second,
            conductivity,
            k,
        );
        y = std::array::from_fn(|i| y[i] + h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]));
    }
    y
}

#[test]
fn g1_both_sound_modes_and_bulk_attenuation_converge_to_independent_four_field_ode() {
    let n = 12;
    let initial = [0.02, 0.001, 0.001, 0.003];
    let duration = 0.02;
    let k = 2.0 * (std::f64::consts::PI / (2.0 * n as f64)).sin() * n as f64;
    for alpha in [-0.03, 0.03] {
        let state = material(0.03, alpha);
        let second = [0.02, 0.2, 0.06];
        let law = LongitudinalTransport::try_new(
            state,
            SecondViscosities::try_new(second[0], second[1], second[2]).unwrap(),
            0.4,
        )
        .unwrap();
        let expected = reference(initial, state, second, 0.4, k, duration);
        let mut errors = [0.0; 2];
        for (trial, steps) in [20, 40].into_iter().enumerate() {
            let mut duct = channel(state, n, false);
            duct.configure_closed_transport(law, ClosedTransportControls::default(), || false)
                .unwrap();
            initialize(&mut duct, initial);
            for _ in 0..steps {
                let report = duct
                    .step(
                        duration / steps as f64,
                        BoundaryVelocities::CLOSED,
                        MutualFriction::None,
                        || false,
                    )
                    .unwrap();
                assert!(report.longitudinal_viscous_heat_j > 0.0);
                assert!(report.conductive_availability_loss_j > 0.0);
                assert_eq!(report.transport_substeps, 2);
            }
            for i in 0..n {
                let shape = (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos();
                let c = duct.cell(i).unwrap();
                errors[trial] += ((c.pressure_perturbation_pa - expected[0] * shape) / 0.02)
                    .powi(2)
                    + ((c.temperature_perturbation_k - expected[1] * shape) / 0.001).powi(2);
            }
            for i in 1..n {
                let shape = (std::f64::consts::PI * i as f64 / n as f64).sin();
                let f = duct.face(i).unwrap();
                errors[trial] += ((f.barycentric_velocity_m_s - expected[2] * shape) / 0.001)
                    .powi(2)
                    + ((f.counterflow_velocity_m_s - expected[3] * shape) / 0.003).powi(2);
            }
            errors[trial] = (errors[trial] / n as f64).sqrt();
        }
        assert!(errors[1] < 1e-4, "errors={errors:?}");
        assert!(
            (3.7..4.3).contains(&(errors[0] / errors[1])),
            "errors={errors:?}"
        );
    }
}

#[test]
fn g1_conduction_cannot_be_deposited_as_fictitious_caloric_heat() {
    let state = material(0.0, 0.03);
    let mut duct = channel(state, 12, false);
    duct.configure_closed_transport(
        LongitudinalTransport::try_new(state, SecondViscosities::ZERO, 10.0).unwrap(),
        ClosedTransportControls::default(),
        || false,
    )
    .unwrap();
    initialize(&mut duct, [0.0, 0.002, 0.0, 0.0]);
    let before = duct.wave_energy_j().unwrap();
    let (report, ()) = duct
        .step_with_heat(
            0.001,
            BoundaryVelocities::CLOSED,
            MutualFriction::None,
            || false,
            |report, heat| {
                assert_eq!(report.dissipated_heat_j(), 0.0);
                assert!(heat.iter().all(|&q| q == 0.0));
                Ok::<(), ()>(())
            },
        )
        .unwrap();
    assert!(report.conductive_availability_loss_j > 0.0);
    near(
        before - duct.wave_energy_j().unwrap(),
        report.conductive_availability_loss_j,
        before,
    );
    assert_eq!(duct.cumulative_dissipated_heat_j(), 0.0);
    assert_eq!(
        duct.cumulative_conductive_availability_loss_j(),
        report.conductive_availability_loss_j
    );
}

#[test]
fn g1_bulk_traction_and_heat_compose_with_normal_wall_friction_and_body_force() {
    let state = material(0.03, -0.02);
    let mut duct = channel(state, 10, true);
    duct.configure_closed_transport(
        LongitudinalTransport::try_new(
            state,
            SecondViscosities::try_new(-0.03, 0.2, 0.1).unwrap(),
            2.0,
        )
        .unwrap(),
        ClosedTransportControls::default(),
        || false,
    )
    .unwrap();
    initialize(&mut duct, [0.02, 0.001, 0.004, 0.002]);
    let old_energy = duct.wave_energy_j().unwrap();
    let (forced, deposited) = duct
        .step_forced_closed_with_heat(
            0.001,
            BodyForce::Uniform(ComponentForceDensity {
                normal_n_m3: 0.003,
                superfluid_n_m3: -0.001,
            }),
            MutualFriction::Linear {
                gamma_kg_m3_s: 0.03,
            },
            || false,
            |_, heat| Ok::<f64, ()>(heat.iter().sum()),
        )
        .unwrap();
    let report = forced.channel;
    assert!(report.longitudinal_viscous_heat_j > 0.0);
    assert!(report.normal_viscous_heat_j > 0.0);
    assert!(report.mutual_friction_heat_j > 0.0);
    assert!(
        report
            .boundary_longitudinal_impulse_n_s
            .iter()
            .any(|&p| p != 0.0)
    );
    near(
        deposited,
        report.dissipated_heat_j(),
        report.dissipated_heat_j(),
    );
    near(duct.last_dissipated_heat_j().sum(), deposited, deposited);
    near(
        old_energy + report.body_force_work_j - duct.wave_energy_j().unwrap(),
        deposited + report.conductive_availability_loss_j,
        old_energy,
    );
    near(
        forced.momentum_after_kg_m_s - forced.momentum_before_kg_m_s,
        forced.body_impulse_n_s
            + forced.wall_impulse_n_s.iter().sum::<f64>()
            + report.interior_wall_impulse_n_s,
        forced.momentum_before_kg_m_s.abs(),
    );
}

#[test]
fn g3_zero_transport_retains_original_wave_arithmetic() {
    let state = material(0.0, 0.02);
    let mut original = channel(state, 8, false);
    let mut augmented = channel(state, 8, false);
    augmented
        .configure_closed_transport(
            LongitudinalTransport::try_new(state, SecondViscosities::ZERO, 0.0).unwrap(),
            ClosedTransportControls::default(),
            || false,
        )
        .unwrap();
    initialize(&mut original, [0.01, 0.001, 0.001, 0.002]);
    initialize(&mut augmented, [0.01, 0.001, 0.001, 0.002]);
    for _ in 0..4 {
        let old = original
            .step(
                0.001,
                BoundaryVelocities::CLOSED,
                MutualFriction::None,
                || false,
            )
            .unwrap();
        let new = augmented
            .step(
                0.001,
                BoundaryVelocities::CLOSED,
                MutualFriction::None,
                || false,
            )
            .unwrap();
        assert_eq!(old, new);
        assert_eq!(snapshot(&original), snapshot(&augmented));
    }
}

#[test]
fn g4_refused_heat_cancellation_and_diffusion_budget_publish_nothing() {
    let state = material(0.03, 0.02);
    let law = LongitudinalTransport::try_new(
        state,
        SecondViscosities::try_new(0.01, 0.2, 0.1).unwrap(),
        10.0,
    )
    .unwrap();
    let mut duct = channel(state, 8, false);
    duct.configure_closed_transport(law, ClosedTransportControls::default(), || false)
        .unwrap();
    initialize(&mut duct, [0.01, 0.001, 0.001, 0.002]);
    let before = snapshot(&duct);
    let refused = duct.step_with_heat(
        0.001,
        BoundaryVelocities::CLOSED,
        MutualFriction::None,
        || false,
        |_, _| Err::<(), _>("caloric refusal"),
    );
    assert!(matches!(
        refused,
        Err(HeatStepError::Heat("caloric refusal"))
    ));
    assert_eq!(before, snapshot(&duct));
    let count = Cell::new(0);
    let cancelled = duct.step(
        0.001,
        BoundaryVelocities::CLOSED,
        MutualFriction::None,
        || {
            count.set(count.get() + 1);
            count.get() >= 15
        },
    );
    assert_eq!(cancelled, Err(ChannelError::Cancelled));
    assert_eq!(before, snapshot(&duct));
    let stiff = LongitudinalTransport::try_new(state, SecondViscosities::ZERO, 1e8).unwrap();
    duct.configure_closed_transport(
        stiff,
        ClosedTransportControls {
            max_diffusive_courant: 0.5,
            max_substeps: 1,
        },
        || false,
    )
    .unwrap();
    assert!(
        duct.step(
            0.001,
            BoundaryVelocities::CLOSED,
            MutualFriction::None,
            || false
        )
        .is_err()
    );
    assert_eq!(before, snapshot(&duct));
    duct.configure_closed_transport(law, ClosedTransportControls::default(), || false)
        .unwrap();
    let accepted = duct
        .step(
            0.001,
            BoundaryVelocities::CLOSED,
            MutualFriction::None,
            || false,
        )
        .unwrap();
    near(
        duct.cumulative_longitudinal_viscous_heat_j(),
        accepted.longitudinal_viscous_heat_j,
        accepted.longitudinal_viscous_heat_j,
    );
}

#[test]
fn g0_open_ports_and_wrong_material_refuse_explicitly() {
    let state = material(0.03, 0.02);
    let mut duct = channel(state, 8, false);
    let wrong =
        LongitudinalTransport::try_new(material(0.04, 0.02), SecondViscosities::ZERO, 1.0).unwrap();
    assert!(
        duct.configure_closed_transport(wrong, ClosedTransportControls::default(), || false)
            .is_err()
    );
    let law = LongitudinalTransport::try_new(state, SecondViscosities::ZERO, 1.0).unwrap();
    duct.configure_closed_transport(law, ClosedTransportControls::default(), || false)
        .unwrap();
    let before = snapshot(&duct);
    assert!(
        duct.step(
            0.001,
            BoundaryVelocities {
                left: EndVelocities {
                    normal_m_s: 0.001,
                    superfluid_m_s: 0.001
                },
                right: EndVelocities::default(),
            },
            MutualFriction::None,
            || false
        )
        .is_err()
    );
    assert!(
        duct.thermomechanical_ends_response(
            0.001,
            [[0.0; 2]; 2],
            BodyForce::None,
            MutualFriction::None,
            || false
        )
        .is_err()
    );
    assert_eq!(before, snapshot(&duct));
}

#[test]
fn g0_single_end_response_refuses_bulk_transport_without_publishing() {
    let state = material(0.03, 0.02);
    let mut duct = channel(state, 8, false);
    initialize(&mut duct, [0.01, 0.001, 0.001, 0.002]);
    for end in [ChannelEnd::Left, ChannelEnd::Right] {
        assert!(
            duct.impermeable_end_response(0.001, end, MutualFriction::None, || false)
                .is_ok()
        );
    }
    duct.configure_closed_transport(
        LongitudinalTransport::try_new(state, SecondViscosities::ZERO, 1.0).unwrap(),
        ClosedTransportControls::default(),
        || false,
    )
    .unwrap();
    // Refusing a response must preserve already accepted heat and fields too.
    let accepted = duct
        .step(
            0.001,
            BoundaryVelocities::CLOSED,
            MutualFriction::None,
            || false,
        )
        .unwrap();
    assert!(accepted.longitudinal_viscous_heat_j > 0.0);
    assert!(accepted.conductive_availability_loss_j > 0.0);
    let before = snapshot(&duct);
    for end in [ChannelEnd::Left, ChannelEnd::Right] {
        assert!(matches!(
            duct.impermeable_end_response(0.001, end, MutualFriction::None, || false),
            Err(ChannelError::Invalid(
                "configured bulk transport excludes endpoint response queries"
            ))
        ));
        assert_eq!(before, snapshot(&duct));
    }
}
