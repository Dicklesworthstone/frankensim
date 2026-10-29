//! Manufactured two-body/air reference, with actual finite source calorimetry.
//! The ordinary owner follows the complete continuous ODE. Prepared air uses
//! its documented second-order averaged interconnection; refinement is checked
//! against the SAME ODE, not mistaken for an augmented exact exponential.
use super::*;
use fs_couple::frame_vibration::thermoelastic::FrameThermoelasticSystem;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeModel};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection,
};
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;

const COMPANION_MASS: f64 = 0.05;
const COMPANION_K: f64 = 12.0;
const COMPANION_C: f64 = 0.08;
const SUPPORT_K: f64 = 2.0;
const SUPPORT_C: f64 = 0.02;
const SUPPORT_REST: f64 = 0.001;
const BRIDGE_K: f64 = 4.0;
const BRIDGE_C: f64 = 0.07;
const BRIDGE_REST: f64 = -0.0015;
const AIR_OMEGA: f64 = 25.0;
const AIR_ZETA: f64 = 0.08;
const AIR_FRAME: f64 = 0.6;
const AIR_COMPANION: f64 = 0.4;
const MAXWELL_K: f64 = 3.0;
const MAXWELL_TAU: f64 = 0.04;
const EXTERNAL: [f64; 2] = [0.01, 0.003]; // frame, companion [N]
const START: usize = 1;

fn assembly_config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_connections: 11,
        max_setup_terms: 32_768,
        ..config()
    }
}
fn source_specs(live: bool) -> [Spec; 2] {
    let mut s = specimens();
    s[0].modulus_slope = if live { -2.0 } else { 0.0 };
    s[1].modulus_slope = if live { -3.0 } else { 0.0 };
    s[1].curved_eta = live;
    for spec in &mut s {
        spec.cp *= 0.01;
    }
    s
}
fn initial(memory: bool) -> [f64; 7] {
    [
        0.01,
        0.03,
        -0.004,
        0.02,
        0.0,
        0.0,
        if memory { 0.002 } else { 0.0 },
    ]
}
fn attachment(component: usize, shape: f64) -> ModalAttachment {
    ModalAttachment {
        component,
        shapes: vec![shape],
    }
}
fn companion(rate: u32) -> ModalAcousticTimeModel {
    let omega = (COMPANION_K / COMPANION_MASS).sqrt();
    let mut model = ModalAcousticTimeModel::try_new(
        rate,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: COMPANION_C / (2.0 * COMPANION_MASS * omega),
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: -0.004 * COMPANION_MASS.sqrt(),
            velocity_m_sqrt_kg_per_s: 0.02 * COMPANION_MASS.sqrt(),
        }])
        .unwrap();
    model
}
fn support() -> ModalConnection {
    ModalConnection {
        left: attachment(0, 1.0 / COMPANION_MASS.sqrt()),
        right: attachment(0, 0.0),
        stiffness_n_m: SUPPORT_K,
        damping_n_s_m: SUPPORT_C,
        rest_extension_m: SUPPORT_REST,
    }
}
fn bridge(shape: f64) -> ModalConnection {
    ModalConnection {
        left: attachment(1, shape),
        right: attachment(0, 1.0 / COMPANION_MASS.sqrt()),
        stiffness_n_m: BRIDGE_K,
        damping_n_s_m: BRIDGE_C,
        rest_extension_m: BRIDGE_REST,
    }
}
fn assembled<'n, 'a>(
    component: FrameVibrationComponent,
    thermal: &DepositionSimulation<'_>,
    source: &'n MaterialNetwork<'a>,
    rate: u32,
    cached: bool,
    memory: bool,
) -> FrameThermoelasticSystem<'n, 'a> {
    let gate = CancelGate::new_clock_free();
    let mut source_component = component
        .basis()
        .coupled_component(
            rate,
            &[0.0],
            ModalAcousticTimeBudget::audible_reference(),
            18,
            || false,
        )
        .unwrap();
    source_component.restore_states(component.states()).unwrap();
    let (binding, model, ports) = source_component
        .into_thermoelastic_parts(
            thermal,
            source,
            &DESTINATIONS,
            0.05,
            1,
            START,
            assembly_config(),
            &gate,
        )
        .unwrap();
    let shape = binding.basis().modes()[0].displacement[3];
    let mut connections = vec![support()];
    connections.extend(ports);
    connections.push(bridge(shape));
    let models = vec![companion(rate), model];
    let mut mechanical = if cached {
        CoupledModalSystem::new_prepared(models, connections, assembly_config(), prepared(), &gate)
    } else {
        CoupledModalSystem::new(models, connections, assembly_config(), &gate)
    }
    .unwrap();
    let air = Model {
        ports: 2,
        poles: vec![Pole {
            omega: AIR_OMEGA,
            zeta: AIR_ZETA,
            coupling: vec![AIR_COMPANION / COMPANION_MASS.sqrt(), AIR_FRAME * shape],
        }],
    };
    mechanical = mechanical
        .with_radiation_load(&air, &[0, 1], 1, &gate)
        .unwrap();
    if memory {
        mechanical = mechanical
            .with_material_relaxation(
                vec![RelaxationBranch {
                    projection: vec![-0.4 / COMPANION_MASS.sqrt(), 0.0, shape, 0.0],
                    stiffness: MAXWELL_K,
                    relaxation_time_s: MAXWELL_TAU,
                }],
                Some(vec![0.002]),
                &gate,
            )
            .unwrap();
    }
    let before = mechanical.total_energy_j().unwrap();
    let owner = binding.into_system(mechanical, thermal, &gate).unwrap();
    assert_eq!(owner.component_index(), 1);
    assert_eq!(owner.connection_start(), START);
    assert_eq!(
        owner.mechanical().total_energy_j().unwrap().to_bits(),
        before.to_bits()
    );
    owner
}
fn physical(owner: &FrameThermoelasticSystem<'_, '_>) -> [f64; 7] {
    let shape = owner.basis().modes()[0].displacement[3];
    let m = owner.mechanical();
    let frame = m.components()[1].states()[0];
    let companion = m.components()[0].states()[0];
    [
        shape * frame.displacement_m_sqrt_kg,
        shape * frame.velocity_m_sqrt_kg_per_s,
        companion.displacement_m_sqrt_kg / COMPANION_MASS.sqrt(),
        companion.velocity_m_sqrt_kg_per_s / COMPANION_MASS.sqrt(),
        m.radiation_state()[0],
        m.radiation_state()[1],
        m.material_memory().first().copied().unwrap_or(0.0),
    ]
}
fn forces(owner: &FrameThermoelasticSystem<'_, '_>) -> [f64; 2] {
    [
        EXTERNAL[1] / COMPANION_MASS.sqrt(),
        EXTERNAL[0] * owner.basis().modes()[0].displacement[3],
    ]
}
fn internal_energy(x: [f64; 7], specs: [Spec; 2], memory: bool) -> f64 {
    let mass = RHO * areas().iter().sum::<f64>() * LENGTH / 3.0;
    let source = (0..2)
        .map(|i| 0.5 * (YOUNG - 300.0 * specs[i].modulus_slope) * areas()[i] / LENGTH * x[0] * x[0])
        .sum::<f64>();
    0.5 * mass * x[1] * x[1]
        + source
        + 0.5 * COMPANION_MASS * x[3] * x[3]
        + 0.5 * COMPANION_K * x[2] * x[2]
        + 0.5 * SUPPORT_K * (x[2] - SUPPORT_REST).powi(2)
        + 0.5 * BRIDGE_K * (x[0] - x[2] - BRIDGE_REST).powi(2)
        + 0.5 * ((AIR_OMEGA * x[4]).powi(2) + x[5] * x[5])
        + if memory {
            0.5 * (MAXWELL_K.sqrt() * (x[0] - 0.4 * x[2]) - x[6]).powi(2)
        } else {
            0.0
        }
}
fn caloric(thermal: &DepositionSimulation<'_>) -> f64 {
    thermal
        .cells()
        .values()
        .map(|c| c.mass_kg * c.specific_h_j_kg)
        .sum()
}

// Independent physical coordinates and consistent frame inertia, with no
// native projected K/C, transition matrices, ledgers or endpoint energy read.
// Last six coordinates integrate source loss, foreign connector loss, other
// component loss, air loss, Maxwell loss, and supplied physical work.
fn reference(x: [f64; 7], t: [f64; 2], specs: [Spec; 2], dt: f64, memory: bool) -> [f64; 13] {
    let mass = RHO * areas().iter().sum::<f64>() * LENGTH / 3.0;
    let stiffness = (0..2)
        .map(|i| (YOUNG + specs[i].modulus_slope * (t[i] - 300.0)) * areas()[i] / LENGTH)
        .sum::<f64>();
    let damping = (0..2)
        .map(|i| {
            specs[i].axial * areas()[i] / LENGTH
                * if specs[i].curved_eta {
                    0.005 * t[i] - 0.5
                } else {
                    1.0
                }
        })
        .sum::<f64>();
    let derivative = |x: [f64; 13]| -> [f64; 13] {
        let mut d = [0.0; 13];
        let bridge = BRIDGE_K * (x[0] - x[2] - BRIDGE_REST) + BRIDGE_C * (x[1] - x[3]);
        let mismatch = if memory {
            MAXWELL_K.sqrt() * (x[0] - 0.4 * x[2]) - x[6]
        } else {
            0.0
        };
        let stress = if memory {
            MAXWELL_K.sqrt() * mismatch
        } else {
            0.0
        };
        d[0] = x[1];
        d[1] = (-stiffness * x[0] - damping * x[1] - bridge - AIR_FRAME * x[5] - stress
            + EXTERNAL[0])
            / mass;
        d[2] = x[3];
        d[3] = (-COMPANION_K * x[2]
            - COMPANION_C * x[3]
            - SUPPORT_K * (x[2] - SUPPORT_REST)
            - SUPPORT_C * x[3]
            + bridge
            - AIR_COMPANION * x[5]
            + 0.4 * stress
            + EXTERNAL[1])
            / COMPANION_MASS;
        d[4] = x[5];
        d[5] = -AIR_OMEGA * AIR_OMEGA * x[4] - 2.0 * AIR_ZETA * AIR_OMEGA * x[5]
            + AIR_FRAME * x[1]
            + AIR_COMPANION * x[3];
        d[6] = mismatch / MAXWELL_TAU;
        d[7] = damping * x[1] * x[1];
        d[8] = SUPPORT_C * x[3] * x[3] + BRIDGE_C * (x[1] - x[3]).powi(2);
        d[9] = COMPANION_C * x[3] * x[3];
        d[10] = 2.0 * AIR_ZETA * AIR_OMEGA * x[5] * x[5];
        d[11] = mismatch * mismatch / MAXWELL_TAU;
        d[12] = EXTERNAL[0] * x[1] + EXTERNAL[1] * x[3];
        d
    };
    let mut state = [0.0; 13];
    state[..7].copy_from_slice(&x);
    let h = dt / 64.0;
    for _ in 0..64 {
        let a = derivative(state);
        let b = derivative(std::array::from_fn(|i| state[i] + 0.5 * h * a[i]));
        let c = derivative(std::array::from_fn(|i| state[i] + 0.5 * h * b[i]));
        let d = derivative(std::array::from_fn(|i| state[i] + h * c[i]));
        state =
            std::array::from_fn(|i| state[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    state
}

#[test]
fn live_source_frame_drives_other_body_and_reciprocal_air_without_restarting_any_owner() {
    let gate = CancelGate::new_clock_free();
    let mut endpoints = Vec::new();
    for live in [false, true] {
        let mut prepared_errors = Vec::new();
        for (cached, rate, memory) in [
            (false, RATE, false),
            (false, RATE, true),
            (true, RATE, false),
            (true, 2 * RATE, false),
        ] {
            let specs = source_specs(live);
            let endpoint = fixture(
                specs,
                false,
                false,
                0.0,
                |component, mut thermal, source| {
                    let mut owner = assembled(component, &thermal, source, rate, cached, memory);
                    let schema = owner.mechanical().port_schemas().to_vec();
                    let unchanged = [
                        owner.mechanical().connection_coefficients(0).unwrap(),
                        owner.mechanical().connection_coefficients(10).unwrap(),
                    ];
                    let mut x = initial(memory);
                    let mut temperature = [300.0; 2];
                    let mut maximum_error = 0.0_f64;
                    let mut observed_external_loss = 0.0;
                    let mut total = internal_energy(x, specs, memory) + caloric(&thermal);
                    for sample in 0..20 * rate / RATE {
                        let next = reference(x, temperature, specs, 1.0 / f64::from(rate), memory);
                        let input = forces(&owner);
                        let report = owner
                            .step(&mut thermal, &input, solve(), &[0.0; 2], &gate)
                            .unwrap();
                        let actual = physical(&owner);
                        for i in 0..7 {
                            let error = (actual[i] - next[i]).abs();
                            maximum_error = maximum_error.max(error);
                            if !cached {
                                near(actual[i], next[i], 2e-8);
                            }
                        }
                        let c: [f64; 2] = std::array::from_fn(|i| {
                            specs[i].axial * areas()[i] / LENGTH
                                * if specs[i].curved_eta {
                                    0.005 * temperature[i] - 0.5
                                } else {
                                    1.0
                                }
                        });
                        for i in 0..2 {
                            let heat = next[7] * c[i] / c.iter().sum::<f64>();
                            let reversible =
                                0.5 * temperature[i] * specs[i].modulus_slope * areas()[i] / LENGTH
                                    * (next[0] * next[0] - x[0] * x[0]);
                            temperature[i] +=
                                (heat + reversible) / (RHO * areas()[i] * LENGTH * specs[i].cp);
                            if !cached {
                                near(report.layers[i].material_heat_j, heat, 2e-11);
                                near(report.layers[i].reversible_heat_j, reversible, 2e-11);
                                near(
                                    thermal.temperature_k(CELLS[i]).unwrap(),
                                    temperature[i],
                                    2e-7,
                                );
                            }
                        }
                        if !cached {
                            near(report.deposited_material_heat_j, next[7], 2e-11);
                            near(
                                report.external_mechanical_loss_j,
                                next[8] + next[9] + next[10] + next[11],
                                2e-11,
                            );
                            near(report.mechanical.external_work_j, next[12], 2e-11);
                            near(report.mechanical.radiation_dissipation_j, next[10], 2e-11);
                            near(report.mechanical.material_dissipation_j, next[11], 2e-11);
                        }
                        observed_external_loss += report.external_mechanical_loss_j;
                        total +=
                            report.mechanical.external_work_j - report.external_mechanical_loss_j;
                        near(
                            owner.mechanical_internal_energy_j().unwrap(),
                            internal_energy(actual, specs, memory),
                            2e-10,
                        );
                        near(
                            internal_energy(actual, specs, memory) + caloric(&thermal),
                            total,
                            2e-9,
                        );
                        assert_eq!(owner.mechanical().samples_rendered(), u64::from(sample + 1));
                        assert_eq!(owner.mechanical().port_schemas(), schema.as_slice());
                        assert_eq!(
                            owner.mechanical().connection_coefficients(0).unwrap(),
                            unchanged[0]
                        );
                        assert_eq!(
                            owner.mechanical().connection_coefficients(10).unwrap(),
                            unchanged[1]
                        );
                        near(
                            report.total_energy_residual_j,
                            0.0,
                            solve().total_energy_tolerance_j,
                        );
                        near(report.entropy.accounting_residual_j_k, 0.0, 2e-12);
                        x.copy_from_slice(&next[..7]);
                    }
                    assert!(observed_external_loss > 1e-5);
                    assert!(owner.mechanical().radiation_energy_j().unwrap() > 0.0);
                    if memory {
                        assert_ne!(owner.mechanical().material_memory()[0], 0.002);
                    }
                    if cached {
                        prepared_errors.push(maximum_error);
                    }
                    physical(&owner)
                },
            );
            if !cached && !memory {
                endpoints.push(endpoint);
            }
        }
        assert!(prepared_errors[0] < 2e-5);
        assert!(
            prepared_errors[1] < 0.4 * prepared_errors[0],
            "{prepared_errors:?}"
        );
    }
    // Only source data changed. Its effect crosses the reciprocal bridge and
    // reaches the original air state as well as the source frame motion.
    assert!((endpoints[0][0] - endpoints[1][0]).abs() > 1e-7);
    assert!((endpoints[0][2] - endpoints[1][2]).abs() > 1e-8);
    assert!((endpoints[0][5] - endpoints[1][5]).abs() > 1e-8);
}

#[test]
fn assembled_heat_refusal_and_cancellation_preserve_all_histories_for_clean_retry() {
    fixture(
        source_specs(true),
        false,
        false,
        0.0,
        |component, mut thermal, source| {
            let gate = CancelGate::new_clock_free();
            let basis = component.basis().clone();
            let states = component.states().to_vec();
            let mut other = basis
                .coupled_component(
                    RATE,
                    &[0.0],
                    ModalAcousticTimeBudget::audible_reference(),
                    18,
                    || false,
                )
                .unwrap();
            other.restore_states(&states).unwrap();
            let mut clean_thermal = thermal.clone();
            let mut owner = assembled(component, &thermal, source, RATE, false, true);
            let mut clean = assembled(other, &clean_thermal, source, RATE, false, true);
            let input = forces(&owner);
            for _ in 0..4 {
                owner
                    .step(&mut thermal, &input, solve(), &[0.0; 2], &gate)
                    .unwrap();
                clean
                    .step(&mut clean_thermal, &input, solve(), &[0.0; 2], &gate)
                    .unwrap();
            }
            let x = physical(&owner);
            let frame = owner.mechanical().last_frame().cloned();
            let transfers = owner.mechanical().last_window_transfers().unwrap().to_vec();
            let coefficients: Vec<_> = (0..11)
                .map(|i| owner.mechanical().connection_coefficients(i).unwrap())
                .collect();
            let caloric = thermal.cells().clone();
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            for controls in [
                FrameThermalSolve {
                    maximum_material_trials: 1,
                    ..solve()
                },
                solve(),
            ] {
                let selected = if controls.maximum_material_trials == 1 {
                    &gate
                } else {
                    &cancelled
                };
                assert!(
                    owner
                        .step(&mut thermal, &input, controls, &[0.0; 2], selected)
                        .is_err()
                );
                assert_eq!(physical(&owner), x);
                assert_eq!(owner.mechanical().last_frame(), frame.as_ref());
                assert_eq!(
                    owner.mechanical().last_window_transfers().unwrap(),
                    transfers.as_slice()
                );
                assert_eq!(
                    (0..11)
                        .map(|i| owner.mechanical().connection_coefficients(i).unwrap())
                        .collect::<Vec<_>>(),
                    coefficients
                );
                assert_eq!(thermal.cells(), &caloric);
                assert_eq!(thermal.time_s().to_bits(), clean_thermal.time_s().to_bits());
            }
            for _ in 0..8 {
                owner
                    .step(&mut thermal, &input, solve(), &[0.0; 2], &gate)
                    .unwrap();
                clean
                    .step(&mut clean_thermal, &input, solve(), &[0.0; 2], &gate)
                    .unwrap();
                assert_eq!(physical(&owner), physical(&clean));
                assert_eq!(
                    owner.mechanical().last_frame(),
                    clean.mechanical().last_frame()
                );
                assert_eq!(thermal.cells(), clean_thermal.cells());
            }
        },
    );
}

#[test]
fn source_binding_refuses_foreign_columns_modified_laws_and_unobserved_warm_history() {
    for defect in 0..5 {
        fixture(
            source_specs(true),
            false,
            false,
            0.0,
            |component, thermal, source| {
                let gate = CancelGate::new_clock_free();
                let (binding, mut model, mut ports) = component
                    .into_thermoelastic_parts(
                        &thermal,
                        source,
                        &DESTINATIONS,
                        0.05,
                        1,
                        START,
                        assembly_config(),
                        &gate,
                    )
                    .unwrap();
                match defect {
                    0 => ports[0].right = attachment(0, 0.01),
                    1 => ports[3].stiffness_n_m *= 1.01,
                    2 => ports[3].rest_extension_m += 0.0001,
                    3 => model
                        .restore_states(&[ModalAcousticState::default()])
                        .unwrap(),
                    _ => (),
                }
                let mut links = vec![support()];
                links.extend(ports);
                let mut mechanical = CoupledModalSystem::new(
                    vec![companion(RATE), model],
                    links,
                    assembly_config(),
                    &gate,
                )
                .unwrap();
                if defect == 4 {
                    mechanical.step_under_gate(&[0.0; 2], &gate).unwrap();
                }
                assert!(binding.into_system(mechanical, &thermal, &gate).is_err());
                assert_eq!(thermal.time_s(), 0.0);
            },
        );
    }
}
