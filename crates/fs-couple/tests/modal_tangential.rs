//! G1/G3: stateful interface work on the actual two-way mechanical network.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::friction::tangential::{
    TangentialModalDrive, TangentialModalSystem,
};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, FrictionModalSystem, ModalFrictionConfig, ModalFrictionPort,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalCouplingConfig,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::mindlin::{MindlinDeresiewiczPoint, MindlinPointParameters};
use fs_tribo::{FrictionLaw, TangentialLaw};

#[path = "modal_tangential/heat.rs"]
mod heat;

fn budgets() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 4,
        max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-8,
    }
}

fn config() -> ModalFrictionConfig {
    ModalFrictionConfig {
        max_iterations: 80,
        maximum_force_n: 100.0,
        maximum_slip_speed_m_s: 10.0,
        force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-10,
        sticking_tolerance_m: 1e-12,
    }
}

fn image(omega: f64, velocity: f64) -> ModalAcousticTimeModel {
    image_at_rate(48000, omega, velocity)
}

fn image_at_rate(rate_hz: u32, omega: f64, velocity: f64) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new(
        rate_hz,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::new(1.0, 0.0),
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: velocity,
        }])
        .unwrap();
    model
}

fn attachment(component: usize) -> ModalAttachment {
    let mass: f64 = if component == 0 { 0.002 } else { 0.01 };
    ModalAttachment {
        component,
        shapes: vec![mass.sqrt().recip()],
    }
}

fn network(velocity: f64, budget: ModalCouplingConfig) -> CoupledModalSystem {
    CoupledModalSystem::new(
        vec![image(800.0, velocity), image(1100.0, 0.0)],
        vec![],
        budget,
        &CancelGate::new(),
    )
    .unwrap()
}

fn layer(yield_force_n: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new(2e5, yield_force_n, 0.1, 300.0, 0.01).unwrap()
}

fn system<L: TangentialLaw>(law: L) -> TangentialModalSystem<L> {
    TangentialModalSystem::new(
        network(0.0, budgets()),
        attachment(0),
        Some(attachment(1)),
        law,
        config(),
        &CancelGate::new(),
    )
    .unwrap()
}

fn drive(speed_m_s: f64) -> TangentialModalDrive {
    TangentialModalDrive {
        speed_m_s,
        normal_force_n: 1.0,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
    }
}

fn states<L: TangentialLaw>(s: &TangentialModalSystem<L>) -> Vec<ModalAcousticState> {
    s.network()
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}

#[test]
fn g1_elastic_contact_uses_both_body_responses_and_stores_work() {
    let mut s = system(layer(10.0));
    let dt = s.network().sample_period_s();
    let speed = 0.1;
    // Independent forced-oscillator compliance, summed for the two reactions.
    let mut compliance = 0.0;
    for (i, omega) in [800.0_f64, 1100.0].into_iter().enumerate() {
        let b = attachment(i).shapes[0];
        compliance += b * b * (1.0 - (omega * dt).cos()) / (omega * omega);
    }
    let expected_force = 1e5 * speed * dt / (1.0 + 1e5 * compliance);
    let f = *s.step(&[0.0; 2], drive(speed)).unwrap();
    assert!((f.traction_n - expected_force).abs() < 1e-12);
    assert!(s.network().components()[0].states()[0].displacement_m_sqrt_kg > 0.0);
    assert!(s.network().components()[1].states()[0].displacement_m_sqrt_kg < 0.0);
    let elastic = s.state().elastic_slip_m()[0];
    assert!((elastic - f.slip_distance_m).abs() < 1e-14);
    assert!((f.contact_energy_j - 1e5 * elastic * elastic).abs() < 1e-18);
    assert!(f.contact_energy_j > 0.0);
    assert_eq!(f.friction_dissipation_j, 0.0);
    assert!((s.total_energy_j().unwrap() - f.drive_work_j).abs() < 1e-14);
}

#[test]
fn g1_contact_elastic_limit_converges_to_independent_coupled_frequency() {
    let mass = 0.002;
    let stiffness = 2e5;
    let initial_velocity = 0.01;
    let omega: f64 = 800.0;
    let coupled_omega = (omega * omega + stiffness / mass).sqrt();
    let duration = 0.001;
    let mut errors = Vec::new();
    for rate in [24000, 48000, 96000] {
        let network = CoupledModalSystem::new(
            vec![image_at_rate(rate, omega, initial_velocity)],
            vec![],
            budgets(),
            &CancelGate::new(),
        )
        .unwrap();
        let mut s = TangentialModalSystem::new(
            network,
            attachment(0),
            None,
            layer(10.0),
            config(),
            &CancelGate::new(),
        )
        .unwrap();
        for _ in 0..rate / 1000 {
            let f = *s.step(&[0.0], drive(0.0)).unwrap();
            assert_eq!(f.friction_dissipation_j, 0.0);
        }
        let state = s.network().components()[0].states()[0];
        // This oscillator's exact frequency includes the physical interface
        // stiffness; the stepper's original bare frequency remains 800 rad/s.
        let phase = coupled_omega * duration;
        let error = (coupled_omega * state.displacement_m_sqrt_kg - initial_velocity * phase.sin())
            .hypot(state.velocity_m_sqrt_kg_per_s - initial_velocity * phase.cos());
        errors.push(error);
        assert!((s.total_energy_j().unwrap() - 0.5 * initial_velocity.powi(2)).abs() < 1e-12);
    }
    for adjacent in errors.windows(2) {
        let order = (adjacent[0] / adjacent[1]).log2();
        assert!(
            (1.9..2.1).contains(&order),
            "elastic contact time order {order}: {errors:?}"
        );
    }
    assert!(errors[2] < 1e-4, "finest state error: {errors:?}");
}

#[test]
fn g1_reversing_viscoplastic_contact_closes_storage_and_physical_heat() {
    let mut s = system(layer(0.03));
    let mut work = 0.0;
    let mut heat = 0.0;
    let mut released_storage = false;
    for k in 0..400 {
        let f = *s
            .step(&[0.02, -0.03], drive(if k < 200 { 0.3 } else { -0.3 }))
            .unwrap();
        work += f.external_work_j + f.drive_work_j;
        heat += f.friction_dissipation_j + f.network_dissipation_j;
        released_storage |= f.contact_stored_change_j < -1e-10;
        assert!(f.friction_dissipation_j >= 0.0);
        assert!(f.contact_work_residual_j.abs() <= f.energy_tolerance_j);
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
    }
    assert!(
        released_storage,
        "reversal must retain and release elastic contact history"
    );
    assert!(heat > 1e-6);
    assert!((s.total_energy_j().unwrap() + heat - work).abs() < 1e-11);
}

#[test]
fn g3_coulomb_trait_consumer_reproduces_original_dry_coupling() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.8,
        kinetic_mu: 0.4,
    };
    let mut general = system(law.clone());
    let mut original = FrictionModalSystem::new(
        network(0.0, budgets()),
        ModalFrictionPort {
            left: attachment(0),
            right: Some(attachment(1)),
            law,
        },
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    for k in 0..128 {
        let speed = if k < 64 { 0.3 } else { -0.3 };
        let f = *general.step(&[0.0; 2], drive(speed)).unwrap();
        let reference = *original
            .step(
                &[0.0; 2],
                FrictionDrive {
                    speed_m_s: speed,
                    normal_force_n: 1.0,
                },
            )
            .unwrap();
        assert_eq!(f.traction_n.to_bits(), reference.traction_n.to_bits());
        // The trait integrates its solved slip; the original dry adapter uses
        // the independently propagated body displacement. Equivalent affine
        // response expressions can differ by roundoff, unlike the shared
        // traction and accepted body states checked bitwise below.
        assert!(
            (f.friction_dissipation_j - reference.friction_dissipation_j).abs() < 1e-15
        );
        assert_eq!(f.contact_energy_j, 0.0);
        for (a, b) in general
            .network()
            .components()
            .iter()
            .zip(original.components())
        {
            assert_eq!(a.states(), b.states());
        }
    }
}

#[test]
fn g3_cancelled_and_unsupported_steps_preserve_loaded_history_and_retry() {
    let mut s = system(layer(0.03));
    let mut direct = system(layer(0.03));
    for _ in 0..12 {
        s.step(&[0.0; 2], drive(0.1)).unwrap();
        direct.step(&[0.0; 2], drive(0.1)).unwrap();
    }
    let before = states(&s);
    let state = s.state().clone();
    let frame = s.last_frame().copied();
    let gate = CancelGate::new();
    gate.request();
    assert!(s.step_under_gate(&[0.0; 2], drive(-0.1), &gate).is_err());
    for bad in [
        TangentialModalDrive {
            contact_temperature_k: Some(301.0),
            ..drive(-0.1)
        },
        TangentialModalDrive {
            normal_force_n: 0.0,
            ..drive(-0.1)
        },
        TangentialModalDrive {
            contact_temperature_k: None,
            ..drive(-0.1)
        },
    ] {
        assert!(s.step(&[0.0; 2], bad).is_err());
        assert_eq!(s.state(), &state);
        assert_eq!(states(&s), before);
        assert_eq!(s.last_frame().copied(), frame);
    }
    s.step(&[0.0; 2], drive(-0.1)).unwrap();
    direct.step(&[0.0; 2], drive(-0.1)).unwrap();
    assert_eq!(s.state(), direct.state());
    assert_eq!(states(&s), states(&direct));
    assert_eq!(s.last_frame(), direct.last_frame());
}

#[test]
fn g0_total_energy_budget_includes_reversible_contact_storage() {
    let mut reference = system(layer(10.0));
    let f = *reference.step(&[0.0; 2], drive(0.1)).unwrap();
    let mut budget = budgets();
    budget.maximum_total_energy_j = f.network_energy_j + 0.5 * f.contact_energy_j;
    let mut bounded = TangentialModalSystem::new(
        network(0.0, budget),
        attachment(0),
        Some(attachment(1)),
        layer(10.0),
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let before = states(&bounded);
    let state = bounded.state().clone();
    assert!(bounded.step(&[0.0; 2], drive(0.1)).is_err());
    assert_eq!(states(&bounded), before);
    assert_eq!(bounded.state(), &state);
    assert!(bounded.last_frame().is_none());
    assert_eq!(bounded.network().samples_rendered(), 0);
}

fn mindlin() -> MindlinDeresiewiczPoint {
    MindlinDeresiewiczPoint::new(MindlinPointParameters {
        normal_force_n: 10.0,
        contact_radius_m: 1e-3,
        shear_moduli_pa: [30e9; 2],
        poisson_ratios: [0.25; 2],
        friction_coefficient: 0.6,
        maximum_history_intervals: 128,
    })
    .unwrap()
}

#[test]
fn g1_mindlin_actual_contact_motion_matches_hertz_backbone_and_work_effort() {
    let mut s = system(mindlin());
    let dt = s.network().sample_period_s();
    // Independently derived from the two bodies and Hertz contact radius.
    let k0 = 8.0 * 1e-3 / (2.0 * (2.0 - 0.25) / 30e9);
    let capacity = 0.6 * 10.0;
    let gross_displacement = 1.5 * capacity / k0;
    let fraction: f64 = 0.35;
    let displacement = fraction * gross_displacement;
    let endpoint_force = capacity * (1.0 - (1.0 - fraction).powf(1.5));
    let average_force =
        capacity * (1.0 - (2.0 / (5.0 * fraction)) * (1.0 - (1.0 - fraction).powf(2.5)));
    let mut compliance = 0.0;
    for (i, omega) in [800.0_f64, 1100.0].into_iter().enumerate() {
        let b = attachment(i).shapes[0];
        compliance += b * b * (1.0 - (omega * dt).cos()) / (omega * omega);
    }
    let speed = (displacement + compliance * average_force) / dt;
    let f = *s
        .step(
            &[0.0; 2],
            TangentialModalDrive {
                normal_force_n: 10.0,
                ..drive(speed)
            },
        )
        .unwrap();
    assert!((f.slip_distance_m - displacement).abs() < 1e-14);
    assert!((s.state().slip_m() - displacement).abs() < 1e-14);
    assert!((s.state().endpoint_force_n() - endpoint_force).abs() < 1e-8);
    assert!((f.traction_n - average_force).abs() < 1e-8);
    assert!(
        f.traction_n < 0.75 * s.state().endpoint_force_n(),
        "the held mechanical effort must be the path average, not endpoint force"
    );
    assert!(f.contact_energy_j > 0.0);
    assert!(
        f.friction_dissipation_j > 0.0,
        "partial slip dissipates before gross sliding"
    );
    assert!(
        (f.contact_energy_j + f.friction_dissipation_j - average_force * displacement).abs()
            < 1e-14
    );
    assert!(
        (s.total_energy_j().unwrap() + f.friction_dissipation_j - f.drive_work_j).abs() < 1e-12
    );
}

#[test]
fn g3_mindlin_loaded_reversal_keeps_history_and_refuses_normal_load_changes() {
    let mut s = system(mindlin());
    let mut direct = system(mindlin());
    let input = TangentialModalDrive {
        normal_force_n: 10.0,
        ..drive(0.04)
    };
    let mut work = 0.0;
    let mut heat = 0.0;
    for _ in 0..16 {
        let f = *s.step(&[0.0; 2], input).unwrap();
        direct.step(&[0.0; 2], input).unwrap();
        work += f.drive_work_j;
        heat += f.friction_dissipation_j + f.network_dissipation_j;
    }
    let before = states(&s);
    let state = s.state().clone();
    let last = s.last_frame().copied();
    assert!(
        s.step(
            &[0.0; 2],
            TangentialModalDrive {
                normal_force_n: 11.0,
                ..input
            }
        )
        .is_err()
    );
    assert_eq!(states(&s), before);
    assert_eq!(s.state(), &state);
    assert_eq!(s.last_frame().copied(), last);
    let mut returned_storage = false;
    for _ in 0..32 {
        let reverse = TangentialModalDrive {
            speed_m_s: -0.04,
            ..input
        };
        let f = *s.step(&[0.0; 2], reverse).unwrap();
        direct.step(&[0.0; 2], reverse).unwrap();
        work += f.drive_work_j;
        heat += f.friction_dissipation_j + f.network_dissipation_j;
        returned_storage |= f.contact_stored_change_j < -1e-10;
        assert!(f.friction_dissipation_j >= 0.0);
    }
    assert!(returned_storage);
    assert!(heat > 0.0);
    assert!((s.total_energy_j().unwrap() + heat - work).abs() < 1e-10);
    assert_eq!(s.state(), direct.state());
    assert_eq!(states(&s), states(&direct));
    assert_eq!(s.last_frame(), direct.last_frame());
}
