//! Actual connected native bodies driven by sourced event-resolved adhesion.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;
use common::{Z, close, sources};
use fs_contact::normal_patch::dmt::elastic_events::DmtElasticEventBudget;
use fs_contact::normal_patch::dmt::events::{
    DmtSlipEventBudget, DmtSlipEventKind, MAX_SLIP_EVENTS,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::adhesive_elastic_events::*;
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_ports::{ModuleWindowAudit, PortedModule};

const MA: f64 = 2e-9;
const MB: f64 = 3e-9;
const MU: f64 = MA * MB / (MA + MB);
const K: f64 = 0.12;
const DT: f64 = 1e-4;

fn config() -> AdhesiveElasticEventModalConfig {
    AdhesiveElasticEventModalConfig {
        events: DmtElasticEventBudget {
            motion: DmtSlipEventBudget {
                maximum_step_s: DT,
                minimum_indentation_m: -5.0 * Z,
                maximum_normal_increment_m: Z,
                maximum_accumulated_slip_m: 1e-3,
                maximum_tangent_speed_m_s: 1.0,
                maximum_force_n: 1e-2,
                maximum_events: MAX_SLIP_EVENTS,
                impulse_absolute_tolerance_n_s: 1e-24,
                impulse_relative_tolerance: 2e-11,
                energy_absolute_tolerance_j: 1e-26,
                energy_relative_tolerance: 2e-11,
            },
            maximum_spring_extension_m: 1e-3,
            maximum_stage_speed_m_s: 1.0,
            maximum_segments: 1024,
            maximum_root_iterations: 64,
            event_time_absolute_tolerance_s: 1e-16,
            maximum_numerical_work_uncertainty_j: 1e-22,
        },
        coordinate_absolute_tolerance_m: 2e-19,
        coordinate_relative_tolerance: 3e-11,
        velocity_absolute_tolerance_m_s: 2e-15,
        velocity_relative_tolerance: 3e-11,
    }
}
fn column() -> [f64; 2] {
    [1.0 / det::sqrt(MA), -1.0 / det::sqrt(MB)]
}
fn network(
    x: f64,
    v: f64,
    vc: f64,
    damping: f64,
    displacement_limit: f64,
    modal_elastic: bool,
    links: usize,
) -> CoupledModalSystem {
    network_with_energy_tolerance(
        x,
        v,
        vc,
        damping,
        displacement_limit,
        modal_elastic,
        links,
        1e-26,
        3e-11,
    )
}

#[allow(clippy::too_many_arguments)]
fn network_with_energy_tolerance(
    x: f64,
    v: f64,
    vc: f64,
    damping: f64,
    displacement_limit: f64,
    modal_elastic: bool,
    links: usize,
    absolute: f64,
    relative: f64,
) -> CoupledModalSystem {
    let models = [MA, MB]
        .into_iter()
        .enumerate()
        .map(|(i, mass)| {
            let factor = if i == 0 {
                MB / (MA + MB)
            } else {
                -MA / (MA + MB)
            };
            let mut budget = ModalAcousticTimeBudget::audible_reference();
            budget.maximum_abs_displacement_m_sqrt_kg = displacement_limit;
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                10_000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: if modal_elastic { 100.0 } else { 0.0 },
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                budget,
            )
            .unwrap();
            model
                .restore_states(&[ModalAcousticState {
                    displacement_m_sqrt_kg: det::sqrt(mass) * factor * x,
                    velocity_m_sqrt_kg_per_s: det::sqrt(mass) * (vc + factor * v),
                }])
                .unwrap();
            model
        })
        .collect();
    let connections = (0..links)
        .map(|_| ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0 / det::sqrt(MA)],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![1.0 / det::sqrt(MB)],
            },
            stiffness_n_m: K,
            damping_n_s_m: damping,
            rest_extension_m: 0.0,
        })
        .collect();
    CoupledModalSystem::new(
        models,
        connections,
        ModalCouplingConfig {
            max_modes: 2,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e-9,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: 1e-2,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: absolute,
            energy_relative_tolerance: relative,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn make<'a>(
    n: &'a fs_contact::normal_patch::dmt::SourceDmtSphereContact,
    s: &'a fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear,
    x: f64,
    v: f64,
    vc: f64,
    damping: f64,
) -> AdhesiveElasticEventModalSystem<'a> {
    let materials = n.normal_state().elastic().interface().surface_materials();
    AdhesiveElasticEventModalSystem::new(
        network(x, v, vc, damping, 1.0, false, 1),
        n,
        s,
        [&materials[0], &materials[1]],
        0.5 * Z,
        x,
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn states(s: &AdhesiveElasticEventModalSystem<'_>) -> Vec<ModalAcousticState> {
    s.components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn momentum(s: &AdhesiveElasticEventModalSystem<'_>) -> f64 {
    [MA, MB]
        .into_iter()
        .zip(s.components())
        .map(|(m, b)| det::sqrt(m) * b.states()[0].velocity_m_sqrt_kg_per_s)
        .sum()
}
fn center(s: &AdhesiveElasticEventModalSystem<'_>) -> f64 {
    [MA, MB]
        .into_iter()
        .zip(s.components())
        .map(|(m, b)| det::sqrt(m) * b.states()[0].displacement_m_sqrt_kg)
        .sum::<f64>()
        / (MA + MB)
}
fn check(a: f64, b: f64) {
    close(a, b, 3e-23, 2e-9);
}

#[test]
fn original_connection_and_common_motion_follow_independent_oscillator() {
    let (n, s) = sources(1e-5, 0.0, 0.0);
    let (x0, v0, vc) = (2e-7, 0.001, 0.002);
    let mut owner = make(&n, &s, x0, v0, vc, 0.0);
    let initial = owner.total_energy_j().unwrap();
    let w = det::sqrt(K / MU);
    let f = *owner.step(&[0.0; 2]).unwrap();
    let x = x0 * det::cos(w * DT) + v0 / w * det::sin(w * DT);
    let v = -w * x0 * det::sin(w * DT) + v0 * det::cos(w * DT);
    check(owner.state().tangent_displacement_m(), x);
    check(owner.state().tangent_velocity_m_s(), v);
    check(owner.effective_mass_kg(), MU);
    check(center(&owner), vc * DT);
    check(momentum(&owner), (MA + MB) * vc);
    check(f.spring_energy_j, 0.5 * K * x * x);
    check(f.connection_force_n, -K * x);
    check(owner.total_energy_j().unwrap(), initial);
    assert_eq!(f.friction_dissipation_j, 0.0);
    assert_eq!(f.damping_dissipation_j, 0.0);
    let native = owner.network().last_frame().unwrap();
    check(native.connection_energy_j, f.spring_energy_j);
    check(native.connection_forces_n[0], f.connection_force_n);
    assert_eq!(native.sample, 1);
    let ModuleWindowAudit::Available(audit) = owner.network().audit_last_window().unwrap() else {
        panic!("original port audit")
    };
    assert!(audit.violations().is_empty());
}

#[test]
fn actual_loading_spring_supports_stop_static_hold_and_loaded_restart() {
    let (n, s) = sources(1e-5, 2e7, 1e7);
    let fkin = 1e7 * n.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let mut owner = make(&n, &s, 0.0, 0.0005, 0.002, 0.0);
    let initial = owner.total_energy_j().unwrap();
    let f = *owner.step(&[0.0; 2]).unwrap();
    assert_eq!(f.path.events[0].unwrap().kind, DmtSlipEventKind::Stopped);
    assert_eq!(owner.state().tangent_velocity_m_s(), 0.0);
    check(
        owner.total_energy_j().unwrap() + f.friction_dissipation_j,
        initial,
    );
    let rest = *owner.step(&[0.0; 2]).unwrap();
    assert_eq!(rest.path.accumulated_slip_m, 0.0);
    let before = owner.total_energy_j().unwrap();
    let drive = column().map(|b| b * 3.0 * fkin);
    let restart = *owner.step(&drive).unwrap();
    assert_eq!(
        restart.path.events[0].unwrap().kind,
        DmtSlipEventKind::SlipStarted
    );
    assert!(restart.path.accumulated_slip_m > 0.0 && restart.friction_dissipation_j > 0.0);
    check(
        owner.total_energy_j().unwrap() + restart.friction_dissipation_j - restart.external_work_j,
        before,
    );
    check(momentum(&owner), (MA + MB) * 0.002);
    assert_eq!(owner.samples_rendered(), 3);
    assert_eq!(owner.state().steps(), 3);
}

#[test]
fn original_dashpot_and_driven_center_of_mass_share_actual_work() {
    let (n, s) = sources(1e-5, 0.0, 0.0);
    let (x0, v0, vc, c) = (1e-7, 0.001, 0.002, 4.8e-6);
    let (fa, fb) = (1e-8, -2e-8);
    let drive = [fa / det::sqrt(MA), fb / det::sqrt(MB)];
    let q = MU * (fa / MA - fb / MB);
    let mut owner = make(&n, &s, x0, v0, vc, c);
    let before = owner.total_energy_j().unwrap();
    let alpha = c / (2.0 * MU);
    let wd = det::sqrt(K / MU - alpha * alpha);
    let a = x0 - q / K;
    let b = (v0 + alpha * a) / wd;
    let exact = |t: f64| {
        let decay = det::exp(-alpha * t);
        let (sn, cs) = (det::sin(wd * t), det::cos(wd * t));
        (
            q / K + decay * (a * cs + b * sn),
            decay * ((-alpha * a + wd * b) * cs + (-alpha * b - wd * a) * sn),
        )
    };
    let panels = 1024;
    let h = DT / panels as f64;
    let integral = (0..=panels)
        .map(|i| {
            let weight = if i == 0 || i == panels {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            weight * exact(i as f64 * h).1.powi(2)
        })
        .sum::<f64>()
        * h
        / 3.0;
    let f = *owner.step(&drive).unwrap();
    let (x, v) = exact(DT);
    check(owner.state().tangent_displacement_m(), x);
    check(owner.state().tangent_velocity_m_s(), v);
    check(f.damping_dissipation_j, c * integral);
    assert!(f.damping_dissipation_j > 0.0);
    check(momentum(&owner), (MA + MB) * vc + (fa + fb) * DT);
    check(
        center(&owner),
        vc * DT + (fa + fb) * DT * DT / (2.0 * (MA + MB)),
    );
    check(
        owner.total_energy_j().unwrap() + f.damping_dissipation_j - f.external_work_j,
        before,
    );
    let native = owner.network().last_frame().unwrap();
    check(native.connection_dissipation_j, f.damping_dissipation_j);
    check(native.external_work_j, f.external_work_j);
    let ModuleWindowAudit::Available(audit) = owner.network().audit_last_window().unwrap() else {
        panic!("original port audit")
    };
    assert!(audit.violations().is_empty());
}

#[test]
fn unsupported_topology_and_original_source_assignment_are_refused() {
    let (n, s) = sources(1e-5, 2e7, 1e7);
    let (_, alien) = sources(2e-5, 2e7, 1e7);
    let materials = n.normal_state().elastic().interface().surface_materials();
    for (elastic, links, swap, wrong_source) in [
        (true, 1, false, false),
        (false, 0, false, false),
        (false, 2, false, false),
        (false, 1, true, false),
        (false, 1, false, true),
    ] {
        let pair = if swap {
            [&materials[1], &materials[0]]
        } else {
            [&materials[0], &materials[1]]
        };
        assert!(
            AdhesiveElasticEventModalSystem::new(
                network(0.0, 0.0, 0.0, 0.0, 1.0, elastic, links),
                &n,
                if wrong_source { &alien } else { &s },
                pair,
                0.5 * Z,
                0.0,
                config(),
                &CancelGate::new_clock_free()
            )
            .is_err()
        );
    }
}

#[test]
fn late_original_component_budget_and_cancellation_preserve_native_retry() {
    let (n, s) = sources(1e-5, 2e7, 1e7);
    let materials = n.normal_state().elastic().interface().surface_materials();
    let create = || {
        AdhesiveElasticEventModalSystem::new(
            network(0.0, 0.0, 0.0, 0.0, 1e-10, false, 1),
            &n,
            &s,
            [&materials[0], &materials[1]],
            0.5 * Z,
            0.0,
            config(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    };
    let (mut owner, mut reference) = (create(), create());
    owner.step(&[0.0; 2]).unwrap();
    reference.step(&[0.0; 2]).unwrap();
    let original = (
        states(&owner),
        owner.state(),
        *owner.last_frame().unwrap(),
        owner.network().last_frame().cloned(),
    );
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert!(matches!(
        owner.step_under_gate(&[0.0; 2], &cancel),
        Err(AdhesiveElasticEventModalError::Network(
            ModalCouplingError::Cancelled
        ))
    ));
    // Equal physical acceleration leaves contact motion zero but exceeds an
    // original body's displacement limit after the candidate has been formed.
    let common = [1000.0 * det::sqrt(MA), 1000.0 * det::sqrt(MB)];
    assert!(matches!(
        owner.step(&common),
        Err(AdhesiveElasticEventModalError::Network(
            ModalCouplingError::Component { .. }
        ))
    ));
    assert_eq!(
        (
            states(&owner),
            owner.state(),
            *owner.last_frame().unwrap(),
            owner.network().last_frame().cloned()
        ),
        original
    );
    assert_eq!(owner.samples_rendered(), 1);
    assert_eq!(
        *owner.step(&[0.0; 2]).unwrap(),
        *reference.step(&[0.0; 2]).unwrap()
    );
    assert_eq!(states(&owner), states(&reference));
}

#[test]
fn conservative_propagation_uncertainty_must_fit_original_native_energy_budget() {
    let (n, s) = sources(1e-5, 0.0, 0.0);
    let mut ordinary = make(&n, &s, 2e-7, 0.001, 0.0, 0.0);
    let frame = *ordinary.step(&[0.0; 2]).unwrap();
    assert!(frame.numerical_work_uncertainty_j > frame.energy_residual_j.abs());
    assert_eq!(frame.friction_dissipation_j, 0.0);
    assert_eq!(frame.damping_dissipation_j, 0.0);
    let materials = n.normal_state().elastic().interface().surface_materials();
    let mut strict = AdhesiveElasticEventModalSystem::new(
        network_with_energy_tolerance(
            2e-7,
            0.001,
            0.0,
            0.0,
            1.0,
            false,
            1,
            0.5 * frame.numerical_work_uncertainty_j,
            1e-30,
        ),
        &n,
        &s,
        [&materials[0], &materials[1]],
        0.5 * Z,
        2e-7,
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let original = (
        states(&strict),
        strict.state(),
        strict.network().last_frame().cloned(),
    );
    assert!(matches!(
        strict.step(&[0.0; 2]),
        Err(AdhesiveElasticEventModalError::EnergyAccuracy { .. })
    ));
    assert_eq!(
        (
            states(&strict),
            strict.state(),
            strict.network().last_frame().cloned()
        ),
        original
    );
    assert_eq!(strict.samples_rendered(), 0);
    assert!(strict.last_frame().is_none());
}
