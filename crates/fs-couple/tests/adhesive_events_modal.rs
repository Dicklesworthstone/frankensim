//! Exact junction events applied to both original bodies, their work and clocks.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;

use common::{Z, close, sources};
use fs_contact::normal_patch::dmt::events::{
    DmtSlipEventBudget, DmtSlipEventKind, DmtSlipRegime, MAX_SLIP_EVENTS,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::adhesive_events::*;
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_ports::{ModuleWindowAudit, PortedModule};

const MA: f64 = 2e-9;
const MB: f64 = 3e-9;
const MU: f64 = MA * MB / (MA + MB);
const DT: f64 = 1e-4;

fn config() -> AdhesiveEventModalConfig {
    AdhesiveEventModalConfig {
        events: DmtSlipEventBudget {
            maximum_step_s: DT,
            minimum_indentation_m: -5.0 * Z,
            maximum_normal_increment_m: 3.0 * Z,
            maximum_accumulated_slip_m: 1e-3,
            maximum_tangent_speed_m_s: 1.0,
            maximum_force_n: 1e-2,
            maximum_events: MAX_SLIP_EVENTS,
            impulse_absolute_tolerance_n_s: 1e-25,
            impulse_relative_tolerance: 2e-12,
            energy_absolute_tolerance_j: 1e-27,
            energy_relative_tolerance: 2e-12,
        },
        coordinate_absolute_tolerance_m: 2e-20,
        coordinate_relative_tolerance: 3e-12,
        velocity_absolute_tolerance_m_s: 2e-16,
        velocity_relative_tolerance: 3e-12,
    }
}
fn column() -> Vec<f64> {
    vec![-1.0 / det::sqrt(MA), 0.0, 1.0 / det::sqrt(MB), 0.0]
}
fn network(
    vrel: f64,
    vc: f64,
    displacement_limit: f64,
    elastic: bool,
    connected: bool,
) -> CoupledModalSystem {
    let models = [MA, MB]
        .into_iter()
        .enumerate()
        .map(|(i, mass)| {
            let mut budget = ModalAcousticTimeBudget::audible_reference();
            budget.maximum_abs_displacement_m_sqrt_kg = displacement_limit;
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                10_000,
                vec![
                    ModalAcousticMode {
                        angular_frequency_rad_s: if elastic { 100.0 } else { 0.0 },
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    };
                    2
                ],
                budget,
            )
            .unwrap();
            let v = vc
                + if i == 0 {
                    -MB / (MA + MB) * vrel
                } else {
                    MA / (MA + MB) * vrel
                };
            model
                .restore_states(&[
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 0.0,
                        velocity_m_sqrt_kg_per_s: det::sqrt(mass) * v,
                    },
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 0.0,
                        velocity_m_sqrt_kg_per_s: det::sqrt(mass) * (i as f64 + 1.0) * 0.02,
                    },
                ])
                .unwrap();
            model
        })
        .collect();
    let links = if connected {
        vec![ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0, 0.0],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![1.0, 0.0],
            },
            stiffness_n_m: 0.0,
            damping_n_s_m: 0.0,
            rest_extension_m: 0.0,
        }]
    } else {
        vec![]
    };
    CoupledModalSystem::new(
        models,
        links,
        ModalCouplingConfig {
            max_modes: 4,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e-9,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: 1.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-26,
            energy_relative_tolerance: 3e-12,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn make<'a>(
    n: &'a fs_contact::normal_patch::dmt::SourceDmtSphereContact,
    s: &'a fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear,
    vrel: f64,
    vc: f64,
) -> AdhesiveEventModalSystem<'a> {
    let materials = n.normal_state().elastic().interface().surface_materials();
    AdhesiveEventModalSystem::new(
        network(vrel, vc, 1.0, false, false),
        n,
        s,
        [&materials[0], &materials[1]],
        0.5 * Z,
        0.0,
        column(),
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn states(s: &AdhesiveEventModalSystem<'_>) -> Vec<ModalAcousticState> {
    s.components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn momentum(s: &AdhesiveEventModalSystem<'_>) -> f64 {
    det::sqrt(MA) * s.components()[0].states()[0].velocity_m_sqrt_kg_per_s
        + det::sqrt(MB) * s.components()[1].states()[0].velocity_m_sqrt_kg_per_s
}
fn center(s: &AdhesiveEventModalSystem<'_>) -> f64 {
    (det::sqrt(MA) * s.components()[0].states()[0].displacement_m_sqrt_kg
        + det::sqrt(MB) * s.components()[1].states()[0].displacement_m_sqrt_kg)
        / (MA + MB)
}
fn relative_force(q: f64) -> Vec<f64> {
    column().into_iter().map(|b| b * q).collect()
}
fn check(a: f64, b: f64) {
    close(a, b, 2e-24, 1e-10);
}

#[test]
fn native_exact_stop_rest_and_restart_preserve_original_bodies_and_common_motion() {
    let (n, s) = sources(1e-5, 2e7, 1e7);
    let k = 1e7 * n.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let stop = DT / 4.0;
    let v0 = k * stop / MU;
    let mut owner = make(&n, &s, v0, 0.01);
    let before = owner.total_energy_j().unwrap();
    let p0 = momentum(&owner);
    let original = states(&owner);
    let frame = *owner.step(&[0.0; 4], 0.5 * Z).unwrap();
    assert_eq!(frame.path.event_count, 1);
    assert_eq!(
        frame.path.events[0].unwrap().kind,
        DmtSlipEventKind::Stopped
    );
    check(frame.path.events[0].unwrap().offset_s, stop);
    check(frame.path.tangent_increment_m, v0 * stop / 2.0);
    check(frame.path.friction_time_moment_n_s2, k * stop * stop / 2.0);
    assert_eq!(owner.state().regime(), DmtSlipRegime::Sticking);
    check(frame.friction_dissipation_j, MU * v0 * v0 / 2.0);
    check(
        owner.total_energy_j().unwrap() + frame.friction_dissipation_j,
        before,
    );
    check(momentum(&owner), p0);
    check(center(&owner), 0.01 * DT);
    for i in 0..4 {
        let b = column()[i];
        let actual = states(&owner)[i];
        check(
            actual.velocity_m_sqrt_kg_per_s,
            original[i].velocity_m_sqrt_kg_per_s - b * k * stop,
        );
        check(
            actual.displacement_m_sqrt_kg,
            DT * original[i].velocity_m_sqrt_kg_per_s - b * (DT * k * stop - k * stop * stop / 2.0),
        );
    }
    let rest = *owner.step(&[0.0; 4], 0.5 * Z).unwrap();
    assert_eq!(rest.path.accumulated_slip_m, 0.0);
    assert_eq!(rest.path.event_count, 0);
    let restart = *owner.step(&relative_force(3.0 * k), 0.5 * Z).unwrap();
    assert_eq!(
        restart.path.events[0].unwrap().kind,
        DmtSlipEventKind::SlipStarted
    );
    check(restart.path.events[0].unwrap().offset_s, 0.0);
    check(owner.state().tangent_velocity_m_s(), 2.0 * k * DT / MU);
    assert_eq!(owner.samples_rendered(), 3);
    assert_eq!(owner.state().steps(), 3);
    check(owner.state().time_s(), 3.0 * DT);
    assert_eq!(owner.network().last_frame().unwrap().sample, 3);
    assert_eq!(restart.observer_pressure_pa, 0.0);
    assert!(std::ptr::eq(owner.normal(), &n) && std::ptr::eq(owner.shear(), &s));
}

#[test]
fn native_zero_net_slip_reversal_keeps_actual_loss_impulses_and_positions() {
    let (n, s) = sources(1e-5, 1.5e7, 1e7);
    let k = 1e7 * n.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let stop = DT / (1.0 + det::sqrt(3.0));
    let v0 = 3.0 * k * stop / MU;
    let mut owner = make(&n, &s, v0, 0.0);
    let before = owner.total_energy_j().unwrap();
    let initial = states(&owner);
    let drive = relative_force(-2.0 * k);
    let f = *owner.step(&drive, 0.5 * Z).unwrap();
    close(f.path.tangent_increment_m, 0.0, 3e-22, 0.0);
    check(f.path.accumulated_slip_m, v0 * stop);
    check(f.friction_dissipation_j, k * v0 * stop);
    assert!(f.friction_dissipation_j > 0.0);
    assert_eq!(f.path.event_count, 2);
    check(owner.state().tangent_velocity_m_s(), -k * (DT - stop) / MU);
    check(momentum(&owner), 0.0);
    check(
        owner.total_energy_j().unwrap() + f.friction_dissipation_j - f.external_work_j,
        before,
    );
    for i in [0, 2] {
        let a = states(&owner)[i];
        check(
            a.displacement_m_sqrt_kg,
            DT * initial[i].velocity_m_sqrt_kg_per_s + DT * DT * drive[i] / 2.0
                - column()[i]
                    * (DT * f.path.friction_impulse_n_s - f.path.friction_time_moment_n_s2),
        );
        check(
            a.velocity_m_sqrt_kg_per_s,
            initial[i].velocity_m_sqrt_kg_per_s + DT * drive[i]
                - column()[i] * f.path.friction_impulse_n_s,
        );
    }
    assert_eq!(owner.state().regime(), DmtSlipRegime::Sliding);
}

#[test]
fn native_driven_common_motion_and_original_boundary_work_exclude_normal_actuator() {
    let (n, s) = sources(1e-5, 1.4e7, 1e7);
    let k = 1e7 * n.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let mut owner = make(&n, &s, k * DT / (4.0 * MU), 0.01);
    let (fa, fb) = (1.1 * k, -2.2 * k);
    let drive = [fa / det::sqrt(MA), 0.0, fb / det::sqrt(MB), 0.0];
    let before = owner.total_energy_j().unwrap();
    let p0 = momentum(&owner);
    let f = *owner.step(&drive, 0.2 * Z).unwrap();
    check(momentum(&owner), p0 + (fa + fb) * DT);
    check(
        center(&owner),
        0.01 * DT + (fa + fb) * DT * DT / (2.0 * (MA + MB)),
    );
    check(
        f.prescribed_normal_work_j,
        n.reversible_work_between(0.5 * Z, 0.2 * Z).unwrap(),
    );
    assert!(f.prescribed_normal_work_j != 0.0);
    check(
        owner.total_energy_j().unwrap() + f.friction_dissipation_j
            - f.external_work_j
            - f.prescribed_normal_work_j,
        before,
    );
    let native = owner.network().last_frame().unwrap();
    check(
        native.external_work_j,
        f.external_work_j + f.friction_work_on_network_j,
    );
    assert_eq!(native.total_dissipation_j(), 0.0);
    check(f.friction_work_on_network_j, -f.friction_dissipation_j);
    assert_eq!(owner.network().last_window_duration_s(), Some(DT));
    let ModuleWindowAudit::Available(audit) = owner.network().audit_last_window().unwrap() else {
        panic!("native port audit unavailable")
    };
    assert!(audit.violations().is_empty());
}

#[test]
fn native_rejects_unsupported_original_dynamics_and_mismatched_sources() {
    let (n, s) = sources(1e-5, 2e7, 1e7);
    let (_, alien) = sources(2e-5, 2e7, 1e7);
    let materials = n.normal_state().elastic().interface().surface_materials();
    for (elastic, connected, wrong_order, wrong_shear) in [
        (true, false, false, false),
        (false, true, false, false),
        (false, false, true, false),
        (false, false, false, true),
    ] {
        let ordered = if wrong_order {
            [&materials[1], &materials[0]]
        } else {
            [&materials[0], &materials[1]]
        };
        assert!(
            AdhesiveEventModalSystem::new(
                network(0.0, 0.0, 1.0, elastic, connected),
                &n,
                if wrong_shear { &alien } else { &s },
                ordered,
                0.5 * Z,
                0.0,
                column(),
                config(),
                &CancelGate::new_clock_free()
            )
            .is_err()
        );
    }
}

#[test]
fn native_cancel_source_range_and_late_body_budget_are_atomic_and_retryable() {
    let (n, s) = sources(1e-5, 2e7, 1e7);
    let materials = n.normal_state().elastic().interface().surface_materials();
    let create = || {
        AdhesiveEventModalSystem::new(
            network(0.0, 0.0, 5e-10, false, false),
            &n,
            &s,
            [&materials[0], &materials[1]],
            0.5 * Z,
            0.0,
            column(),
            config(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    };
    let (mut owner, mut reference) = (create(), create());
    owner.step(&[0.0; 4], 0.5 * Z).unwrap();
    reference.step(&[0.0; 4], 0.5 * Z).unwrap();
    let original = (
        states(&owner),
        owner.state(),
        *owner.last_frame().unwrap(),
        owner.network().last_frame().cloned(),
    );
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert!(matches!(
        owner.step_under_gate(&[0.0; 4], 0.5 * Z, &cancel),
        Err(AdhesiveEventModalError::Network(
            ModalCouplingError::Cancelled
        ))
    ));
    assert!(owner.step(&[0.0; 4], 100.0 * Z).is_err());
    // Only a noncontacted original coordinate violates its component budget.
    assert!(matches!(
        owner.step(&[0.0, 1.0, 0.0, 0.0], 0.5 * Z),
        Err(AdhesiveEventModalError::Network(
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
    assert_eq!(owner.network().last_window_duration_s(), Some(DT));
    let actual = *owner.step(&[0.0; 4], 0.5 * Z).unwrap();
    let expected = *reference.step(&[0.0; 4], 0.5 * Z).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(states(&owner), states(&reference));
}

#[test]
fn native_event_budget_refusal_keeps_initial_state_and_allows_a_valid_retry() {
    let (n, s) = sources(1e-5, 1.5e7, 1e7);
    let k = 1e7 * n.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let v0 = 3.0 * k * DT / ((1.0 + det::sqrt(3.0)) * MU);
    let materials = n.normal_state().elastic().interface().surface_materials();
    let mut bounds = config();
    bounds.events.maximum_events = 1;
    let mut owner = AdhesiveEventModalSystem::new(
        network(v0, 0.01, 1.0, false, false),
        &n,
        &s,
        [&materials[0], &materials[1]],
        0.5 * Z,
        0.0,
        column(),
        bounds,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let original = (states(&owner), owner.state());
    assert!(matches!(
        owner.step(&relative_force(-2.0 * k), 0.5 * Z),
        Err(AdhesiveEventModalError::Contact(
            fs_contact::normal_patch::dmt::events::DmtSlipEventError::EventBudget
        ))
    ));
    assert_eq!((states(&owner), owner.state()), original);
    assert_eq!(owner.samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
    owner.step(&[0.0; 4], 0.5 * Z).unwrap();
    assert_eq!(owner.samples_rendered(), 1);
}
