//! G1/G3/G4: finite bodies coupled to the original dissipative wetting history.
#[path = "../../fs-contact/tests/contact_line_bridge/common.rs"]
mod common;

use common::{GAP, VOLUME, close, source};
use fs_contact::capillary_bridge::contact_line::load_line::ContactLineLoadLineConfig;
use fs_contact::capillary_bridge::contact_line::{ContactLineRegime, SourceContactLineBridge};
use fs_contact::capillary_bridge::load_line::BridgeLoadLineConfig;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::capillary::{
    CapillaryModalError, ContactLineCapillaryModalConfig, ContactLineCapillaryModalSystem,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

const MASS: [f64; 2] = [2e-4, 3e-4];
const RATE: u32 = 50_000;

fn network(connected: bool) -> CoupledModalSystem {
    network_with_energy_budget(connected, 2e-18, 2e-10)
}
fn network_with_energy_budget(connected: bool, absolute: f64, relative: f64) -> CoupledModalSystem {
    let components = MASS
        .into_iter()
        .map(|_| {
            ModalAcousticTimeModel::try_new_with_free_coordinates(
                RATE,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap()
        })
        .collect();
    let connections = if connected {
        vec![ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0 / det::sqrt(MASS[0])],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![1.0 / det::sqrt(MASS[1])],
            },
            stiffness_n_m: 0.1,
            damping_n_s_m: 1e-4,
            rest_extension_m: 0.0,
        }]
    } else {
        vec![]
    };
    CoupledModalSystem::new(
        components,
        connections,
        ModalCouplingConfig {
            max_modes: 2,
            max_connections: 1,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 0.01,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: absolute,
            energy_relative_tolerance: relative,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn controls() -> ContactLineCapillaryModalConfig {
    ContactLineCapillaryModalConfig {
        contact: ContactLineLoadLineConfig {
            mechanical: BridgeLoadLineConfig {
                max_iterations: 32,
                max_evaluations: 161,
                max_backtracks: 16,
                gap_absolute_tolerance_m: 2e-15,
                gap_relative_tolerance: 1e-10,
                maximum_gap_increment_m: 2e-5,
                minimum_compatibility_slope: 0.05,
            },
            force_absolute_tolerance_n: 2e-12,
            force_relative_tolerance: 2e-9,
        },
        force_absolute_tolerance_n: 6e-12,
        force_relative_tolerance: 6e-9,
    }
}
fn column() -> Vec<f64> {
    vec![-1.0 / det::sqrt(MASS[0]), 1.0 / det::sqrt(MASS[1])]
}
fn owner(
    source: &SourceContactLineBridge,
    angle: f64,
    connected: bool,
    fixed: bool,
) -> ContactLineCapillaryModalSystem<'_> {
    ContactLineCapillaryModalSystem::new_contact_line(
        network(connected),
        source,
        [
            &source.original().wetting()[0].surfaces()[0].material,
            &source.original().wetting()[1].surfaces()[0].material,
        ],
        source.initialize(GAP, angle).unwrap(),
        if fixed { vec![0.0; 2] } else { column() },
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn g1_native_surface_relaxation_body_motion_and_line_loss_share_actual_work() {
    let source = source(1.0, 1.0, 0.99, 1.01);
    let mut system = owner(&source, 1.08, false, false);
    let initial = system.total_energy_j().unwrap();
    let mut line_loss = 0.0;
    for _ in 0..120 {
        let v0 = core::array::from_fn::<_, 2, _>(|i| {
            system.components()[i].states()[0].velocity_m_sqrt_kg_per_s
        });
        let frame = *system.step(&[0.0; 2], 0.0).unwrap();
        let v1 = core::array::from_fn::<_, 2, _>(|i| {
            system.components()[i].states()[0].velocity_m_sqrt_kg_per_s
        });
        for i in 0..2 {
            close(
                det::sqrt(MASS[i]) * (v1[i] - v0[i]),
                if i == 0 {
                    -frame.applied_mean_force_n / f64::from(RATE)
                } else {
                    frame.applied_mean_force_n / f64::from(RATE)
                },
                2e-18,
            );
        }
        close(
            det::sqrt(MASS[0]) * v1[0] + det::sqrt(MASS[1]) * v1[1],
            0.0,
            2e-18,
        );
        let physical_gap = GAP
            + system.components()[1].states()[0].displacement_m_sqrt_kg / det::sqrt(MASS[1])
            - system.components()[0].states()[0].displacement_m_sqrt_kg / det::sqrt(MASS[0]);
        close(system.gap_m(), physical_gap, 2e-18);
        assert_eq!(frame.bridge.profile.liquid_volume_m3, VOLUME);
        assert_eq!(frame.topology_release_j, 0.0);
        assert!(frame.contact_line_dissipation_j >= 0.0);
        assert_eq!(
            frame.contact_line_dissipation_j,
            frame.bridge.contact_line_dissipation_j
        );
        assert!(
            frame.compatibility_residual_m.abs() + frame.estimated_compatibility_error_m
                <= frame.compatibility_tolerance_m
        );
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(
            frame.energy_residual_j.abs() + frame.estimated_constitutive_energy_error_j
                <= frame.energy_tolerance_j
        );
        assert!(frame.evaluations <= controls().contact.mechanical.max_evaluations);
        line_loss += frame.contact_line_dissipation_j;
        close(system.total_energy_j().unwrap() + line_loss, initial, 3e-16);
    }
    assert!(line_loss > 0.0);
    assert!(system.gap_m() < GAP);
    assert!(core::ptr::eq(system.source(), &source));
}

#[test]
fn g1_original_line_friction_changes_relaxation_and_pinning_retains_footprint() {
    let fast_source = source(1.0, 0.5, 0.99, 1.01);
    let slow_source = source(1.0, 5.0, 0.99, 1.01);
    let mut fast = owner(&fast_source, 1.1, false, true);
    let mut slow = owner(&slow_source, 1.1, false, true);
    let initial = fast.total_energy_j().unwrap();
    close(initial, slow.total_energy_j().unwrap(), 1e-22);
    let mut fast_loss = 0.0;
    for _ in 0..100 {
        let a = *fast.step(&[0.0; 2], 0.0).unwrap();
        let b = *slow.step(&[0.0; 2], 0.0).unwrap();
        fast_loss += a.contact_line_dissipation_j;
        assert_eq!(
            a.external_work_j + a.prescribed_work_j + a.network_energy_j,
            0.0
        );
        assert_eq!(
            b.external_work_j + b.prescribed_work_j + b.network_energy_j,
            0.0
        );
        close(fast.total_energy_j().unwrap() + fast_loss, initial, 2e-16);
    }
    assert!(fast.state().contact_angle_radians() < slow.state().contact_angle_radians() - 0.005);
    assert!(fast.state().wetted_span_m() > slow.state().wetted_span_m());

    let pinned_source = source(1.0, 1.0, 0.9, 1.1);
    let mut pinned = owner(&pinned_source, 1.0, false, true);
    let footprint = pinned.state().wetted_span_m();
    let frame = *pinned.step(&[0.0; 2], 1e-4).unwrap();
    assert_eq!(frame.bridge.regime, ContactLineRegime::Pinned);
    assert_eq!(pinned.state().wetted_span_m(), footprint);
    assert_eq!(frame.contact_line_dissipation_j, 0.0);
    assert!(frame.prescribed_work_j > 0.0);
}

#[test]
fn g1_native_connection_external_and_prescribed_work_keep_separate_losses() {
    let source = source(1.0, 1.0, 0.99, 1.01);
    let mut system = owner(&source, 1.08, true, false);
    let initial = system.total_energy_j().unwrap();
    let external = [0.0, 3e-6 / det::sqrt(MASS[1])];
    let (mut work, mut line_loss, mut network_loss) = (0.0, 0.0, 0.0);
    for i in 0..100 {
        let frame = *system.step(&external, 5e-5).unwrap();
        work += frame.external_work_j + frame.prescribed_work_j;
        line_loss += frame.contact_line_dissipation_j;
        network_loss += frame.network_dissipation_j;
        close(
            system.total_energy_j().unwrap() + line_loss + network_loss,
            initial + work,
            3e-16,
        );
        let momentum = (0..2)
            .map(|j| {
                det::sqrt(MASS[j]) * system.components()[j].states()[0].velocity_m_sqrt_kg_per_s
            })
            .sum::<f64>();
        close(momentum, f64::from(i + 1) * 3e-6 / f64::from(RATE), 2e-18);
    }
    assert!(line_loss > 0.0 && network_loss > 0.0 && work != 0.0);
}

#[test]
fn g0_g4_original_state_budget_and_cancel_refusals_are_atomic_and_retryable() {
    let source = source(1.0, 1.0, 0.99, 1.01);
    let different = common::source(1.0, 2.0, 0.99, 1.01);
    let materials = [
        &source.original().wetting()[0].surfaces()[0].material,
        &source.original().wetting()[1].surfaces()[0].material,
    ];
    assert!(
        ContactLineCapillaryModalSystem::new_contact_line(
            network(false),
            &source,
            materials,
            different.initialize(GAP, 1.08).unwrap(),
            column(),
            controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    assert!(
        ContactLineCapillaryModalSystem::new_contact_line(
            network(false),
            &source,
            [materials[1], materials[0]],
            source.initialize(GAP, 1.08).unwrap(),
            column(),
            controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let mut system = owner(&source, 1.08, false, false);
    let mut clean = owner(&source, 1.08, false, false);
    let expected = *clean.step(&[0.0; 2], 0.0).unwrap();
    let state = system.state();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(system.step_under_gate(&[0.0; 2], 0.0, &cancelled).is_err());
    assert!(system.step(&[0.0; 2], 1e6).is_err());
    assert_eq!(system.state(), state);
    assert_eq!(system.gap_m(), GAP);
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.last_frame().is_none());
    for component in system.components() {
        assert_eq!(component.states()[0].displacement_m_sqrt_kg, 0.0);
        assert_eq!(component.states()[0].velocity_m_sqrt_kg_per_s, 0.0);
    }
    assert_eq!(*system.step(&[0.0; 2], 0.0).unwrap(), expected);

    let mut small = controls();
    small.contact.mechanical.max_evaluations = 2;
    let mut budgeted = ContactLineCapillaryModalSystem::new_contact_line(
        network(false),
        &source,
        materials,
        source.initialize(GAP, 1.08).unwrap(),
        column(),
        small,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(budgeted.step(&[0.0; 2], 0.0).is_err());
    assert_eq!(budgeted.samples_rendered(), 0);
    assert!(budgeted.last_frame().is_none());
}

#[test]
fn g4_native_energy_budget_includes_original_liquid_numerical_uncertainty() {
    let source = source(1.0, 1.0, 0.99, 1.01);
    let mut reference = owner(&source, 1.08, false, true);
    let accepted = *reference.step(&[0.0; 2], 0.0).unwrap();
    let residual = accepted.energy_residual_j.abs();
    let estimate = accepted.estimated_constitutive_energy_error_j;
    assert!(estimate > residual && estimate > 0.0);
    // This budget admits the measured arithmetic residual but cannot also
    // admit the original geometry/quadrature uncertainty. The combined owner
    // must refuse before publishing either its clock or the liquid footprint.
    let absolute = residual + 0.5 * estimate;
    let initial = source.initialize(GAP, 1.08).unwrap();
    let mut strict = ContactLineCapillaryModalSystem::new_contact_line(
        network_with_energy_budget(false, absolute, 1e-30),
        &source,
        [
            &source.original().wetting()[0].surfaces()[0].material,
            &source.original().wetting()[1].surfaces()[0].material,
        ],
        initial,
        vec![0.0; 2],
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(matches!(
        strict.step(&[0.0; 2], 0.0),
        Err(CapillaryModalError::EnergyAccuracy { .. })
    ));
    assert_eq!(strict.state(), initial);
    assert_eq!(strict.gap_m(), GAP);
    assert_eq!(strict.samples_rendered(), 0);
    assert!(strict.last_frame().is_none());
    for component in strict.components() {
        assert_eq!(component.states()[0].displacement_m_sqrt_kg, 0.0);
        assert_eq!(component.states()[0].velocity_m_sqrt_kg_per_s, 0.0);
    }
}
