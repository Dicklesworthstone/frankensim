//! G1/G3/G4: unlike original circular wetting surfaces move actual native bodies.
#[path = "../../fs-contact/tests/asymmetric_circular_contact_line/common.rs"]
mod common;
use common::{INITIAL_Q, VOLUME, close, source};
use fs_contact::capillary_bridge::axisymmetric::asymmetric_contact_line::SourceAsymmetricCircularContactLineBridge;
use fs_contact::capillary_bridge::contact_line::ContactLineRegime;
use fs_contact::capillary_bridge::contact_line::load_line::ContactLineLoadLineConfig;
use fs_contact::capillary_bridge::load_line::BridgeLoadLineConfig;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::capillary::{
    AsymmetricCircularContactLineCapillaryModalConfig,
    AsymmetricCircularContactLineCapillaryModalSystem, CapillaryModalError,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

const MASS: [f64; 2] = [2e-4, 3e-4];
const RATE: u32 = 50_000;
fn network(connected: bool, absolute: f64, relative: f64) -> CoupledModalSystem {
    let models = MASS
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
    let links = if connected {
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
        models,
        links,
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
fn controls() -> AsymmetricCircularContactLineCapillaryModalConfig {
    AsymmetricCircularContactLineCapillaryModalConfig {
        contact: ContactLineLoadLineConfig {
            mechanical: BridgeLoadLineConfig {
                max_iterations: 32,
                max_evaluations: 161,
                max_backtracks: 16,
                gap_absolute_tolerance_m: 3e-15,
                gap_relative_tolerance: 1e-10,
                maximum_gap_increment_m: 2e-5,
                minimum_compatibility_slope: 0.05,
            },
            force_absolute_tolerance_n: 2e-10,
            force_relative_tolerance: 2e-7,
        },
        force_absolute_tolerance_n: 6e-10,
        force_relative_tolerance: 6e-7,
    }
}
fn column() -> Vec<f64> {
    vec![-1.0 / det::sqrt(MASS[0]), 1.0 / det::sqrt(MASS[1])]
}
fn owner(
    source: &SourceAsymmetricCircularContactLineBridge,
    angles: [f64; 2],
    connected: bool,
) -> AsymmetricCircularContactLineCapillaryModalSystem<'_> {
    owner_with_budget(source, angles, connected, 4e-18, 2e-10)
}
fn owner_with_budget(
    source: &SourceAsymmetricCircularContactLineBridge,
    angles: [f64; 2],
    connected: bool,
    absolute: f64,
    relative: f64,
) -> AsymmetricCircularContactLineCapillaryModalSystem<'_> {
    let gate = CancelGate::new_clock_free();
    AsymmetricCircularContactLineCapillaryModalSystem::new_asymmetric_circular_contact_line(
        network(connected, absolute, relative),
        source,
        [
            &source.original().wetting()[0].surfaces()[0].material,
            &source.original().wetting()[1].surfaces()[0].material,
        ],
        source
            .initialize_at_suction(INITIAL_Q, angles, &gate)
            .unwrap(),
        column(),
        controls(),
        &gate,
    )
    .unwrap()
}
fn states(
    owner: &AsymmetricCircularContactLineCapillaryModalSystem<'_>,
) -> Vec<ModalAcousticState> {
    owner
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}

#[test]
fn one_original_footprint_stays_pinned_while_the_other_advances_and_moves_both_bodies() {
    let source = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let mut system = owner(&source, [1.05, 1.28], false);
    let gap0 = system.gap_m();
    let radii0 = system.state().wetted_radii_m();
    let energy0 = system.total_energy_j().unwrap();
    let mut loss = 0.0;
    for _ in 0..12 {
        let previous = system.state();
        let old = states(&system);
        let frame = *system.step(&[0.0; 2], 0.0).unwrap();
        let next = states(&system);
        assert_eq!(
            frame.bridge.regimes,
            [ContactLineRegime::Pinned, ContactLineRegime::Advancing]
        );
        assert_eq!(system.state().wetted_radii_m()[0], radii0[0]);
        assert!(system.state().wetted_radii_m()[1] > previous.wetted_radii_m()[1]);
        assert_eq!(frame.bridge.per_line_viscous_dissipation_j[0], 0.0);
        assert_eq!(frame.bridge.per_line_pinning_dissipation_j[0], 0.0);
        assert!(frame.bridge.per_line_viscous_dissipation_j[1] > 0.0);
        close(
            frame.contact_line_dissipation_j,
            frame
                .bridge
                .per_line_viscous_dissipation_j
                .iter()
                .sum::<f64>()
                + frame
                    .bridge
                    .per_line_pinning_dissipation_j
                    .iter()
                    .sum::<f64>(),
            1e-25,
        );
        for i in 0..2 {
            close(
                det::sqrt(MASS[i])
                    * (next[i].velocity_m_sqrt_kg_per_s - old[i].velocity_m_sqrt_kg_per_s),
                if i == 0 { -1.0 } else { 1.0 } * frame.applied_mean_force_n / f64::from(RATE),
                3e-18,
            );
            let radius = system.state().wetted_radii_m()[i];
            let residual = frame.bridge.profile.geometric_radius_residuals_m[i].abs();
            let area = frame.bridge.profile.wetted_areas_m2[i];
            close(
                area,
                core::f64::consts::PI * radius * radius,
                core::f64::consts::PI * residual * (2.0 * radius + residual)
                    + 16.0 * f64::EPSILON * area.abs(),
            );
        }
        close(
            det::sqrt(MASS[0]) * next[0].velocity_m_sqrt_kg_per_s
                + det::sqrt(MASS[1]) * next[1].velocity_m_sqrt_kg_per_s,
            0.0,
            3e-18,
        );
        close(
            system.gap_m(),
            gap0 + next[1].displacement_m_sqrt_kg / det::sqrt(MASS[1])
                - next[0].displacement_m_sqrt_kg / det::sqrt(MASS[0]),
            2e-18,
        );
        assert_eq!(frame.bridge.profile.liquid_volume_m3, VOLUME);
        assert_eq!(frame.topology_release_j, 0.0);
        assert!(
            frame.energy_residual_j.abs() + frame.estimated_constitutive_energy_error_j
                <= frame.energy_tolerance_j
        );
        assert_eq!(
            frame.bridge,
            source
                .advance(
                    previous,
                    system.gap_m(),
                    system.sample_period_s(),
                    &CancelGate::new_clock_free()
                )
                .unwrap()
        );
        loss += frame.contact_line_dissipation_j;
        close(system.total_energy_j().unwrap() + loss, energy0, 1e-16);
    }
    assert_ne!(system.gap_m(), gap0);
    assert_ne!(system.state().contact_angles_radians()[0], 1.05);
}

#[test]
fn changing_one_original_line_resistance_changes_its_relaxation_and_native_traction() {
    let slow = source([1.05, 1.15], [8.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let fast = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let mut a = owner(&slow, [1.19, 1.28], false);
    let mut b = owner(&fast, [1.19, 1.28], false);
    let r0 = a.state().wetted_radii_m();
    let fa = *a.step(&[0.0; 2], 0.0).unwrap();
    let fb = *b.step(&[0.0; 2], 0.0).unwrap();
    assert!(b.state().wetted_radii_m()[0] - r0[0] > a.state().wetted_radii_m()[0] - r0[0]);
    assert_ne!(fa.applied_mean_force_n, fb.applied_mean_force_n);
    assert_ne!(states(&a), states(&b));
    assert_eq!(fa.bridge.regimes, [ContactLineRegime::Advancing; 2]);
    assert_eq!(fb.bridge.regimes, [ContactLineRegime::Advancing; 2]);
}

#[test]
fn exchanging_unlike_surfaces_swaps_footprints_and_preserves_native_normal_motion() {
    let lower = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let upper = source([1.15, 1.05], [4.0, 2.0], [1.1, 1.0], [1.2, 1.1]);
    let mut a = owner(&lower, [1.05, 1.28], false);
    let mut b = owner(&upper, [1.28, 1.05], false);
    for _ in 0..5 {
        let fa = *a.step(&[0.0; 2], 0.0).unwrap();
        let fb = *b.step(&[0.0; 2], 0.0).unwrap();
        close(fa.applied_mean_force_n, fb.applied_mean_force_n, 5e-13);
        close(a.gap_m(), b.gap_m(), 2e-15);
        for i in 0..2 {
            close(
                a.state().wetted_radii_m()[i],
                b.state().wetted_radii_m()[1 - i],
                2e-15,
            );
            close(
                fa.bridge.per_line_viscous_dissipation_j[i],
                fb.bridge.per_line_viscous_dissipation_j[1 - i],
                2e-20,
            );
        }
        for (sa, sb) in states(&a).iter().zip(states(&b)) {
            close(sa.displacement_m_sqrt_kg, sb.displacement_m_sqrt_kg, 1e-17);
            close(
                sa.velocity_m_sqrt_kg_per_s,
                sb.velocity_m_sqrt_kg_per_s,
                2e-14,
            );
        }
    }
}

#[test]
fn both_contact_lines_share_original_connections_external_work_and_prescribed_gap_work() {
    let source = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let mut system = owner(&source, [1.19, 1.28], true);
    let initial = system.total_energy_j().unwrap();
    let external = [1e-5 / det::sqrt(MASS[0]), -2e-5 / det::sqrt(MASS[1])];
    let (mut supplied, mut loss) = (0.0, 0.0);
    let mut native_loss = 0.0;
    for i in 0..12 {
        let frame = *system
            .step(&external, if i < 6 { 1e-4 } else { -5e-5 })
            .unwrap();
        supplied += frame.external_work_j + frame.prescribed_work_j;
        loss += frame.network_dissipation_j + frame.contact_line_dissipation_j;
        native_loss += frame.network_dissipation_j;
        close(
            system.total_energy_j().unwrap() + loss - supplied,
            initial,
            2e-16,
        );
        assert!(
            frame.compatibility_residual_m.abs() + frame.estimated_compatibility_error_m
                <= frame.compatibility_tolerance_m
        );
    }
    assert!(native_loss > 0.0);
    assert_ne!(supplied, 0.0);
}

#[test]
fn source_order_cancellation_and_refused_motion_preserve_both_original_footprint_histories() {
    let source = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let (mut system, mut reference) = (
        owner(&source, [1.05, 1.28], true),
        owner(&source, [1.05, 1.28], true),
    );
    system.step(&[0.0; 2], 0.0).unwrap();
    reference.step(&[0.0; 2], 0.0).unwrap();
    let before = (
        states(&system),
        system.state(),
        system.gap_m(),
        *system.last_frame().unwrap(),
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(system.step_under_gate(&[0.0; 2], 0.0, &gate).is_err());
    assert!(system.step(&[0.0; 2], 1000.0).is_err());
    assert_eq!(
        (
            states(&system),
            system.state(),
            system.gap_m(),
            *system.last_frame().unwrap()
        ),
        before
    );
    assert_eq!(system.samples_rendered(), 1);
    assert_eq!(
        *system.step(&[0.0; 2], 0.0).unwrap(),
        *reference.step(&[0.0; 2], 0.0).unwrap()
    );
    assert_eq!(states(&system), states(&reference));
    let gate = CancelGate::new_clock_free();
    assert!(
        AsymmetricCircularContactLineCapillaryModalSystem::new_asymmetric_circular_contact_line(
            network(false, 4e-18, 2e-10),
            &source,
            [
                &source.original().wetting()[1].surfaces()[0].material,
                &source.original().wetting()[0].surfaces()[0].material
            ],
            source
                .initialize_at_suction(INITIAL_Q, [1.05, 1.28], &gate)
                .unwrap(),
            column(),
            controls(),
            &gate
        )
        .is_err()
    );
}

#[test]
fn unresolved_coupled_shape_work_cannot_publish_any_native_or_wetting_state() {
    let source = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let mut ordinary = owner(&source, [1.05, 1.28], false);
    let frame = *ordinary.step(&[0.0; 2], 0.0).unwrap();
    assert!(frame.estimated_constitutive_energy_error_j > 0.0);
    let mut strict = owner_with_budget(
        &source,
        [1.05, 1.28],
        false,
        0.5 * frame.estimated_constitutive_energy_error_j,
        1e-30,
    );
    let before = (states(&strict), strict.state(), strict.gap_m());
    assert!(matches!(
        strict.step(&[0.0; 2], 0.0),
        Err(CapillaryModalError::EnergyAccuracy { .. })
    ));
    assert_eq!((states(&strict), strict.state(), strict.gap_m()), before);
    assert_eq!(strict.samples_rendered(), 0);
    assert!(strict.last_frame().is_none());
}
