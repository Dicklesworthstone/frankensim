//! G1/G3/G4: actual circular wetting footprints and native moving bodies.
#[path = "../../fs-contact/tests/axisymmetric_contact_line/common.rs"]
mod common;
use common::{INITIAL_Q, VOLUME, close, source};
use fs_contact::capillary_bridge::axisymmetric::contact_line::{
    ContactLineRegime, SourceAxisymmetricContactLineBridge,
};
use fs_contact::capillary_bridge::contact_line::load_line::ContactLineLoadLineConfig;
use fs_contact::capillary_bridge::load_line::BridgeLoadLineConfig;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::capillary::{
    AxisymmetricContactLineCapillaryModalConfig, AxisymmetricContactLineCapillaryModalSystem,
    CapillaryModalError,
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
fn controls() -> AxisymmetricContactLineCapillaryModalConfig {
    AxisymmetricContactLineCapillaryModalConfig {
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
    source: &SourceAxisymmetricContactLineBridge,
    angle: f64,
    connected: bool,
) -> AxisymmetricContactLineCapillaryModalSystem<'_> {
    owner_with_budget(source, angle, connected, 4e-18, 2e-10)
}
fn owner_with_budget(
    source: &SourceAxisymmetricContactLineBridge,
    angle: f64,
    connected: bool,
    absolute: f64,
    relative: f64,
) -> AxisymmetricContactLineCapillaryModalSystem<'_> {
    let gate = CancelGate::new_clock_free();
    AxisymmetricContactLineCapillaryModalSystem::new_axisymmetric_contact_line(
        network(connected, absolute, relative),
        source,
        [
            &source.original().wetting()[0].surfaces()[0].material,
            &source.original().wetting()[1].surfaces()[0].material,
        ],
        source
            .initialize_at_suction(INITIAL_Q, angle, &gate)
            .unwrap(),
        column(),
        controls(),
        &gate,
    )
    .unwrap()
}
fn states(owner: &AxisymmetricContactLineCapillaryModalSystem<'_>) -> Vec<ModalAcousticState> {
    owner
        .components()
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}

#[test]
fn circular_line_relaxation_moves_original_bodies_with_reciprocal_work_and_momentum() {
    let source = source(1.1, 2.0, 1.08, 1.12);
    let mut system = owner(&source, 1.2, false);
    let initial_gap = system.gap_m();
    let initial_radius = system.state().wetted_radius_m();
    let initial_energy = system.total_energy_j().unwrap();
    let mut loss = 0.0;
    for _ in 0..24 {
        let previous = system.state();
        let old = states(&system);
        let frame = *system.step(&[0.0; 2], 0.0).unwrap();
        let next = states(&system);
        for i in 0..2 {
            close(
                det::sqrt(MASS[i])
                    * (next[i].velocity_m_sqrt_kg_per_s - old[i].velocity_m_sqrt_kg_per_s),
                if i == 0 { -1.0 } else { 1.0 } * frame.applied_mean_force_n / f64::from(RATE),
                3e-18,
            );
        }
        close(
            det::sqrt(MASS[0]) * next[0].velocity_m_sqrt_kg_per_s
                + det::sqrt(MASS[1]) * next[1].velocity_m_sqrt_kg_per_s,
            0.0,
            3e-18,
        );
        let actual_gap = initial_gap + next[1].displacement_m_sqrt_kg / det::sqrt(MASS[1])
            - next[0].displacement_m_sqrt_kg / det::sqrt(MASS[0]);
        close(system.gap_m(), actual_gap, 2e-18);
        assert_eq!(frame.bridge.profile.liquid_volume_m3, VOLUME);
        let radius = system.state().wetted_radius_m();
        for area in frame.bridge.profile.wetted_areas_m2 {
            let residual = frame.bridge.profile.geometric_radius_residual_m.abs();
            let tolerance = core::f64::consts::PI * residual * (2.0 * radius + residual)
                + 16.0 * f64::EPSILON * area.abs();
            close(area, core::f64::consts::PI * radius * radius, tolerance);
        }
        assert_eq!(frame.topology_release_j, 0.0);
        assert!(frame.contact_line_dissipation_j > 0.0);
        close(
            frame.contact_line_dissipation_j,
            frame.bridge.viscous_line_dissipation_j + frame.bridge.pinning_dissipation_j,
            1e-25,
        );
        assert!(
            frame.energy_residual_j.abs() + frame.estimated_constitutive_energy_error_j
                <= frame.energy_tolerance_j
        );
        // The published circular state represents one physical relaxation,
        // despite all intervening mechanical Newton and shape probes.
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
        close(
            system.total_energy_j().unwrap() + loss,
            initial_energy,
            1e-16,
        );
    }
    assert!(system.state().wetted_radius_m() > initial_radius);
    assert_ne!(system.gap_m(), initial_gap);
    assert_eq!(system.samples_rendered(), 24);
}

#[test]
fn a_pinned_circular_footprint_changes_its_angle_as_the_actual_gap_moves() {
    let source = source(1.1, 2.0, 1.0, 1.2);
    let mut system = owner(&source, 1.1, false);
    let radius = system.state().wetted_radius_m();
    let angle = system.state().contact_angle_radians();
    let gap = system.gap_m();
    let before = system.total_energy_j().unwrap();
    let frame = *system.step(&[0.0; 2], 5e-3).unwrap();
    assert_eq!(frame.bridge.regime, ContactLineRegime::Pinned);
    assert_eq!(system.state().wetted_radius_m(), radius);
    assert_ne!(system.state().contact_angle_radians(), angle);
    assert_ne!(system.gap_m(), gap);
    assert_eq!(frame.contact_line_dissipation_j, 0.0);
    close(
        system.total_energy_j().unwrap() - before,
        frame.prescribed_work_j,
        1e-17,
    );
    assert_ne!(states(&system)[0].velocity_m_sqrt_kg_per_s, 0.0);
}

#[test]
fn original_line_friction_changes_circular_relaxation_and_native_traction() {
    let slow = source(1.1, 8.0, 1.08, 1.12);
    let fast = source(1.1, 2.0, 1.08, 1.12);
    let mut a = owner(&slow, 1.2, false);
    let mut b = owner(&fast, 1.2, false);
    let radius_a = a.state().wetted_radius_m();
    let radius_b = b.state().wetted_radius_m();
    let fa = *a.step(&[0.0; 2], 0.0).unwrap();
    let fb = *b.step(&[0.0; 2], 0.0).unwrap();
    assert!(b.state().wetted_radius_m() - radius_b > a.state().wetted_radius_m() - radius_a);
    assert_ne!(fa.applied_mean_force_n, fb.applied_mean_force_n);
    assert_ne!(states(&a), states(&b));
}

#[test]
fn original_connections_external_forces_and_prescribed_gap_share_the_3d_surface_account() {
    let source = source(1.1, 2.0, 1.08, 1.12);
    let mut system = owner(&source, 1.18, true);
    let initial = system.total_energy_j().unwrap();
    let external = [1e-5 / det::sqrt(MASS[0]), -2e-5 / det::sqrt(MASS[1])];
    let (mut supplied, mut loss) = (0.0, 0.0);
    let mut any_native_loss = false;
    for step in 0..16 {
        let frame = *system
            .step(&external, if step < 8 { 1e-4 } else { -5e-5 })
            .unwrap();
        supplied += frame.external_work_j + frame.prescribed_work_j;
        loss += frame.network_dissipation_j + frame.contact_line_dissipation_j;
        any_native_loss |= frame.network_dissipation_j > 0.0;
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
    assert!(any_native_loss);
    assert_ne!(supplied, 0.0);
}

#[test]
fn source_assignment_and_failed_steps_retain_original_circular_and_native_history() {
    let source = source(1.1, 2.0, 1.08, 1.12);
    let mut system = owner(&source, 1.2, true);
    let mut reference = owner(&source, 1.2, true);
    system.step(&[0.0; 2], 0.0).unwrap();
    reference.step(&[0.0; 2], 0.0).unwrap();
    let old = (
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
        old
    );
    assert_eq!(system.samples_rendered(), 1);
    assert_eq!(
        *system.step(&[0.0; 2], 0.0).unwrap(),
        *reference.step(&[0.0; 2], 0.0).unwrap()
    );
    assert_eq!(states(&system), states(&reference));
    let gate = CancelGate::new_clock_free();
    assert!(
        AxisymmetricContactLineCapillaryModalSystem::new_axisymmetric_contact_line(
            network(false, 4e-18, 2e-10),
            &source,
            [
                &source.original().wetting()[1].surfaces()[0].material,
                &source.original().wetting()[0].surfaces()[0].material
            ],
            source.initialize_at_suction(INITIAL_Q, 1.2, &gate).unwrap(),
            column(),
            controls(),
            &gate,
        )
        .is_err()
    );
}

#[test]
fn unresolved_shape_work_cannot_publish_native_motion_or_circular_footprints() {
    let source = source(1.1, 2.0, 1.08, 1.12);
    let mut ordinary = owner(&source, 1.2, false);
    let frame = *ordinary.step(&[0.0; 2], 0.0).unwrap();
    assert!(frame.estimated_constitutive_energy_error_j > 0.0);
    let mut strict = owner_with_budget(
        &source,
        1.2,
        false,
        0.5 * frame.estimated_constitutive_energy_error_j,
        1e-30,
    );
    let old = (states(&strict), strict.state(), strict.gap_m());
    assert!(matches!(
        strict.step(&[0.0; 2], 0.0),
        Err(CapillaryModalError::EnergyAccuracy { .. })
    ));
    assert_eq!((states(&strict), strict.state(), strict.gap_m()), old);
    assert_eq!(strict.samples_rendered(), 0);
    assert!(strict.last_frame().is_none());
}
