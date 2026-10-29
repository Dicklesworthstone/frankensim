//! G1/G3/G4: unlike wetting surfaces retain their histories under native body motion.
#[path = "../../fs-contact/tests/asymmetric_contact_line/common.rs"]
mod common;

use common::{GAP, VOLUME, close, source};
use fs_contact::capillary_bridge::asymmetric_contact_line::{
    ContactLineRegime, SourceAsymmetricContactLineBridge,
};
use fs_contact::capillary_bridge::contact_line::load_line::ContactLineLoadLineConfig;
use fs_contact::capillary_bridge::load_line::BridgeLoadLineConfig;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::capillary::{
    AsymmetricContactLineCapillaryModalConfig, AsymmetricContactLineCapillaryModalSystem,
    CapillaryModalError,
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
fn controls() -> AsymmetricContactLineCapillaryModalConfig {
    AsymmetricContactLineCapillaryModalConfig {
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
    source: &SourceAsymmetricContactLineBridge,
    angles: [f64; 2],
    connected: bool,
    fixed: bool,
) -> AsymmetricContactLineCapillaryModalSystem<'_> {
    AsymmetricContactLineCapillaryModalSystem::new_asymmetric_contact_line(
        network(connected),
        source,
        [
            &source.wetting()[0].surfaces()[0].material,
            &source.wetting()[1].surfaces()[0].material,
        ],
        source.initialize(GAP, angles).unwrap(),
        if fixed { vec![0.0; 2] } else { column() },
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn g1_unlike_surface_relaxation_and_reciprocal_native_impulses_share_actual_work() {
    let source = source([0.9, 1.1], [1.0, 3.0], [0.89, 1.09], [0.91, 1.11]);
    let mut system = owner(&source, [1.02, 1.23], false, false);
    let initial = system.total_energy_j().unwrap();
    let mut losses = [0.0; 2];
    for _ in 0..120 {
        let previous = system.state();
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
                if i == 0 { -1.0 } else { 1.0 } * frame.applied_mean_force_n / f64::from(RATE),
                2e-18,
            );
            assert!(frame.bridge.contact_line_dissipation_by_surface_j[i] >= 0.0);
            close(
                frame.bridge.contact_line_dissipation_by_surface_j[i],
                frame.bridge.viscous_line_dissipation_by_surface_j[i]
                    + frame.bridge.pinning_dissipation_by_surface_j[i],
                1e-26,
            );
            losses[i] += frame.bridge.contact_line_dissipation_by_surface_j[i];
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
        assert_eq!(
            frame.contact_line_dissipation_j,
            frame
                .bridge
                .contact_line_dissipation_by_surface_j
                .iter()
                .sum::<f64>(),
        );
        assert!(
            frame.compatibility_residual_m.abs() + frame.estimated_compatibility_error_m
                <= frame.compatibility_tolerance_m
        );
        assert!(
            frame.constitutive_residual_n.abs() + frame.estimated_constitutive_force_error_n
                <= frame.force_tolerance_n
        );
        assert!(
            frame.energy_residual_j.abs() + frame.estimated_constitutive_energy_error_j
                <= frame.energy_tolerance_j
        );
        assert!(frame.evaluations <= controls().contact.mechanical.max_evaluations);
        // Numerical probes and the actual-motion audit must represent ONE
        // physical relaxation from the original accepted checkpoint.
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
        close(
            system.total_energy_j().unwrap() + losses.iter().sum::<f64>(),
            initial,
            3e-16,
        );
    }
    assert!(losses.iter().all(|d| *d > 0.0));
    assert!(system.gap_m() < GAP);
    assert!(core::ptr::eq(system.source(), &source));
}

#[test]
fn g1_one_original_surface_stays_pinned_while_the_other_and_both_bodies_move() {
    let source = source([0.9, 1.1], [4.0, 15.0], [0.8, 1.09], [1.0, 1.11]);
    let mut system = owner(&source, [0.94, 1.2], false, false);
    let initial = system.state().wetted_spans_m();
    let energy = system.total_energy_j().unwrap();
    let mut loss = 0.0;
    for _ in 0..100 {
        let frame = *system.step(&[0.0; 2], 0.0).unwrap();
        assert_eq!(
            frame.bridge.regimes,
            [ContactLineRegime::Pinned, ContactLineRegime::Advancing]
        );
        assert_eq!(system.state().wetted_spans_m()[0], initial[0]);
        assert_eq!(frame.bridge.contact_line_dissipation_by_surface_j[0], 0.0);
        assert!(frame.bridge.viscous_line_dissipation_by_surface_j[1] > 0.0);
        assert!(frame.bridge.pinning_dissipation_by_surface_j[1] > 0.0);
        loss += frame.contact_line_dissipation_j;
        close(system.total_energy_j().unwrap() + loss, energy, 3e-16);
    }
    assert!(system.state().wetted_spans_m()[1] > initial[1]);
    assert!(system.gap_m() < GAP);
    assert!(system.components()[0].states()[0].displacement_m_sqrt_kg > 0.0);
    assert!(system.components()[1].states()[0].displacement_m_sqrt_kg < 0.0);
}

#[test]
fn g3_original_source_exchange_and_independent_friction_survive_native_coupling() {
    let slow_source = source([0.9, 1.1], [4.0, 15.0], [0.8, 1.09], [1.0, 1.11]);
    let fast_source = source([0.9, 1.1], [4.0, 3.0], [0.8, 1.09], [1.0, 1.11]);
    assert_eq!(
        slow_source.wetting()[0].resolved().identity(),
        fast_source.wetting()[0].resolved().identity()
    );
    assert_ne!(
        slow_source.wetting()[1].resolved().identity(),
        fast_source.wetting()[1].resolved().identity()
    );
    let mut reverse_wetting = fast_source.wetting().clone();
    reverse_wetting.swap(0, 1);
    let reversed_source = SourceAsymmetricContactLineBridge::new(
        fast_source.liquid().clone(),
        reverse_wetting,
        VOLUME,
        common::WIDTH,
        fast_source.limits(),
        fast_source.config(),
    )
    .unwrap();
    let mut slow = owner(&slow_source, [0.94, 1.2], true, false);
    let mut fast = owner(&fast_source, [0.94, 1.2], true, false);
    let mut reversed = owner(&reversed_source, [1.2, 0.94], true, false);
    for _ in 0..100 {
        slow.step(&[0.0; 2], 0.0).unwrap();
        let a = *fast.step(&[0.0; 2], 0.0).unwrap();
        let b = *reversed.step(&[0.0; 2], 0.0).unwrap();
        close(a.applied_mean_force_n, b.applied_mean_force_n, 2e-13);
        close(fast.gap_m(), reversed.gap_m(), 2e-16);
        close(a.surface_energy_j, b.surface_energy_j, 2e-18);
        for i in 0..2 {
            close(
                fast.state().wetted_spans_m()[i],
                reversed.state().wetted_spans_m()[1 - i],
                2e-15,
            );
            close(
                a.bridge.contact_line_dissipation_by_surface_j[i],
                b.bridge.contact_line_dissipation_by_surface_j[1 - i],
                2e-20,
            );
            close(
                fast.components()[i].states()[0].displacement_m_sqrt_kg,
                reversed.components()[i].states()[0].displacement_m_sqrt_kg,
                2e-18,
            );
        }
    }
    assert!(fast.state().wetted_spans_m()[1] > slow.state().wetted_spans_m()[1]);
    assert!(fast.state().contact_angles_radians()[1] < slow.state().contact_angles_radians()[1]);
    assert_ne!(
        fast.last_frame().unwrap().applied_mean_force_n,
        slow.last_frame().unwrap().applied_mean_force_n
    );
}

#[test]
fn g1_native_spring_damper_external_and_prescribed_work_retain_two_surface_losses() {
    let source = source([0.9, 1.1], [1.0, 3.0], [0.89, 1.09], [0.91, 1.11]);
    let mut system = owner(&source, [1.02, 1.23], true, false);
    let initial = system.total_energy_j().unwrap();
    let external = [0.0, 3e-6 / det::sqrt(MASS[1])];
    let (mut work, mut line_loss, mut network_loss) = (0.0, [0.0; 2], 0.0);
    for i in 0..100 {
        let frame = *system.step(&external, 5e-5).unwrap();
        work += frame.external_work_j + frame.prescribed_work_j;
        for (total, loss) in line_loss
            .iter_mut()
            .zip(frame.bridge.contact_line_dissipation_by_surface_j)
        {
            *total += loss;
        }
        network_loss += frame.network_dissipation_j;
        close(
            system.total_energy_j().unwrap() + line_loss.iter().sum::<f64>() + network_loss,
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
    assert!(line_loss.iter().all(|d| *d > 0.0) && network_loss > 0.0 && work != 0.0);
}

#[test]
fn g0_g4_original_state_source_budget_and_cancel_refusals_preserve_accepted_history() {
    let source = source([0.9, 1.1], [1.0, 3.0], [0.89, 1.09], [0.91, 1.11]);
    let different = common::source([0.9, 1.1], [1.0, 4.0], [0.89, 1.09], [0.91, 1.11]);
    let materials = [
        &source.wetting()[0].surfaces()[0].material,
        &source.wetting()[1].surfaces()[0].material,
    ];
    assert!(
        AsymmetricContactLineCapillaryModalSystem::new_asymmetric_contact_line(
            network(false),
            &source,
            materials,
            different.initialize(GAP, [1.02, 1.23]).unwrap(),
            column(),
            controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    assert!(
        AsymmetricContactLineCapillaryModalSystem::new_asymmetric_contact_line(
            network(false),
            &source,
            [materials[1], materials[0]],
            source.initialize(GAP, [1.02, 1.23]).unwrap(),
            column(),
            controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let mut system = owner(&source, [1.02, 1.23], false, false);
    let mut clean = owner(&source, [1.02, 1.23], false, false);
    assert_eq!(
        *system.step(&[0.0; 2], 0.0).unwrap(),
        *clean.step(&[0.0; 2], 0.0).unwrap()
    );
    let state = system.state();
    let frame = *system.last_frame().unwrap();
    let states = core::array::from_fn::<_, 2, _>(|i| system.components()[i].states()[0]);
    let gap = system.gap_m();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(system.step_under_gate(&[0.0; 2], 0.0, &cancelled).is_err());
    assert!(system.step(&[0.0; 2], 1e6).is_err());
    assert_eq!(system.state(), state);
    assert_eq!(system.gap_m(), gap);
    assert_eq!(system.samples_rendered(), 1);
    assert_eq!(*system.last_frame().unwrap(), frame);
    for (i, old) in states.iter().enumerate() {
        assert_eq!(system.components()[i].states()[0], *old);
    }
    assert_eq!(
        *system.step(&[0.0; 2], 0.0).unwrap(),
        *clean.step(&[0.0; 2], 0.0).unwrap()
    );

    let mut small = controls();
    small.contact.mechanical.max_evaluations = 2;
    let initial = source.initialize(GAP, [1.02, 1.23]).unwrap();
    let mut budgeted = AsymmetricContactLineCapillaryModalSystem::new_asymmetric_contact_line(
        network(false),
        &source,
        materials,
        initial,
        column(),
        small,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(budgeted.step(&[0.0; 2], 0.0).is_err());
    assert_eq!(budgeted.state(), initial);
    assert_eq!(budgeted.samples_rendered(), 0);
    assert!(budgeted.last_frame().is_none());
}

#[test]
fn g4_late_accuracy_refusal_discards_staged_body_motion_and_both_wetting_histories() {
    let source = source([0.9, 1.1], [1.0, 3.0], [0.89, 1.09], [0.91, 1.11]);
    let mut reference = owner(&source, [1.02, 1.23], false, false);
    let accepted = *reference.step(&[0.0; 2], 0.0).unwrap();
    let residual = accepted.energy_residual_j.abs();
    let estimate = accepted.estimated_constitutive_energy_error_j;
    assert!(estimate > residual && estimate > 0.0);
    assert!(reference.components()[0].states()[0].velocity_m_sqrt_kg_per_s > 0.0);
    // Real bodies and both wetting footprints have already been staged by the
    // time this final combined numerical-accuracy gate refuses publication.
    let absolute = residual + 0.5 * estimate;
    let initial = source.initialize(GAP, [1.02, 1.23]).unwrap();
    let mut strict = AsymmetricContactLineCapillaryModalSystem::new_asymmetric_contact_line(
        network_with_energy_budget(false, absolute, 1e-30),
        &source,
        [
            &source.wetting()[0].surfaces()[0].material,
            &source.wetting()[1].surfaces()[0].material,
        ],
        initial,
        column(),
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for _ in 0..2 {
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
}
