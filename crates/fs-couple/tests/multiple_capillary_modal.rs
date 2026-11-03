//! G1/G3/G4: simultaneous original bridges, complete bodies and atomic history.
#[path = "../../fs-contact/tests/asymmetric_circular_contact_line/common.rs"]
mod common;
use common::{INITIAL_Q, close, source};
use fs_contact::capillary_bridge::axisymmetric::asymmetric_contact_line::SourceAsymmetricCircularContactLineBridge;
use fs_contact::capillary_bridge::contact_line::ContactLineRegime;
use fs_contact::capillary_bridge::contact_line::load_line::ContactLineLoadLineConfig;
use fs_contact::capillary_bridge::load_line::BridgeLoadLineConfig;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::capillary::{
    AsymmetricCircularContactLineCapillaryModalConfig,
    AsymmetricCircularContactLineCapillaryModalSystem, CapillaryModalError, multiple::*,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, PreparedModalBudget,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

const MASS: [f64; 3] = [2e-4, 3e-4, 5e-4];
const RATE: u32 = 50_000;

fn controls() -> MultipleCircularWettingModalConfig {
    MultipleCircularWettingModalConfig {
        maximum_bridges: 4,
        maximum_setup_terms: 4096,
        maximum_iterations: 24,
        maximum_backtracks: 16,
        maximum_evaluations: 256,
        gap_absolute_tolerance_m: 3e-15,
        gap_relative_tolerance: 1e-10,
        maximum_gap_increment_m: 2e-5,
        force_absolute_tolerance_n: 6e-10,
        force_relative_tolerance: 6e-7,
        maximum_compatibility_gain: 10.0,
    }
}
fn single_controls() -> AsymmetricCircularContactLineCapillaryModalConfig {
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
fn network(
    elastic: bool,
    absolute: f64,
    relative: f64,
    displacement_limit: f64,
) -> CoupledModalSystem {
    let models = MASS
        .into_iter()
        .enumerate()
        .map(|(i, mass)| {
            let mut budget = ModalAcousticTimeBudget::audible_reference();
            budget.maximum_abs_displacement_m_sqrt_kg = displacement_limit;
            let mut modes = vec![ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            }];
            if elastic {
                modes.push(ModalAcousticMode {
                    angular_frequency_rad_s: [180.0, 230.0, 270.0][i],
                    damping_ratio: 0.03,
                    pressure_per_modal_velocity: C64 { re: 0.5, im: 0.1 },
                });
            }
            let mut model =
                ModalAcousticTimeModel::try_new_with_free_coordinates(RATE, modes, budget).unwrap();
            if elastic {
                model
                    .restore_states(&[
                        ModalAcousticState {
                            displacement_m_sqrt_kg: det::sqrt(mass) * [2e-8, -1e-8, 1e-8][i],
                            velocity_m_sqrt_kg_per_s: 0.0,
                        },
                        ModalAcousticState {
                            displacement_m_sqrt_kg: det::sqrt(mass) * [3e-8, -2e-8, 1e-8][i],
                            velocity_m_sqrt_kg_per_s: det::sqrt(mass)
                                * [0.0003, -0.0005, 0.0002][i],
                        },
                    ])
                    .unwrap();
            }
            model
        })
        .collect();
    let links = if elastic {
        (0..2)
            .map(|i| ModalConnection {
                left: ModalAttachment {
                    component: i,
                    shapes: vec![1.0 / det::sqrt(MASS[i]), 0.2 / det::sqrt(MASS[i])],
                },
                right: ModalAttachment {
                    component: i + 1,
                    shapes: vec![1.0 / det::sqrt(MASS[i + 1]), -0.1 / det::sqrt(MASS[i + 1])],
                },
                stiffness_n_m: [0.04, 0.025][i],
                damping_n_s_m: [2e-7, 1e-7][i],
                rest_extension_m: [1e-8, -0.5e-8][i],
            })
            .collect()
    } else {
        vec![]
    };
    CoupledModalSystem::new_prepared(
        models,
        links,
        ModalCouplingConfig {
            max_modes: 6,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 0.01,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: absolute,
            energy_relative_tolerance: relative,
        },
        PreparedModalBudget {
            max_modes: 6,
            max_dense_bytes: 1024 * 1024,
            max_multiply_terms: 100_000_000,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn columns(elastic: bool) -> Vec<Vec<f64>> {
    let width = if elastic { 2 } else { 1 };
    (0..2)
        .map(|i| {
            let mut b = vec![0.0; 3 * width];
            b[i * width] = -1.0 / det::sqrt(MASS[i]);
            b[(i + 1) * width] = 1.0 / det::sqrt(MASS[i + 1]);
            if elastic {
                b[i * width + 1] = [-0.25, 0.3][i] / det::sqrt(MASS[i]);
                b[(i + 1) * width + 1] = [0.4, -0.2][i] / det::sqrt(MASS[i + 1]);
            }
            b
        })
        .collect()
}
fn bridges<'a>(
    sources: [&'a SourceAsymmetricCircularContactLineBridge; 2],
    elastic: bool,
) -> Vec<CircularWettingBridge<'a>> {
    let gate = CancelGate::new_clock_free();
    columns(elastic)
        .into_iter()
        .enumerate()
        .map(|(i, column)| CircularWettingBridge {
            source: sources[i],
            substrate_materials: [
                &sources[i].original().wetting()[0].surfaces()[0].material,
                &sources[i].original().wetting()[1].surfaces()[0].material,
            ],
            state: sources[i]
                .initialize_at_suction(INITIAL_Q, [[1.05, 1.28], [1.19, 1.28]][i], &gate)
                .unwrap(),
            gap_column: column,
        })
        .collect()
}
fn owner<'a>(
    sources: [&'a SourceAsymmetricCircularContactLineBridge; 2],
    elastic: bool,
) -> MultipleCircularWettingModalSystem<'a> {
    MultipleCircularWettingModalSystem::new(
        network(elastic, 8e-18, 2e-10, 1.0),
        bridges(sources, elastic),
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn states(system: &MultipleCircularWettingModalSystem<'_>) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}

#[test]
fn coupled_bridge_material_updates_are_atomic_and_keep_native_history_and_cross_compliance() {
    use fs_contact::capillary_bridge::BridgeVolumeConstraint;
    use fs_matdb::QueryPoint;
    let gate = CancelGate::new_clock_free();
    let (material, interfaces) = common::temperature_cards(1.0);
    let source = common::temperature_source(&material, &interfaces, 300.0);
    let hot = QueryPoint::new().with("T", 320.0).unwrap();
    let mut system = owner([&source, &source], true);
    let mut reordered = owner([&source, &source], true);
    let initial = system.total_energy_j().unwrap();
    let stale = source
        .refresh_at_query_point(
            system.state(1).unwrap(),
            &material,
            interfaces.each_ref(),
            &hot,
            BridgeVolumeConstraint::FixedVolume,
            &gate,
        )
        .unwrap();
    let frame = *system.step(&[0.0; 6], &[0.0; 2]).unwrap();
    reordered.step(&[0.0; 6], &[0.0; 2]).unwrap();
    let mut loss = frame.network_dissipation_j + frame.contact_line_dissipation_j;
    let before = (
        states(&system),
        [system.state(0), system.state(1)],
        [system.gap_m(0), system.gap_m(1)],
        *system.last_frame().unwrap(),
        [
            *system.bridge_frame(0).unwrap(),
            *system.bridge_frame(1).unwrap(),
        ],
        system.compliance_m_per_n().to_vec(),
        system.samples_rendered(),
    );
    let refresh0 = source
        .refresh_at_query_point(
            system.state(0).unwrap(),
            &material,
            interfaces.each_ref(),
            &hot,
            BridgeVolumeConstraint::FixedVolume,
            &gate,
        )
        .unwrap();
    let refresh1 = source
        .refresh_at_query_point(
            system.state(1).unwrap(),
            &material,
            interfaces.each_ref(),
            &hot,
            BridgeVolumeConstraint::FixedVolume,
            &gate,
        )
        .unwrap();
    assert!(
        system
            .apply_source_refreshes(&[(0, &refresh0), (1, &stale)], &gate)
            .is_err()
    );
    assert!(
        system
            .apply_source_refreshes(&[(0, &refresh0), (0, &refresh0)], &gate)
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        system
            .apply_source_refreshes(&[(0, &refresh0), (1, &refresh1)], &cancelled)
            .is_err()
    );
    assert_eq!(
        (
            states(&system),
            [system.state(0), system.state(1)],
            [system.gap_m(0), system.gap_m(1)],
            *system.last_frame().unwrap(),
            [
                *system.bridge_frame(0).unwrap(),
                *system.bridge_frame(1).unwrap()
            ],
            system.compliance_m_per_n().to_vec(),
            system.samples_rendered()
        ),
        before
    );
    assert_eq!(system.source(0).unwrap().identity(), source.identity());
    assert_eq!(system.source(1).unwrap().identity(), source.identity());
    let energy_before = system.total_energy_j().unwrap();
    let work = system
        .apply_source_refreshes(&[(0, &refresh0), (1, &refresh1)], &gate)
        .unwrap();
    let other_work = reordered
        .apply_source_refreshes(&[(1, &refresh1), (0, &refresh0)], &gate)
        .unwrap();
    assert_eq!(work, other_work);
    assert_eq!(
        work,
        refresh0.parameter_work_j() + refresh1.parameter_work_j()
    );
    assert_eq!(states(&system), before.0);
    for i in 0..2 {
        assert_eq!(
            system.state(i).unwrap().wetted_radii_m(),
            before.1[i].unwrap().wetted_radii_m()
        );
        assert_eq!(system.gap_m(i), before.2[i]);
        assert_eq!(*system.bridge_frame(i).unwrap(), before.4[i]);
    }
    assert_eq!(*system.last_frame().unwrap(), before.3);
    assert_eq!(system.compliance_m_per_n(), before.5);
    assert_eq!(system.samples_rendered(), before.6);
    close(
        system.total_energy_j().unwrap() - energy_before,
        work,
        4e-18,
    );
    for _ in 0..4 {
        let a = *system.step(&[0.0; 6], &[0.0; 2]).unwrap();
        let b = *reordered.step(&[0.0; 6], &[0.0; 2]).unwrap();
        assert_eq!(a, b);
        assert_eq!(states(&system), states(&reordered));
        loss += a.network_dissipation_j + a.contact_line_dissipation_j;
        close(
            system.total_energy_j().unwrap() + loss - work,
            initial,
            5e-16,
        );
    }
    assert_eq!(
        system.bridge_frame(0).unwrap().bridge.regimes[0],
        ContactLineRegime::Advancing
    );

    // A local material change leaves the other physical bridge untouched.
    let warm = QueryPoint::new().with("T", 310.0).unwrap();
    let local = system
        .source(0)
        .unwrap()
        .refresh_at_query_point(
            system.state(0).unwrap(),
            &material,
            interfaces.each_ref(),
            &warm,
            BridgeVolumeConstraint::FixedVolume,
            &gate,
        )
        .unwrap();
    let other = (
        system.source(1).unwrap().identity(),
        system.state(1),
        system.gap_m(1),
    );
    let bodies = states(&system);
    system
        .apply_source_refreshes(&[(0, &local)], &gate)
        .unwrap();
    assert_eq!(
        (
            system.source(1).unwrap().identity(),
            system.state(1),
            system.gap_m(1)
        ),
        other
    );
    assert_eq!(states(&system), bodies);
    assert_eq!(system.samples_rendered(), 5);
}

#[test]
fn multiple_bridge_source_refresh_respects_the_combined_native_energy_budget() {
    use fs_contact::capillary_bridge::BridgeVolumeConstraint;
    use fs_matdb::QueryPoint;
    let gate = CancelGate::new_clock_free();
    let (material, interfaces) = common::temperature_cards(1.0);
    let source = common::temperature_source(&material, &interfaces, 300.0);
    let mut system = MultipleCircularWettingModalSystem::new(
        network(false, 1e-30, 1e-30, 1.0),
        bridges([&source, &source], false),
        controls(),
        &gate,
    )
    .unwrap();
    let hot = QueryPoint::new().with("T", 320.0).unwrap();
    let refresh0 = source
        .refresh_at_query_point(
            system.state(0).unwrap(),
            &material,
            interfaces.each_ref(),
            &hot,
            BridgeVolumeConstraint::FixedVolume,
            &gate,
        )
        .unwrap();
    let refresh1 = source
        .refresh_at_query_point(
            system.state(1).unwrap(),
            &material,
            interfaces.each_ref(),
            &hot,
            BridgeVolumeConstraint::FixedVolume,
            &gate,
        )
        .unwrap();
    let old = (
        states(&system),
        system.state(0),
        system.state(1),
        system.total_energy_j().unwrap(),
    );
    assert!(matches!(
        system.apply_source_refreshes(&[(0, &refresh0), (1, &refresh1)], &gate),
        Err(MultipleCircularWettingModalError::Native(
            CapillaryModalError::EnergyAccuracy { .. }
        ))
    ));
    assert_eq!(
        (
            states(&system),
            system.state(0),
            system.state(1),
            system.total_energy_j().unwrap()
        ),
        old
    );
    assert_eq!(system.source(0).unwrap().identity(), source.identity());
    assert_eq!(system.source(1).unwrap().identity(), source.identity());
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.last_frame().is_none());
}
fn sources() -> [SourceAsymmetricCircularContactLineBridge; 2] {
    [
        source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]),
        source([1.05, 1.15], [3.0, 5.0], [1.0, 1.1], [1.1, 1.2]),
    ]
}

#[test]
fn shared_middle_body_couples_both_original_bridges_and_preserves_work_and_momentum() {
    let sources = sources();
    let mut system = owner([&sources[0], &sources[1]], false);
    let s = system.compliance_m_per_n();
    let dt = system.sample_period_s();
    close(s[1], -0.5 * dt * dt / MASS[1], 1e-20);
    close(s[1], s[2], 1e-20);
    let initial = system.total_energy_j().unwrap();
    let gaps = [system.gap_m(0).unwrap(), system.gap_m(1).unwrap()];
    let pinned_radius = system.state(0).unwrap().wetted_radii_m()[0];
    let mut loss = 0.0;
    for _ in 0..8 {
        let previous = [system.state(0).unwrap(), system.state(1).unwrap()];
        let old = states(&system);
        let frame = *system.step(&[0.0; 3], &[0.0; 2]).unwrap();
        let next = states(&system);
        let forces = [
            system.bridge_frame(0).unwrap().applied_mean_force_n,
            system.bridge_frame(1).unwrap().applied_mean_force_n,
        ];
        for i in 0..3 {
            let force = [-forces[0], forces[0] - forces[1], forces[1]][i];
            close(
                det::sqrt(MASS[i])
                    * (next[i].velocity_m_sqrt_kg_per_s - old[i].velocity_m_sqrt_kg_per_s),
                force * dt,
                3e-18,
            );
            close(
                det::sqrt(MASS[i])
                    * (next[i].displacement_m_sqrt_kg - old[i].displacement_m_sqrt_kg),
                det::sqrt(MASS[i]) * old[i].velocity_m_sqrt_kg_per_s * dt + 0.5 * force * dt * dt,
                3e-22,
            );
        }
        close(
            next.iter()
                .zip(MASS)
                .map(|(s, m)| det::sqrt(m) * s.velocity_m_sqrt_kg_per_s)
                .sum(),
            0.0,
            5e-18,
        );
        for i in 0..2 {
            let actual_gap = gaps[i] + next[i + 1].displacement_m_sqrt_kg / det::sqrt(MASS[i + 1])
                - next[i].displacement_m_sqrt_kg / det::sqrt(MASS[i]);
            close(system.gap_m(i).unwrap(), actual_gap, 2e-18);
            let point = system.bridge_frame(i).unwrap();
            assert_eq!(
                point.bridge,
                sources[i]
                    .advance(
                        previous[i],
                        system.gap_m(i).unwrap(),
                        dt,
                        &CancelGate::new_clock_free()
                    )
                    .unwrap()
            );
            assert!(
                point.compatibility_residual_m.abs() + point.estimated_compatibility_error_m
                    <= point.compatibility_tolerance_m
            );
            assert!(
                point.constitutive_residual_n.abs() + point.estimated_constitutive_force_error_n
                    <= point.force_tolerance_n
            );
        }
        assert_eq!(
            system.bridge_frame(0).unwrap().bridge.regimes,
            [ContactLineRegime::Pinned, ContactLineRegime::Advancing]
        );
        assert_eq!(system.state(0).unwrap().wetted_radii_m()[0], pinned_radius);
        assert_eq!(
            system.bridge_frame(1).unwrap().bridge.regimes,
            [ContactLineRegime::Advancing; 2]
        );
        assert!(
            frame.energy_residual_j.abs() + frame.estimated_constitutive_energy_error_j
                <= frame.energy_tolerance_j
        );
        assert!(frame.evaluations <= controls().maximum_evaluations);
        loss += frame.contact_line_dissipation_j;
        close(system.total_energy_j().unwrap() + loss, initial, 2e-16);
    }
}

#[test]
fn complete_elastic_network_preserves_original_modes_pressure_connections_and_prescribed_work() {
    let sources = sources();
    let mut system = owner([&sources[0], &sources[1]], true);
    let initial = system.total_energy_j().unwrap();
    let mut loss = 0.0;
    let mut work = 0.0;
    for k in 0..8 {
        let external = [2e-6, -1e-6, 0.0, 2e-6, -1e-6, 3e-6];
        let rates = [1e-5, -2e-5];
        let old = states(&system);
        let frame = *system.step(&external, &rates).unwrap();
        let next = states(&system);
        let expected_external: f64 = external
            .iter()
            .zip(next.iter().zip(&old))
            .map(|(f, (n, o))| f * (n.displacement_m_sqrt_kg - o.displacement_m_sqrt_kg))
            .sum();
        close(frame.external_work_j, expected_external, 1e-25);
        let expected_prescribed: f64 = (0..2)
            .map(|i| {
                -system.bridge_frame(i).unwrap().applied_mean_force_n
                    * rates[i]
                    * system.sample_period_s()
            })
            .sum();
        close(frame.prescribed_work_j, expected_prescribed, 1e-25);
        assert_ne!(frame.observer_pressure_pa, 0.0);
        assert_eq!(system.mode_count(), 6);
        let native = system.native_frame().unwrap();
        assert_eq!(native.connection_forces_n.len(), 2);
        assert!(native.connection_forces_n.iter().any(|x| *x != 0.0));
        assert!(native.component_dissipation_j > 0.0);
        assert!(native.connection_dissipation_j > 0.0);
        assert_eq!(native.numerical_dissipation_j, 0.0);
        assert_eq!(frame.sample, k + 1);
        loss += frame.network_dissipation_j + frame.contact_line_dissipation_j;
        work += frame.external_work_j + frame.prescribed_work_j;
        close(
            system.total_energy_j().unwrap() + loss - work,
            initial,
            2e-16,
        );
    }
}

#[test]
fn one_bridge_reduces_to_existing_original_owner() {
    let sources = sources();
    let mut binding = bridges([&sources[0], &sources[1]], false);
    binding.truncate(1);
    let state = binding[0].state;
    let mut single =
        AsymmetricCircularContactLineCapillaryModalSystem::new_asymmetric_circular_contact_line(
            network(false, 8e-18, 2e-10, 1.0),
            &sources[0],
            binding[0].substrate_materials,
            state,
            binding[0].gap_column.clone(),
            single_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut multiple = MultipleCircularWettingModalSystem::new(
        network(false, 8e-18, 2e-10, 1.0),
        binding,
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for _ in 0..8 {
        let a = *single.step(&[0.0; 3], 1e-5).unwrap();
        let b = *multiple.step(&[0.0; 3], &[1e-5]).unwrap();
        close(single.gap_m(), multiple.gap_m(0).unwrap(), 3e-15);
        close(
            a.applied_mean_force_n,
            multiple.bridge_frame(0).unwrap().applied_mean_force_n,
            6e-10,
        );
        close(a.surface_energy_j, b.surface_energy_j, 1e-16);
        for (x, y) in single.components().iter().zip(multiple.components()) {
            for (s, t) in x.states().iter().zip(y.states()) {
                close(s.displacement_m_sqrt_kg, t.displacement_m_sqrt_kg, 1e-16);
                close(
                    s.velocity_m_sqrt_kg_per_s,
                    t.velocity_m_sqrt_kg_per_s,
                    1e-11,
                );
            }
        }
    }
}

#[test]
fn joint_solution_is_covariant_under_bridge_permutation() {
    let sources = sources();
    let mut a = owner([&sources[0], &sources[1]], true);
    let mut bindings = bridges([&sources[0], &sources[1]], true);
    bindings.reverse();
    let mut b = MultipleCircularWettingModalSystem::new(
        network(true, 8e-18, 2e-10, 1.0),
        bindings,
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for _ in 0..6 {
        a.step(&[0.0; 6], &[1e-5, -2e-5]).unwrap();
        b.step(&[0.0; 6], &[-2e-5, 1e-5]).unwrap();
        for (x, y) in states(&a).iter().zip(states(&b)) {
            close(x.displacement_m_sqrt_kg, y.displacement_m_sqrt_kg, 2e-16);
            close(
                x.velocity_m_sqrt_kg_per_s,
                y.velocity_m_sqrt_kg_per_s,
                1e-11,
            );
        }
        for i in 0..2 {
            close(a.gap_m(i).unwrap(), b.gap_m(1 - i).unwrap(), 3e-15);
            for (x, y) in a
                .state(i)
                .unwrap()
                .wetted_radii_m()
                .iter()
                .zip(b.state(1 - i).unwrap().wetted_radii_m())
            {
                close(*x, y, 3e-15);
            }
        }
    }
}

#[test]
fn changing_one_source_changes_shared_motion_and_the_other_bridge_traction() {
    let sources = sources();
    let changed = source([1.05, 1.15], [0.3, 0.5], [1.0, 1.1], [1.1, 1.2]);
    let mut a = owner([&sources[0], &sources[1]], false);
    let mut b = owner([&sources[0], &changed], false);
    for _ in 0..8 {
        a.step(&[0.0; 3], &[0.0; 2]).unwrap();
        b.step(&[0.0; 3], &[0.0; 2]).unwrap();
    }
    assert!(
        (a.state(1).unwrap().wetted_radii_m()[1] - b.state(1).unwrap().wetted_radii_m()[1]).abs()
            > 1e-10
    );
    assert!((a.gap_m(0).unwrap() - b.gap_m(0).unwrap()).abs() > 1e-14);
    assert!(
        (a.bridge_frame(0).unwrap().applied_mean_force_n
            - b.bridge_frame(0).unwrap().applied_mean_force_n)
            .abs()
            > 1e-13
    );
}

#[test]
fn source_order_work_budget_and_cancellation_refuse_without_partial_history() {
    let sources = sources();
    let mut bad = bridges([&sources[0], &sources[1]], false);
    bad[1].substrate_materials.swap(0, 1);
    assert!(
        MultipleCircularWettingModalSystem::new(
            network(false, 8e-18, 2e-10, 1.0),
            bad,
            controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let mut a = owner([&sources[0], &sources[1]], true);
    let mut reference = owner([&sources[0], &sources[1]], true);
    a.step(&[0.0; 6], &[0.0; 2]).unwrap();
    reference.step(&[0.0; 6], &[0.0; 2]).unwrap();
    let old = (
        states(&a),
        a.state(0),
        a.state(1),
        a.gap_m(0),
        a.gap_m(1),
        *a.last_frame().unwrap(),
        *a.bridge_frame(0).unwrap(),
        *a.bridge_frame(1).unwrap(),
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(a.step_under_gate(&[0.0; 6], &[0.0; 2], &gate).is_err());
    // Only the second prescribed coordinate is outside its physical domain.
    assert!(a.step(&[0.0; 6], &[0.0, 1000.0]).is_err());
    assert_eq!(
        (
            states(&a),
            a.state(0),
            a.state(1),
            a.gap_m(0),
            a.gap_m(1),
            *a.last_frame().unwrap(),
            *a.bridge_frame(0).unwrap(),
            *a.bridge_frame(1).unwrap()
        ),
        old
    );
    assert_eq!(
        *a.step(&[0.0; 6], &[0.0; 2]).unwrap(),
        *reference.step(&[0.0; 6], &[0.0; 2]).unwrap()
    );
    assert_eq!(states(&a), states(&reference));
    let mut short = controls();
    short.maximum_evaluations = 4;
    let mut limited = MultipleCircularWettingModalSystem::new(
        network(false, 8e-18, 2e-10, 1.0),
        bridges([&sources[0], &sources[1]], false),
        short,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let old = (states(&limited), limited.state(0), limited.state(1));
    assert!(matches!(
        limited.step(&[0.0; 3], &[0.0; 2]),
        Err(MultipleCircularWettingModalError::Budget)
    ));
    assert_eq!((states(&limited), limited.state(0), limited.state(1)), old);
    assert_eq!(limited.samples_rendered(), 0);
    let mut conditioning = controls();
    conditioning.maximum_compatibility_gain = 0.1;
    let mut limited = MultipleCircularWettingModalSystem::new(
        network(false, 8e-18, 2e-10, 1.0),
        bridges([&sources[0], &sources[1]], false),
        conditioning,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let old = (states(&limited), limited.state(0), limited.state(1));
    assert!(matches!(
        limited.step(&[0.0; 3], &[0.0; 2]),
        Err(MultipleCircularWettingModalError::CompatibilityGain { .. })
    ));
    assert_eq!((states(&limited), limited.state(0), limited.state(1)), old);
    assert_eq!(limited.samples_rendered(), 0);
}

#[test]
fn original_shape_uncertainty_and_late_native_limits_keep_every_history_atomic() {
    let sources = sources();
    let mut ordinary = owner([&sources[0], &sources[1]], false);
    let frame = *ordinary.step(&[0.0; 3], &[0.0; 2]).unwrap();
    let mut strict = MultipleCircularWettingModalSystem::new(
        network(
            false,
            0.5 * frame.estimated_constitutive_energy_error_j,
            1e-30,
            1.0,
        ),
        bridges([&sources[0], &sources[1]], false),
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let old = (states(&strict), strict.state(0), strict.state(1));
    assert!(matches!(
        strict.step(&[0.0; 3], &[0.0; 2]),
        Err(MultipleCircularWettingModalError::Native(
            CapillaryModalError::EnergyAccuracy { .. }
        ))
    ));
    assert_eq!((states(&strict), strict.state(0), strict.state(1)), old);
    assert_eq!(strict.samples_rendered(), 0);
    assert!(strict.last_frame().is_none());
    assert!(strict.bridge_frame(0).is_none() && strict.bridge_frame(1).is_none());
    // Common acceleration changes no bridge gap, so failure occurs in original
    // native state admission after the simultaneous liquid solve succeeds.
    let mut limited = MultipleCircularWettingModalSystem::new(
        network(false, 8e-18, 2e-10, 1e-10),
        bridges([&sources[0], &sources[1]], false),
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let before = (states(&limited), limited.state(0), limited.state(1));
    let acceleration = MASS.map(|m| det::sqrt(m) * 1000.0);
    assert!(limited.step(&acceleration, &[0.0; 2]).is_err());
    assert_eq!(
        (states(&limited), limited.state(0), limited.state(1)),
        before
    );
    assert_eq!(limited.samples_rendered(), 0);
    assert_eq!(*limited.step(&[0.0; 3], &[0.0; 2]).unwrap(), frame);
}
