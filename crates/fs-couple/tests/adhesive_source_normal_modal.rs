//! Signed source adhesion on the original complete modal/Maxwell network.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;
#[path = "../../fs-contact/tests/normal_work_derivative/fixtures.rs"]
mod normal_fixtures;

use fs_contact::normal_patch::{
    NormalLawSource,
    adhesive::{JkrEnvelope, SourceJkrSphereContact},
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig,
    contact::source_normal::{
        SourceNormalModalConfig, SourceNormalModalError, SourceNormalModalSystem,
    },
};
use fs_couple::unilateral_contact::binding::BoundNormalContact;
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_phs::RelaxationBranch;

const MA: f64 = 2e-5;
const MB: f64 = 3e-5;
const RATE: u32 = 100_000;

fn source() -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    SourceJkrSphereContact::new(
        common::bind(&a, &b, None, 0.05).normal_state().clone(),
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}
fn column() -> Vec<f64> {
    vec![
        1.0 / det::sqrt(MA),
        0.25 / det::sqrt(MA),
        -1.0 / det::sqrt(MB),
        -0.18 / det::sqrt(MB),
    ]
}
fn controls(source: &SourceJkrSphereContact) -> SourceNormalModalConfig {
    let h = source.indentation_scale_m();
    SourceNormalModalConfig {
        maximum_iterations: 96,
        maximum_force_n: 0.05,
        minimum_approach_m: source.minimum_attached_indentation_m() + 1e-3 * h,
        maximum_approach_m: 8.0 * h,
        force_absolute_tolerance_n: 1e-13,
        force_relative_tolerance: 1e-10,
        coordinate_absolute_tolerance_m: 2e-10 * h,
        coordinate_relative_tolerance: 2e-10,
        contact_energy_absolute_tolerance_j: 1e-20,
        contact_energy_relative_tolerance: 5e-10,
    }
}
fn network(source: &SourceJkrSphereContact, memory: bool, velocity: f64) -> CoupledModalSystem {
    let h = source.indentation_scale_m();
    let models = [MA, MB]
        .into_iter()
        .enumerate()
        .map(|(i, mass)| {
            let mut model = ModalAcousticTimeModel::try_new(
                RATE,
                [1400.0 + 300.0 * i as f64, 2200.0 + 400.0 * i as f64]
                    .into_iter()
                    .map(|w| ModalAcousticMode {
                        angular_frequency_rad_s: w,
                        damping_ratio: 0.01,
                        pressure_per_modal_velocity: C64::ZERO,
                    })
                    .collect(),
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            let v = if i == 0 { MB } else { -MA } / (MA + MB) * velocity;
            model
                .restore_states(&[
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 0.0,
                        velocity_m_sqrt_kg_per_s: det::sqrt(mass) * v,
                    },
                    ModalAcousticState::default(),
                ])
                .unwrap();
            model
        })
        .collect();
    let attachment = |component, a, b| ModalAttachment {
        component,
        shapes: vec![a, b],
    };
    let connections = vec![
        ModalConnection {
            left: attachment(0, 1.0 / det::sqrt(MA), 0.12 / det::sqrt(MA)),
            right: attachment(1, 1.0 / det::sqrt(MB), -0.10 / det::sqrt(MB)),
            stiffness_n_m: 180.0,
            damping_n_s_m: 0.003,
            rest_extension_m: 0.0,
        },
        ModalConnection {
            left: attachment(0, 0.25 / det::sqrt(MA), 0.5 / det::sqrt(MA)),
            right: attachment(1, 0.08 / det::sqrt(MB), 0.35 / det::sqrt(MB)),
            stiffness_n_m: 300.0,
            damping_n_s_m: 0.002,
            rest_extension_m: 0.05 * h,
        },
    ];
    let gate = CancelGate::new_clock_free();
    let network = CoupledModalSystem::new(
        models,
        connections,
        ModalCouplingConfig {
            max_modes: 4,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e-4,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: 1.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-20,
            energy_relative_tolerance: 5e-10,
        },
        &gate,
    )
    .unwrap();
    if memory {
        network
            .with_material_relaxation(
                vec![RelaxationBranch {
                    projection: column().into_iter().flat_map(|b| [b, 0.0]).collect(),
                    stiffness: 2500.0,
                    relaxation_time_s: 3.0 / RATE as f64,
                }],
                Some(vec![det::sqrt(2500.0) * 0.5 * h]),
                &gate,
            )
            .unwrap()
    } else {
        network
    }
}
fn owner(
    source: &SourceJkrSphereContact,
    memory: bool,
    velocity: f64,
    approach: f64,
) -> SourceNormalModalSystem {
    SourceNormalModalSystem::new(
        network(source, memory, velocity),
        BoundNormalContact::new(
            NormalLawSource::JkrAttached(Box::new(source.clone())),
            column(),
            -approach,
        )
        .unwrap(),
        controls(source),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn states(network: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    network
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:e} != {expected:e}; normalized error {:e}",
        (actual - expected).abs() / scale
    );
}

#[test]
fn tensile_source_drives_the_original_complete_modal_network_and_work_ports() {
    let source = source();
    let h = source.indentation_scale_m();
    let mut owner = owner(&source, false, 0.0, -0.05 * h);
    let mut reference = network(&source, false, 0.0);
    let initial = owner.total_energy_j().unwrap();
    let before = states(owner.network());
    let external = [1e-6, -2e-7, 3e-7, 1e-7];
    let frame = owner.step(&external).unwrap().clone();
    assert!(frame.normal_force_n < 0.0);
    assert!(frame.endpoint_normal_force_n < 0.0);
    assert!(frame.approach_after_m > frame.approach_before_m);
    let original = source
        .attached_path_average(frame.approach_before_m, frame.approach_after_m)
        .unwrap();
    close(frame.normal_force_n, original, 0.01, 2e-10);

    // Independent original native propagation with the complete reciprocal row.
    let held: Vec<_> = external
        .iter()
        .zip(column())
        .map(|(f, b)| f - b * frame.normal_force_n)
        .collect();
    let native = reference.step(&held).unwrap().clone();
    for (actual, expected) in states(owner.network()).iter().zip(states(&reference)) {
        close(
            actual.displacement_m_sqrt_kg,
            expected.displacement_m_sqrt_kg,
            h * det::sqrt(MA),
            2e-10,
        );
        close(
            actual.velocity_m_sqrt_kg_per_s,
            expected.velocity_m_sqrt_kg_per_s,
            h * det::sqrt(MA) * RATE as f64,
            2e-10,
        );
    }
    let after = states(owner.network());
    let work: f64 = external
        .iter()
        .zip(before.iter().zip(&after))
        .map(|(f, (a, b))| f * (b.displacement_m_sqrt_kg - a.displacement_m_sqrt_kg))
        .sum();
    let contact_work = -frame.normal_force_n * (frame.approach_after_m - frame.approach_before_m);
    close(frame.external_work_j, work, 1e-9, 1e-10);
    close(frame.contact_work_on_network_j, contact_work, 1e-9, 1e-10);
    close(
        frame.network_dissipation_j,
        native.total_dissipation_j(),
        1e-9,
        1e-10,
    );
    close(
        owner.total_energy_j().unwrap() + frame.network_dissipation_j - frame.external_work_j,
        initial,
        1e-9,
        1e-9,
    );
    assert!(native.component_dissipation_j > 0.0);
    assert!(native.connection_dissipation_j > 0.0);
    assert_eq!(native.connection_forces_n.len(), 2);
    // The physical translation entries have equal/opposite contact force even
    // while all elastic modes and both independent connection rows participate.
    close(
        -column()[0] * frame.normal_force_n * det::sqrt(MA),
        column()[2] * frame.normal_force_n * det::sqrt(MB),
        0.01,
        1e-14,
    );
    assert!(frame.contact_work_residual_j.abs() <= frame.contact_work_tolerance_j);
}

#[test]
fn original_maxwell_history_and_independent_relaxation_loss_survive_adhesion() {
    let source = source();
    let h = source.indentation_scale_m();
    let mut owner = owner(&source, true, 0.0, -0.05 * h);
    let mut reference = network(&source, true, 0.0);
    let initial_memory = owner.network().material_memory().to_vec();
    let initial = owner.total_energy_j().unwrap();
    let mut total_loss = 0.0;
    let mut material_loss = 0.0;
    for _ in 0..20 {
        let frame = owner.step(&[0.0; 4]).unwrap().clone();
        let held: Vec<_> = column().iter().map(|b| -b * frame.normal_force_n).collect();
        let native = reference.step(&held).unwrap().clone();
        total_loss += frame.network_dissipation_j;
        material_loss += native.material_dissipation_j;
        assert_eq!(native.numerical_dissipation_j, 0.0);
        close(
            owner.network().material_energy_j().unwrap(),
            reference.material_energy_j().unwrap(),
            1e-9,
            1e-10,
        );
        for (actual, expected) in owner
            .network()
            .material_memory()
            .iter()
            .zip(reference.material_memory())
        {
            close(*actual, *expected, det::sqrt(2500.0) * h, 2e-10);
        }
        for (actual, expected) in states(owner.network()).iter().zip(states(&reference)) {
            close(
                actual.displacement_m_sqrt_kg,
                expected.displacement_m_sqrt_kg,
                h * det::sqrt(MA),
                2e-9,
            );
        }
        close(
            owner.total_energy_j().unwrap() + total_loss,
            initial,
            1e-9,
            2e-9,
        );
    }
    assert!(material_loss > 0.0);
    assert_ne!(owner.network().material_memory(), initial_memory);
    assert_eq!(owner.network().material_relaxation_branches().len(), 1);
    assert!(owner.network().last_frame().unwrap().material_dissipation_j > 0.0);
}

#[test]
fn attached_fold_and_receiving_refusals_preserve_original_modal_and_material_state() {
    let source = source();
    let h = source.indentation_scale_m();
    let mut failed = owner(
        &source,
        true,
        -0.1,
        source.minimum_attached_indentation_m() + 0.02 * h,
    );
    let initial = states(failed.network());
    let memory = failed.network().material_memory().to_vec();
    let energy = failed.total_energy_j().unwrap();
    assert!(failed.step(&[0.0; 4]).is_err());
    assert_eq!(states(failed.network()), initial);
    assert_eq!(failed.network().material_memory(), memory);
    assert_eq!(failed.total_energy_j().unwrap(), energy);
    assert_eq!(failed.samples_rendered(), 0);
    assert!(failed.last_frame().is_none());

    // This refusal occurs after the original native body/memory candidate and
    // contact work have all passed, exercising the final publication boundary.
    let mut owner = owner(&source, true, 0.0, -0.05 * h);
    owner.step(&[0.0; 4]).unwrap();
    let initial = states(owner.network());
    let memory = owner.network().material_memory().to_vec();
    let energy = owner.total_energy_j().unwrap();
    let frame = owner.last_frame().copied();
    let native = owner.network().last_frame().cloned();
    let refused: Result<(), SourceNormalModalError> =
        owner.step_with_admission(&[0.0; 4], &CancelGate::new_clock_free(), |_, _| {
            Err(SourceNormalModalError::Unresolved(
                "finite receiver refused",
            ))
        });
    assert!(refused.is_err());
    assert_eq!(states(owner.network()), initial);
    assert_eq!(owner.network().material_memory(), memory);
    assert_eq!(owner.total_energy_j().unwrap(), energy);
    assert_eq!(owner.last_frame(), frame.as_ref());
    assert_eq!(owner.network().last_frame(), native.as_ref());
    assert_eq!(owner.samples_rendered(), 1);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(owner.step_under_gate(&[0.0; 4], &cancelled).is_err());
    assert_eq!(states(owner.network()), initial);
    assert_eq!(owner.network().material_memory(), memory);
    assert_eq!(owner.samples_rendered(), 1);
    owner.step(&[0.0; 4]).unwrap();
    assert_eq!(owner.samples_rendered(), 2);
}

#[test]
fn native_signed_work_preserves_hertz_table_onset_release_and_dmt_attraction() {
    let jkr = source();
    let h = jkr.indentation_scale_m();
    let mut config = controls(&jkr);
    config.minimum_approach_m = -2.0 * h;
    let mut cases = Vec::new();
    for (from, velocity) in [(-0.1 * h, 0.03), (0.1 * h, -0.03)] {
        cases.push((
            NormalLawSource::HertzPoint(Box::new(normal_fixtures::point(None))),
            from,
            velocity,
            true,
        ));
        // Original retained-table owner with explicitly manufactured nodes;
        // this checks its interpolation/work, not a finite-gap BEM producer.
        cases.push((
            NormalLawSource::FiniteGap(Box::new(normal_fixtures::curve(10.0 * h, 0.005))),
            from,
            velocity,
            true,
        ));
    }
    cases.push((
        NormalLawSource::DmtSphere(Box::new(normal_fixtures::dmt())),
        -0.25 * h,
        0.0,
        false,
    ));
    for (source, from, velocity, crosses_touch) in cases {
        let energy = |delta: f64| {
            source
                .energy_and_force(if source.supports_signed_approach() {
                    delta
                } else {
                    delta.max(0.0)
                })
                .unwrap()
                .0
        };
        let mut owner = SourceNormalModalSystem::new(
            network(&jkr, false, velocity),
            BoundNormalContact::new(source.clone(), column(), -from).unwrap(),
            config,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let initial = owner.total_energy_j().unwrap();
        let frame = *owner.step(&[0.0; 4]).unwrap();
        let increment = frame.approach_after_m - frame.approach_before_m;
        let expected = source
            .mean_force_between(frame.approach_before_m, frame.approach_after_m)
            .unwrap()
            .unwrap();
        close(frame.normal_force_n, expected, 0.01, 2e-10);
        close(
            energy(frame.approach_after_m) - energy(frame.approach_before_m),
            expected * increment,
            1e-9,
            2e-10,
        );
        close(
            frame.contact_work_on_network_j,
            -expected * increment,
            1e-9,
            2e-10,
        );
        close(
            owner.total_energy_j().unwrap() + frame.network_dissipation_j,
            initial,
            1e-8,
            1e-9,
        );
        if crosses_touch {
            assert!(frame.approach_before_m * frame.approach_after_m < 0.0);
            assert!(frame.normal_force_n > 0.0);
        } else {
            assert!(frame.normal_force_n < 0.0);
            assert!(frame.approach_after_m < 0.0);
            assert!(increment > 0.0);
        }
    }
    // A nonzero original contact loss cannot be silently discarded by this
    // conservative wrapper, even though the native network has its own losses.
    assert!(
        SourceNormalModalSystem::new(
            network(&jkr, false, 0.0),
            BoundNormalContact::new(
                NormalLawSource::HertzPoint(Box::new(normal_fixtures::point(Some(0.2)))),
                column(),
                0.0
            )
            .unwrap(),
            config,
            &CancelGate::new_clock_free(),
        )
        .is_err()
    );
}
