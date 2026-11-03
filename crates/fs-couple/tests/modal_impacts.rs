//! G1/G3/G4: one joint physical impact on continuing original modal states.
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy, MultiImpactError};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::impulse::impact::ModalImpactContact;
use fs_couple::render::schedule::force::coupled::impulse::impact::multiple::{
    ModalMultiImpactError, PreparedModalImpacts,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;

#[path = "modal_impacts/cavity.rs"]
mod cavity;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 16,
        max_connections: 16,
        max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 1e5,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-13,
        energy_relative_tolerance: 1e-9,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-12 + 5e-10 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn states(system: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}
fn free_network(rate: u32, masses: &[f64], velocities: &[f64]) -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let models = masses
        .iter()
        .zip(velocities)
        .map(|(&mass, &velocity)| {
            ModalAcousticTimeModel::try_free_mass(
                rate,
                mass,
                0.0,
                velocity,
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap()
        })
        .collect();
    CoupledModalSystem::new(models, vec![], config(), &gate).unwrap()
}
fn contact(
    system: &CoupledModalSystem,
    left: usize,
    right: usize,
    masses: &[f64],
    restitution: f64,
) -> ModalImpactContact {
    let left_shape = 1.0 / masses[left].sqrt();
    let right_shape = 1.0 / masses[right].sqrt();
    let accepted = states(system);
    ModalImpactContact {
        first: ModalAttachment {
            component: right,
            shapes: vec![right_shape],
        },
        second: ModalAttachment {
            component: left,
            shapes: vec![left_shape],
        },
        reference_gap_m: left_shape * accepted[left].displacement_m_sqrt_kg
            - right_shape * accepted[right].displacement_m_sqrt_kg
            + 1e-12,
        activation_tolerance_m: 2e-12,
        restitution,
    }
}
fn chain(system: &CoupledModalSystem, masses: &[f64], restitution: f64) -> Vec<ModalImpactContact> {
    vec![
        contact(system, 0, 1, masses, restitution),
        contact(system, 1, 2, masses, restitution),
    ]
}

#[test]
fn simultaneous_finite_mass_chain_is_permutation_invariant_and_continues() {
    let gate = CancelGate::new_clock_free();
    let masses = [4.0, 1.0, 9.0];
    for rate in [1024, 2048] {
        for restitution in [0.0, 0.5, 1.0] {
            for reverse in [false, true] {
                let mut system = free_network(rate, &masses, &[1.0, 0.0, 0.0]);
                system.step_under_gate(&[0.0; 3], &gate).unwrap();
                let before = states(&system);
                let frame = system.last_frame().unwrap().clone();
                let mut contacts = chain(&system, &masses, restitution);
                if reverse {
                    contacts.reverse();
                }
                let mut prepared = PreparedModalImpacts::new(
                    &system,
                    contacts,
                    MultiImpactConfig::default(),
                    &gate,
                )
                .unwrap();
                assert_eq!(prepared.contact_count(), 2);
                close(prepared.delassus_kg_inv()[1], -1.0);
                close(prepared.delassus_kg_inv()[2], -1.0);
                let receipt = system.apply_normal_impacts(&mut prepared, &gate).unwrap();
                assert!(matches!(
                    receipt.contacts.determinacy,
                    MultiImpactDeterminacy::UniqueImpulses
                ));
                let expected_impulses = [
                    20.0 / 7.0 * (1.0 + restitution),
                    18.0 / 7.0 * (1.0 + restitution),
                ];
                for i in 0..2 {
                    close(
                        receipt.contacts.normal_impulses[i],
                        expected_impulses[if reverse { 1 - i } else { i }],
                    );
                }
                // W=[[5/4,-1],[-1,10/9]]. Both initially touching rows
                // participate, although the second initially has zero speed.
                // Sequential independent pair solves produce another motion.
                let expected = [
                    1.0 - 5.0 / 7.0 * (1.0 + restitution),
                    2.0 / 7.0 * (1.0 + restitution),
                    2.0 / 7.0 * (1.0 + restitution),
                ];
                let after = states(&system);
                for i in 0..3 {
                    close(
                        after[i].velocity_m_sqrt_kg_per_s / masses[i].sqrt(),
                        expected[i],
                    );
                    assert_eq!(
                        after[i].displacement_m_sqrt_kg,
                        before[i].displacement_m_sqrt_kg
                    );
                }
                close(
                    after
                        .iter()
                        .zip(masses)
                        .map(|(state, mass)| state.velocity_m_sqrt_kg_per_s * mass.sqrt())
                        .sum(),
                    4.0,
                );
                let loss = 10.0 / 7.0 * (1.0 - restitution * restitution);
                close(receipt.contacts.dissipated_energy, loss);
                close(receipt.kinetic_energy_change_j, -loss);
                close(receipt.contact_work_j, -loss);
                close(receipt.energy_residual_j, 0.0);
                close(receipt.contact_work_residual_j, 0.0);
                assert_eq!(receipt.sample, 1);
                assert_eq!(system.samples_rendered(), 1);
                assert_eq!(system.last_frame(), Some(&frame));
                system.step_under_gate(&[0.0; 3], &gate).unwrap();
                for (next, reset) in states(&system).iter().zip(&after) {
                    close(
                        next.velocity_m_sqrt_kg_per_s,
                        reset.velocity_m_sqrt_kg_per_s,
                    );
                    close(
                        next.displacement_m_sqrt_kg,
                        reset.displacement_m_sqrt_kg + reset.velocity_m_sqrt_kg_per_s / rate as f64,
                    );
                }
            }
        }
    }
}

#[test]
fn redundant_contacts_determine_motion_without_inventing_individual_loads() {
    let gate = CancelGate::new_clock_free();
    let masses = [4.0, 1.0];
    for restitution in [0.0, 0.5, 1.0] {
        let mut system = free_network(2000, &masses, &[1.0, 0.0]);
        let row = contact(&system, 0, 1, &masses, restitution);
        let mut prepared = PreparedModalImpacts::new(
            &system,
            vec![row.clone(), row],
            MultiImpactConfig::default(),
            &gate,
        )
        .unwrap();
        for &entry in prepared.delassus_kg_inv() {
            close(entry, 1.25);
        }
        let receipt = system.apply_normal_impacts(&mut prepared, &gate).unwrap();
        let MultiImpactDeterminacy::UniqueVelocity(witness) = &receipt.contacts.determinacy else {
            panic!("duplicate rows do not determine a unique local load allocation");
        };
        assert_eq!((witness.dimension, witness.rank), (2, 1));
        let total_impulse = 0.8 * (1.0 + restitution);
        close(receipt.contacts.normal_impulses.iter().sum(), total_impulse);
        let after = states(&system);
        close(
            after[0].velocity_m_sqrt_kg_per_s * 0.5,
            1.0 - total_impulse * 0.25,
        );
        close(after[1].velocity_m_sqrt_kg_per_s, total_impulse);
        close(
            2.0 * after[0].velocity_m_sqrt_kg_per_s + after[1].velocity_m_sqrt_kg_per_s,
            4.0,
        );
        close(
            receipt.contacts.dissipated_energy,
            0.4 * (1.0 - restitution * restitution),
        );
        close(receipt.energy_residual_j, 0.0);
    }
}

#[test]
fn active_separating_rows_and_inactive_cross_response_retain_global_admission() {
    let gate = CancelGate::new_clock_free();
    let masses = [1.0; 3];
    let mut system = free_network(2000, &masses, &[1.0, 0.0, 0.25]);
    let before = states(&system);
    let mut rebound = PreparedModalImpacts::new(
        &system,
        chain(&system, &masses, 1.0),
        MultiImpactConfig::default(),
        &gate,
    )
    .unwrap();
    // Original joint Newton targets require impulses [5/4,1/2], giving
    // velocities [-1/4,3/4,3/4] and an inadmissible 1/16 J energy gain.
    let error = system
        .apply_normal_impacts(&mut rebound, &gate)
        .unwrap_err();
    let ModalMultiImpactError::Contact(MultiImpactError::EnergyIncrease { work, .. }) = error
    else {
        panic!("wrong global passivity refusal: {error:?}");
    };
    close(work, 1.0 / 16.0);
    assert_eq!(states(&system), before);
    assert_eq!(system.samples_rendered(), 0);
    let mut inelastic = PreparedModalImpacts::new(
        &system,
        chain(&system, &masses, 0.0),
        MultiImpactConfig::default(),
        &gate,
    )
    .unwrap();
    let receipt = system.apply_normal_impacts(&mut inelastic, &gate).unwrap();
    close(receipt.contacts.normal_impulses[0], 7.0 / 12.0);
    close(receipt.contacts.normal_impulses[1], 1.0 / 6.0);
    for state in states(&system) {
        close(state.velocity_m_sqrt_kg_per_s, 5.0 / 12.0);
    }
    close(receipt.contacts.dissipated_energy, 13.0 / 48.0);

    let mut system = free_network(2000, &masses, &[1.0, 0.0, 0.0]);
    let mut contacts = chain(&system, &masses, 1.0);
    contacts[1].reference_gap_m = 0.1;
    let mut prepared =
        PreparedModalImpacts::new(&system, contacts, MultiImpactConfig::default(), &gate).unwrap();
    let receipt = system.apply_normal_impacts(&mut prepared, &gate).unwrap();
    assert_eq!(receipt.contacts.active_indices, vec![0]);
    assert_eq!(receipt.contacts.normal_impulses[1], 0.0);
    close(receipt.relative_normal_speeds_after_m_s[0], 1.0);
    close(receipt.relative_normal_speeds_after_m_s[1], -1.0);
    for (state, expected) in states(&system).iter().zip([0.0, 1.0, 0.0]) {
        close(state.velocity_m_sqrt_kg_per_s, expected);
    }
}

#[test]
fn refused_candidates_cancelled_calls_and_penetration_publish_nothing_then_retry() {
    let gate = CancelGate::new_clock_free();
    let masses = [1.0; 3];
    let mut capped = ModalAcousticTimeBudget::audible_reference();
    capped.maximum_abs_velocity_m_sqrt_kg_per_s = 0.5;
    let mut system = CoupledModalSystem::new(
        vec![
            ModalAcousticTimeModel::try_free_mass(
                2000,
                1.0,
                0.0,
                1.0,
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
            ModalAcousticTimeModel::try_free_mass(2000, 1.0, 0.0, 0.0, capped).unwrap(),
            ModalAcousticTimeModel::try_free_mass(2000, 1.0, 0.0, 0.0, capped).unwrap(),
        ],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    let before = states(&system);
    let mut prepared = PreparedModalImpacts::new(
        &system,
        chain(&system, &masses, 1.0),
        MultiImpactConfig::default(),
        &gate,
    )
    .unwrap();
    assert!(matches!(
        system.apply_normal_impacts(&mut prepared, &gate),
        Err(ModalMultiImpactError::Coupling(
            ModalCouplingError::Component { component: 1, .. }
        ))
    ));
    assert_eq!(states(&system), before);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.apply_normal_impacts(&mut prepared, &cancelled),
        Err(ModalMultiImpactError::Coupling(
            ModalCouplingError::Cancelled
        ))
    ));
    assert_eq!(states(&system), before);
    let mut contacts = chain(&system, &masses, 0.0);
    contacts[1].reference_gap_m = -1e-15;
    let mut penetrating =
        PreparedModalImpacts::new(&system, contacts, MultiImpactConfig::default(), &gate).unwrap();
    assert!(matches!(
        system.apply_normal_impacts(&mut penetrating, &gate),
        Err(ModalMultiImpactError::Contact(
            MultiImpactError::InvalidContact { index: 1, .. }
        ))
    ));
    assert_eq!(states(&system), before);
    let row = contact(&system, 0, 1, &masses, 0.0);
    assert!(matches!(
        PreparedModalImpacts::new(&system, vec![row; 65], MultiImpactConfig::default(), &gate),
        Err(ModalMultiImpactError::Contact(MultiImpactError::Invalid(_)))
    ));
    let mut empty =
        PreparedModalImpacts::new(&system, vec![], MultiImpactConfig::default(), &gate).unwrap();
    let noop = system.apply_normal_impacts(&mut empty, &gate).unwrap();
    assert!(noop.contacts.normal_impulses.is_empty());
    assert_eq!(states(&system), before);
    // The earlier legal first component was staged before component 1 refused.
    // A fresh complete inelastic event must not reuse that partial candidate.
    let mut retry = PreparedModalImpacts::new(
        &system,
        chain(&system, &masses, 0.0),
        MultiImpactConfig::default(),
        &gate,
    )
    .unwrap();
    system.apply_normal_impacts(&mut retry, &gate).unwrap();
    for state in states(&system) {
        close(state.velocity_m_sqrt_kg_per_s, 1.0 / 3.0);
    }
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.last_frame().is_none());
}

#[test]
fn shared_modal_momentum_cannot_be_lost_inside_a_large_kinetic_background() {
    let gate = CancelGate::new_clock_free();
    let large = 2.0_f64.powi(55);
    let budget = ModalAcousticTimeBudget {
        maximum_abs_velocity_m_sqrt_kg_per_s: 1e20,
        maximum_total_energy_j: 1e40,
        ..ModalAcousticTimeBudget::audible_reference()
    };
    let models = [large + 8.0, large, large]
        .into_iter()
        .map(|velocity| {
            ModalAcousticTimeModel::try_free_mass(2000, 1.0, 0.0, velocity, budget).unwrap()
        })
        .collect();
    let mut limits = config();
    limits.maximum_total_energy_j = 1e40;
    let mut system = CoupledModalSystem::new(models, vec![], limits, &gate).unwrap();
    let before = states(&system);
    let mut prepared = PreparedModalImpacts::new(
        &system,
        chain(&system, &[1.0; 3], 0.0),
        MultiImpactConfig::default(),
        &gate,
    )
    .unwrap();
    assert!(matches!(
        system.apply_normal_impacts(&mut prepared, &gate),
        Err(ModalMultiImpactError::Coupling(
            ModalCouplingError::ImpulseResidual { .. }
        ))
    ));
    assert_eq!(states(&system), before);
}

fn live_network() -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let mut middle = ModalAcousticTimeModel::try_new(
        2000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 40.0,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    middle
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.01,
            velocity_m_sqrt_kg_per_s: 0.0,
        }])
        .unwrap();
    CoupledModalSystem::new(
        vec![
            ModalAcousticTimeModel::try_free_mass(
                2000,
                4.0,
                0.0,
                0.5,
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
            middle,
            ModalAcousticTimeModel::try_free_mass(
                2000,
                9.0,
                0.0,
                -0.2,
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
        ],
        vec![],
        config(),
        &gate,
    )
    .unwrap()
    .with_material_relaxation(
        vec![fs_phs::RelaxationBranch {
            projection: vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            stiffness: 200.0,
            relaxation_time_s: 0.05,
        }],
        Some(vec![0.02]),
        &gate,
    )
    .unwrap()
    .with_radiation_load(
        &Model {
            ports: 1,
            poles: vec![Pole {
                omega: 60.0,
                zeta: 0.1,
                coupling: vec![10.0],
            }],
        },
        &[1],
        1,
        &gate,
    )
    .unwrap()
}

#[test]
fn simultaneous_impacts_preserve_live_material_radiation_and_continuing_dynamics() {
    let gate = CancelGate::new_clock_free();
    let mut system = live_network();
    let mut oracle = live_network();
    for _ in 0..3 {
        system.step_under_gate(&[0.0; 3], &gate).unwrap();
        oracle.step_under_gate(&[0.0; 3], &gate).unwrap();
    }
    let before = states(&system);
    let material = system.material_memory().to_vec();
    let radiation = system.radiation_state().to_vec();
    let frame = system.last_frame().unwrap().clone();
    assert!(material.iter().any(|&value| value != 0.0));
    assert!(radiation.iter().any(|&value| value != 0.0));
    let mut prepared = PreparedModalImpacts::new(
        &system,
        chain(&system, &[4.0, 1.0, 9.0], 0.25),
        MultiImpactConfig::default(),
        &gate,
    )
    .unwrap();
    let receipt = system.apply_normal_impacts(&mut prepared, &gate).unwrap();
    let impulses = &receipt.contacts.normal_impulses;
    oracle
        .apply_modal_impulse(
            &[
                -0.5 * impulses[0],
                impulses[0] - impulses[1],
                (1.0 / 3.0) * impulses[1],
            ],
            &gate,
        )
        .unwrap();
    for (after, old) in states(&system).iter().zip(&before) {
        assert_eq!(after.displacement_m_sqrt_kg, old.displacement_m_sqrt_kg);
    }
    assert_eq!(system.material_memory(), material);
    assert_eq!(system.radiation_state(), radiation);
    assert_eq!(system.last_frame(), Some(&frame));
    assert_eq!(system.samples_rendered(), 3);
    assert!(
        receipt
            .contacts
            .normal_impulses
            .iter()
            .all(|&value| value > 0.0)
    );
    close(receipt.contact_work_j, receipt.kinetic_energy_change_j);
    for _ in 0..16 {
        let actual = system.step_under_gate(&[0.0; 3], &gate).unwrap();
        let expected = oracle.step_under_gate(&[0.0; 3], &gate).unwrap();
        close(actual.total_energy_j(), expected.total_energy_j());
        for (a, b) in states(&system).iter().zip(states(&oracle)) {
            close(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
            close(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
        }
    }
}
