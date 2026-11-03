//! G1/G3/G4: original simultaneous Coulomb laws on continuing modal momenta.
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, FrictionalContactError,
    NormalVelocityLaw,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::impulse::impact::frictional::{
    ModalContactFrameTrace, ModalFrictionalImpactContact, ModalFrictionalImpactError,
    PreparedModalFrictionalImpacts,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 32,
        max_connections: 16,
        max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 1e5,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-9,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 3e-9 + 3e-9 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn model(rate: u32, velocities: &[f64], budget: ModalAcousticTimeBudget) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
        rate,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO
            };
            velocities.len()
        ],
        budget,
    )
    .unwrap();
    model
        .restore_states(
            &velocities
                .iter()
                .map(|&velocity| ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: velocity,
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
    model
}
fn body(rate: u32, mass: f64, velocity: [f64; 3]) -> ModalAcousticTimeModel {
    model(
        rate,
        &velocity.map(|v| v * mass.sqrt()),
        ModalAcousticTimeBudget::audible_reference(),
    )
}
fn states(system: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}
fn axes() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
fn trace(component: usize, mass: f64, axes: [[f64; 3]; 3]) -> ModalContactFrameTrace {
    ModalContactFrameTrace {
        component,
        shapes: axes.map(|axis| axis.map(|value| value / mass.sqrt()).to_vec()),
    }
}
fn contact(
    id: u64,
    first: ModalContactFrameTrace,
    second: ModalContactFrameTrace,
    mu: f64,
    restitution: f64,
) -> ModalFrictionalImpactContact {
    ModalFrictionalImpactContact {
        id,
        first,
        second,
        reference_gap_m: 0.0,
        activation_tolerance_m: 0.0,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        },
    }
}
fn touch(system: &CoupledModalSystem, contact: &mut ModalFrictionalImpactContact) {
    let mut displacement = 0.0;
    for (frame, sign) in [(&contact.first, 1.0), (&contact.second, -1.0)] {
        for (&shape, state) in frame.shapes[0]
            .iter()
            .zip(system.components()[frame.component].states())
        {
            displacement += sign * shape * state.displacement_m_sqrt_kg;
        }
    }
    contact.reference_gap_m = -displacement + 1e-12;
    contact.activation_tolerance_m = 2e-12;
}
fn prepare(
    system: &CoupledModalSystem,
    contacts: Vec<ModalFrictionalImpactContact>,
    gate: &CancelGate,
) -> PreparedModalFrictionalImpacts {
    PreparedModalFrictionalImpacts::new(system, contacts, FrictionalContactConfig::default(), gate)
        .unwrap()
}

#[test]
fn finite_mass_stick_and_slip_retain_circular_friction_under_tangent_rotation() {
    let gate = CancelGate::new_clock_free();
    for rate in [1024, 2048] {
        for (mu, tangent_impulse, loss) in [(0.25, -0.6, 3.375), (2.0, -3.2, 7.6)] {
            for angle in [0.0_f64, 0.61] {
                let (s, c) = angle.sin_cos();
                let frame_axes = [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]];
                let mut system = CoupledModalSystem::new(
                    vec![body(rate, 4.0, [-2.0, 4.0, 0.0]), body(rate, 1.0, [0.0; 3])],
                    vec![],
                    config(),
                    &gate,
                )
                .unwrap();
                system.step_under_gate(&[0.0; 6], &gate).unwrap();
                let before = states(&system);
                let old_frame = system.last_frame().unwrap().clone();
                let mut row = contact(
                    17,
                    trace(0, 4.0, frame_axes),
                    trace(1, 1.0, frame_axes),
                    mu,
                    0.5,
                );
                touch(&system, &mut row);
                let mut prepared = prepare(&system, vec![row], &gate);
                close(prepared.delassus_kg_inv()[0], 1.25);
                let receipt = system
                    .apply_frictional_impacts(&mut prepared, &gate)
                    .unwrap();
                assert_eq!(receipt.contact_ids, vec![17]);
                close(receipt.contacts.impulses[0][0], 2.4);
                close(receipt.contacts.impulses[0][1], tangent_impulse * c);
                close(receipt.contacts.impulses[0][2], -tangent_impulse * s);
                let after = states(&system);
                let impulse_world = [2.4, tangent_impulse, 0.0];
                for axis in 0..3 {
                    close(
                        after[axis].velocity_m_sqrt_kg_per_s,
                        before[axis].velocity_m_sqrt_kg_per_s + 0.5 * impulse_world[axis],
                    );
                    close(
                        after[3 + axis].velocity_m_sqrt_kg_per_s,
                        -impulse_world[axis],
                    );
                    close(
                        2.0 * after[axis].velocity_m_sqrt_kg_per_s
                            + after[3 + axis].velocity_m_sqrt_kg_per_s,
                        2.0 * before[axis].velocity_m_sqrt_kg_per_s,
                    );
                }
                close(receipt.contacts.dissipated_energy, loss);
                close(receipt.contact_work_j, -loss);
                close(receipt.kinetic_energy_change_j, -loss);
                close(receipt.energy_residual_j, 0.0);
                close(receipt.contact_work_residual_j, 0.0);
                assert!(receipt.actual_residuals.tangential_duality_gap <= 1e-9);
                assert_eq!(system.last_frame(), Some(&old_frame));
                assert_eq!(system.samples_rendered(), 1);
                for (old, new) in before.iter().zip(&after) {
                    assert_eq!(old.displacement_m_sqrt_kg, new.displacement_m_sqrt_kg);
                }
                system.step_under_gate(&[0.0; 6], &gate).unwrap();
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
fn off_center_spin_couples_normal_and_tangent_impulses_and_preserves_angular_momentum() {
    let gate = CancelGate::new_clock_free();
    let root_mass = 2.0_f64.sqrt();
    // Body A: mass 2, polar inertia 1, contact arm r=(1,-1,0).
    // Body B: finite mass 2 at that contact position, with translation modes.
    let mut system = CoupledModalSystem::new(
        vec![
            model(
                2000,
                &[-root_mass, 3.0 * root_mass, 0.0, 0.0],
                ModalAcousticTimeBudget::audible_reference(),
            ),
            body(2000, 2.0, [0.0; 3]),
        ],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    let first = ModalContactFrameTrace {
        component: 0,
        shapes: [
            vec![1.0 / root_mass, 0.0, 0.0, 1.0],
            vec![0.0, 1.0 / root_mass, 0.0, 1.0],
            vec![0.0, 0.0, 1.0 / root_mass, 0.0],
        ],
    };
    let mut prepared = prepare(
        &system,
        vec![contact(3, first, trace(1, 2.0, axes()), 0.4, 0.2)],
        &gate,
    );
    close(prepared.delassus_kg_inv()[0], 2.0);
    close(prepared.delassus_kg_inv()[1], 1.0);
    close(prepared.delassus_kg_inv()[4], 2.0);
    let receipt = system
        .apply_frictional_impacts(&mut prepared, &gate)
        .unwrap();
    // Joint sliding equations: 2 pn + pt = 1.2, pt = -0.4 pn.
    // A normal-first clip would use pn=0.6 and miss the rebound target.
    close(receipt.contacts.impulses[0][0], 0.75);
    close(receipt.contacts.impulses[0][1], -0.3);
    close(receipt.relative_velocities_after_m_s[0][0], 0.2);
    close(receipt.relative_velocities_after_m_s[0][1], 3.15);
    let after = states(&system);
    close(after[3].velocity_m_sqrt_kg_per_s, 0.45);
    close(after[0].velocity_m_sqrt_kg_per_s / root_mass, -0.625);
    close(after[1].velocity_m_sqrt_kg_per_s / root_mass, 2.85);
    close(after[4].velocity_m_sqrt_kg_per_s / root_mass, -0.375);
    close(after[5].velocity_m_sqrt_kg_per_s / root_mass, 0.15);
    // Central spin plus the orbital impulse at B's actual arm cancel.
    close(
        after[3].velocity_m_sqrt_kg_per_s
            + root_mass * after[5].velocity_m_sqrt_kg_per_s
            + root_mass * after[4].velocity_m_sqrt_kg_per_s,
        0.0,
    );
    close(
        root_mass * (after[0].velocity_m_sqrt_kg_per_s + after[4].velocity_m_sqrt_kg_per_s),
        -2.0,
    );
    close(
        root_mass * (after[1].velocity_m_sqrt_kg_per_s + after[5].velocity_m_sqrt_kg_per_s),
        6.0,
    );
    close(receipt.contacts.dissipated_energy, 1.2225);
    close(receipt.kinetic_energy_change_j, -1.2225);
}

#[test]
fn shared_finite_body_contacts_use_full_cross_frames_and_stable_id_reduction() {
    let gate = CancelGate::new_clock_free();
    let angle = 0.6_f64;
    let (s, c) = angle.sin_cos();
    let rotated = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
    let initial = [
        [-1.0, 2.25, 0.0],
        [(1.0 - 0.8 * c) / 10.0, (-0.25 - 0.8 * s) / 10.0, 0.02],
        [0.4 * c, 0.4 * s, -1.1],
    ];
    let masses = [1.0, 10.0, 2.0];
    let mut reference = None;
    for reverse in [false, true] {
        let mut system = CoupledModalSystem::new(
            (0..3).map(|i| body(2000, masses[i], initial[i])).collect(),
            vec![],
            config(),
            &gate,
        )
        .unwrap();
        let mut descriptions = vec![
            contact(31, trace(0, 1.0, axes()), trace(1, 10.0, axes()), 0.25, 0.0),
            contact(
                7,
                trace(1, 10.0, rotated),
                trace(2, 2.0, rotated),
                0.25,
                0.0,
            ),
        ];
        if reverse {
            descriptions.reverse();
        }
        let mut prepared = prepare(&system, descriptions, &gate);
        if !reverse {
            close(prepared.delassus_kg_inv()[4], s / 10.0);
        }
        let receipt = system
            .apply_frictional_impacts(&mut prepared, &gate)
            .unwrap();
        for (index, &id) in receipt.contact_ids.iter().enumerate() {
            let expected = if id == 31 {
                [1.0, -0.25, 0.0]
            } else {
                [0.8, 0.0, -0.2]
            };
            for axis in 0..3 {
                close(receipt.contacts.impulses[index][axis], expected[axis]);
            }
        }
        let after = states(&system);
        let expected_velocities = [[0.0, 2.0, 0.0], [0.0; 3], [0.0, 0.0, -1.0]];
        for body in 0..3 {
            for axis in 0..3 {
                close(
                    after[3 * body + axis].velocity_m_sqrt_kg_per_s / masses[body].sqrt(),
                    expected_velocities[body][axis],
                );
            }
        }
        for axis in 0..3 {
            close(
                (0..3)
                    .map(|body| {
                        masses[body].sqrt() * after[3 * body + axis].velocity_m_sqrt_kg_per_s
                    })
                    .sum(),
                [0.0, 2.0, -2.0][axis],
            );
        }
        let initial_energy: f64 = (0..3)
            .map(|body| {
                0.5 * masses[body] * initial[body].iter().map(|value| value * value).sum::<f64>()
            })
            .sum();
        close(receipt.contacts.dissipated_energy, initial_energy - 3.0);
        if let Some((old_states, old_work)) = &reference {
            assert_eq!(&after, old_states);
            assert_eq!(receipt.contact_work_j.to_bits(), *old_work);
        } else {
            reference = Some((after, receipt.contact_work_j.to_bits()));
        }
    }
}

#[test]
fn locked_tangents_keep_native_velocity_determinacy_and_invalid_frames_refuse() {
    let gate = CancelGate::new_clock_free();
    let mut system = CoupledModalSystem::new(
        vec![
            model(2000, &[-2.0], ModalAcousticTimeBudget::audible_reference()),
            model(2000, &[0.0], ModalAcousticTimeBudget::audible_reference()),
        ],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    let first = ModalContactFrameTrace {
        component: 0,
        shapes: [vec![0.5], vec![0.0], vec![0.0]],
    };
    let second = ModalContactFrameTrace {
        component: 1,
        shapes: [vec![1.0], vec![0.0], vec![0.0]],
    };
    let description = contact(5, first, second, 0.4, 0.5);
    let mut bad = description.clone();
    bad.first.shapes[1].clear();
    assert!(matches!(
        PreparedModalFrictionalImpacts::new(&system, vec![bad], Default::default(), &gate),
        Err(ModalFrictionalImpactError::Frame {
            index: 0,
            axis: 1,
            ..
        })
    ));
    assert!(matches!(
        PreparedModalFrictionalImpacts::new(
            &system,
            vec![description.clone(); 2],
            Default::default(),
            &gate
        ),
        Err(ModalFrictionalImpactError::Contact(
            FrictionalContactError::Invalid("duplicate contact ID")
        ))
    ));
    let mut bad = description.clone();
    bad.first.shapes[0][0] = 0.0;
    bad.second.shapes[0][0] = 0.0;
    assert!(matches!(
        PreparedModalFrictionalImpacts::new(&system, vec![bad], Default::default(), &gate),
        Err(ModalFrictionalImpactError::Coupling(
            ModalCouplingError::Invalid { .. }
        ))
    ));
    let mut prepared = prepare(&system, vec![description], &gate);
    assert_eq!(
        prepared.delassus_kg_inv(),
        &[1.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
    let receipt = system
        .apply_frictional_impacts(&mut prepared, &gate)
        .unwrap();
    let FrictionalContactDeterminacy::UniqueVelocity(witness) = &receipt.contacts.determinacy
    else {
        panic!("locked tangent traction remains a gauge, not a unique local load");
    };
    assert_eq!((witness.dimension, witness.rank), (3, 1));
    assert!(receipt.contacts.impulse_work_witness.is_some());
    close(receipt.contacts.impulses[0][0], 1.2);
    close(states(&system)[0].velocity_m_sqrt_kg_per_s, -1.4);
    close(states(&system)[1].velocity_m_sqrt_kg_per_s, -1.2);
    close(receipt.contacts.dissipated_energy, 0.3);
}

#[test]
fn energy_growth_cancel_penetration_and_late_state_failure_remain_atomic() {
    let gate = CancelGate::new_clock_free();
    let mut capped = ModalAcousticTimeBudget::audible_reference();
    capped.maximum_abs_velocity_m_sqrt_kg_per_s = 0.6;
    let mut system = CoupledModalSystem::new(
        vec![
            body(2000, 1.0, [-1.0, 4.0, 0.0]),
            model(2000, &[0.0; 3], capped),
        ],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    let before = states(&system);
    let description = contact(1, trace(0, 1.0, axes()), trace(1, 1.0, axes()), 2.0, 0.0);
    let mut prepared = prepare(&system, vec![description.clone()], &gate);
    assert!(matches!(
        system.apply_frictional_impacts(&mut prepared, &gate),
        Err(ModalFrictionalImpactError::Coupling(
            ModalCouplingError::Component { component: 1, .. }
        ))
    ));
    assert_eq!(states(&system), before);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.apply_frictional_impacts(&mut prepared, &cancelled),
        Err(ModalFrictionalImpactError::Coupling(
            ModalCouplingError::Cancelled
        ))
    ));
    let mut penetrating = description.clone();
    penetrating.reference_gap_m = -1e-15;
    penetrating.activation_tolerance_m = 1e-6;
    let mut penetrating = prepare(&system, vec![penetrating], &gate);
    assert!(matches!(
        system.apply_frictional_impacts(&mut penetrating, &gate),
        Err(ModalFrictionalImpactError::Contact(
            FrictionalContactError::PenetratingContact { index: 0, .. }
        ))
    ));
    let mut source = description.clone();
    source.friction_coefficient = 0.0;
    source.normal_law = NormalVelocityLaw::Target { speed: 2.0 };
    let mut source = prepare(&system, vec![source], &gate);
    assert!(matches!(
        system.apply_frictional_impacts(&mut source, &gate),
        Err(ModalFrictionalImpactError::Contact(
            FrictionalContactError::EnergyIncrease { .. }
        ))
    ));
    assert_eq!(states(&system), before);
    let mut empty = prepare(&system, vec![], &gate);
    assert!(
        system
            .apply_frictional_impacts(&mut empty, &gate)
            .unwrap()
            .contacts
            .impulses
            .is_empty()
    );
    assert_eq!(states(&system), before);
    let mut retry = description;
    retry.friction_coefficient = 0.5;
    let mut retry = prepare(&system, vec![retry], &gate);
    let receipt = system.apply_frictional_impacts(&mut retry, &gate).unwrap();
    close(receipt.contacts.impulses[0][0], 0.5);
    close(receipt.contacts.impulses[0][1], -0.25);
    let expected = [-0.5, 3.75, 0.0, -0.5, 0.25, 0.0];
    for (state, expected) in states(&system).iter().zip(expected) {
        close(state.velocity_m_sqrt_kg_per_s, expected);
    }
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.last_frame().is_none());
}

fn live_network() -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let mut receiver = ModalAcousticTimeModel::try_new(
        2000,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 40.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            };
            3
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    receiver
        .restore_states(&[
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.01,
                velocity_m_sqrt_kg_per_s: 0.0,
            },
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: 0.0,
            },
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: 0.0,
            },
        ])
        .unwrap();
    let mut projection = vec![0.0; 12];
    projection[6] = 1.0;
    CoupledModalSystem::new(
        vec![body(2000, 4.0, [-0.2, 0.4, 0.0]), receiver],
        vec![],
        config(),
        &gate,
    )
    .unwrap()
    .with_material_relaxation(
        vec![fs_phs::RelaxationBranch {
            projection,
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
        &[3],
        1,
        &gate,
    )
    .unwrap()
}

#[test]
fn coulomb_reset_preserves_material_acoustic_and_forced_continuation_histories() {
    let gate = CancelGate::new_clock_free();
    let mut system = live_network();
    let mut oracle = live_network();
    let held = [1e-3, -2e-4, 0.0, 0.0, 0.0, 0.0];
    for _ in 0..3 {
        system.step_under_gate(&held, &gate).unwrap();
        oracle.step_under_gate(&held, &gate).unwrap();
    }
    let before = states(&system);
    let material = system.material_memory().to_vec();
    let radiation = system.radiation_state().to_vec();
    let frame = system.last_frame().unwrap().clone();
    assert!(material.iter().any(|&value| value != 0.0));
    assert!(radiation.iter().any(|&value| value != 0.0));
    let mut row = contact(9, trace(0, 4.0, axes()), trace(1, 1.0, axes()), 0.3, 0.25);
    touch(&system, &mut row);
    let mut prepared = prepare(&system, vec![row], &gate);
    let receipt = system
        .apply_frictional_impacts(&mut prepared, &gate)
        .unwrap();
    // Independent diagonal mobility 5/4: normal J=-v_n and sliding Jt=-0.3J.
    let normal = before[3].velocity_m_sqrt_kg_per_s - 0.5 * before[0].velocity_m_sqrt_kg_per_s;
    let tangent = -0.3 * normal;
    close(receipt.contacts.impulses[0][0], normal);
    close(receipt.contacts.impulses[0][1], tangent);
    assert!(receipt.relative_velocities_after_m_s[0][1] > 0.0);
    oracle
        .apply_modal_impulse(
            &[0.5 * normal, 0.5 * tangent, 0.0, -normal, -tangent, 0.0],
            &gate,
        )
        .unwrap();
    for (old, new) in before.iter().zip(states(&system)) {
        assert_eq!(old.displacement_m_sqrt_kg, new.displacement_m_sqrt_kg);
    }
    assert_eq!(system.material_memory(), material);
    assert_eq!(system.radiation_state(), radiation);
    assert_eq!(system.last_frame(), Some(&frame));
    assert_eq!(system.samples_rendered(), 3);
    close(receipt.contact_work_j, receipt.kinetic_energy_change_j);
    for _ in 0..16 {
        let actual = system.step_under_gate(&held, &gate).unwrap();
        let expected = oracle.step_under_gate(&held, &gate).unwrap();
        close(actual.total_energy_j(), expected.total_energy_j());
        close(actual.external_work_j, expected.external_work_j);
        for (a, b) in states(&system).iter().zip(states(&oracle)) {
            close(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
            close(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
        }
    }
}
