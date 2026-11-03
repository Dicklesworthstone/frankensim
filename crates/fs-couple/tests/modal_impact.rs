//! G1/G3/G4: reciprocal finite-body impacts on continuing native modal/cavity states.
use fs_contact::impact::NormalImpactError;
use fs_couple::cavity_modal::shared::{
    CavitySurface, CylindricalCavityThermalOptions, SharedCavity, ViscousNeckOptions,
};
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::impulse::impact::{
    ModalImpactContact, ModalImpactError, PreparedModalImpact,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_scenario::{HelmholtzCavity, NeckEnd};

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
fn mass(
    rate: u32,
    mass: f64,
    velocity: f64,
    budget: ModalAcousticTimeBudget,
) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_free_mass(rate, mass, 0.0, velocity, budget).unwrap()
}
fn states(system: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-13 + 3e-10 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn contact(reference_gap_m: f64, restitution: f64) -> ModalImpactContact {
    ModalImpactContact {
        first: ModalAttachment {
            component: 0,
            shapes: vec![0.5],
        },
        second: ModalAttachment {
            component: 1,
            shapes: vec![1.0],
        },
        reference_gap_m,
        activation_tolerance_m: 0.0,
        restitution,
    }
}

#[test]
fn finite_masses_exchange_momentum_and_continue_without_a_force_pulse() {
    let gate = CancelGate::new_clock_free();
    for rate in [1024, 2048] {
        for restitution in [0.0, 0.4, 1.0] {
            let mut system = CoupledModalSystem::new(
                vec![
                    mass(
                        rate,
                        4.0,
                        -1.0,
                        ModalAcousticTimeBudget::audible_reference(),
                    ),
                    mass(rate, 1.0, 0.5, ModalAcousticTimeBudget::audible_reference()),
                ],
                vec![],
                config(),
                &gate,
            )
            .unwrap();
            // A real accepted sample precedes the collision. Reference gap
            // chooses the contact location reached by this free-flight state.
            system.step_under_gate(&[0.0, 0.0], &gate).unwrap();
            let before = states(&system);
            let frame = system.last_frame().unwrap().clone();
            let prepared =
                PreparedModalImpact::new(&system, contact(1.5 / rate as f64, restitution), &gate)
                    .unwrap();
            let receipt = system.apply_normal_impact(&prepared, &gate).unwrap();
            let after = states(&system);
            let impulse = 1.2 * (1.0 + restitution);
            assert!(receipt.contact.impacted);
            close(receipt.inverse_effective_mass_kg_inv, 1.25);
            close(receipt.contact.normal_impulse, impulse);
            close(
                after[0].velocity_m_sqrt_kg_per_s * 0.5,
                -1.0 + impulse / 4.0,
            );
            close(after[1].velocity_m_sqrt_kg_per_s, 0.5 - impulse);
            close(
                2.0 * after[0].velocity_m_sqrt_kg_per_s + after[1].velocity_m_sqrt_kg_per_s,
                -3.5,
            );
            close(receipt.relative_normal_speed_after_m_s, 1.5 * restitution);
            close(
                receipt.contact.dissipated_energy,
                0.9 * (1.0 - restitution * restitution),
            );
            close(
                receipt.kinetic_energy_change_j,
                -receipt.contact.dissipated_energy,
            );
            close(receipt.energy_residual_j, 0.0);
            assert_eq!(receipt.sample, 1);
            assert_eq!(system.samples_rendered(), 1);
            assert_eq!(system.last_frame(), Some(&frame));
            assert_eq!(
                after[0].displacement_m_sqrt_kg,
                before[0].displacement_m_sqrt_kg
            );
            assert_eq!(
                after[1].displacement_m_sqrt_kg,
                before[1].displacement_m_sqrt_kg
            );
            system.step_under_gate(&[0.0, 0.0], &gate).unwrap();
            let continued = states(&system);
            for i in 0..2 {
                close(
                    continued[i].velocity_m_sqrt_kg_per_s,
                    after[i].velocity_m_sqrt_kg_per_s,
                );
                close(
                    continued[i].displacement_m_sqrt_kg,
                    after[i].displacement_m_sqrt_kg
                        + after[i].velocity_m_sqrt_kg_per_s / rate as f64,
                );
            }
        }
    }
}

#[test]
fn overlapping_ports_use_full_relative_mobility_and_shared_mode_cross_terms() {
    let gate = CancelGate::new_clock_free();
    let mut component = ModalAcousticTimeModel::try_new_with_free_coordinates(
        1024,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO
            };
            2
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    component
        .restore_states(&[
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: 3.0,
            },
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: -1.0,
            },
        ])
        .unwrap();
    let mut system = CoupledModalSystem::new(vec![component], vec![], config(), &gate).unwrap();
    let prepared = PreparedModalImpact::new(
        &system,
        ModalImpactContact {
            first: ModalAttachment {
                component: 0,
                shapes: vec![1.0, 2.0],
            },
            second: ModalAttachment {
                component: 0,
                shapes: vec![1.0, 0.0],
            },
            reference_gap_m: 0.0,
            activation_tolerance_m: 0.0,
            restitution: 1.0,
        },
        &gate,
    )
    .unwrap();
    // ||[1,2]-[1,0]||²=4; adding two independent mobilities would give 6.
    close(prepared.inverse_effective_mass_kg_inv(), 4.0);
    let receipt = system.apply_normal_impact(&prepared, &gate).unwrap();
    close(receipt.contact.normal_impulse, 1.0);
    close(receipt.relative_normal_speed_before_m_s, -2.0);
    close(receipt.relative_normal_speed_after_m_s, 2.0);
    let actual = states(&system);
    close(actual[0].velocity_m_sqrt_kg_per_s, 3.0);
    close(actual[1].velocity_m_sqrt_kg_per_s, 1.0);
    close(receipt.kinetic_energy_change_j, 0.0);
}

#[test]
fn separating_inactive_penetrating_cancelled_and_late_refused_impacts_publish_nothing() {
    let gate = CancelGate::new_clock_free();
    let mut system = CoupledModalSystem::new(
        vec![
            mass(
                1024,
                4.0,
                -1.0,
                ModalAcousticTimeBudget::audible_reference(),
            ),
            mass(
                1024,
                1.0,
                0.0,
                ModalAcousticTimeBudget {
                    maximum_abs_velocity_m_sqrt_kg_per_s: 0.5,
                    ..ModalAcousticTimeBudget::audible_reference()
                },
            ),
        ],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    let before = states(&system);
    let inactive = PreparedModalImpact::new(&system, contact(0.1, 1.0), &gate).unwrap();
    let receipt = system.apply_normal_impact(&inactive, &gate).unwrap();
    assert!(!receipt.contact.active && !receipt.contact.impacted);
    let penetrating = PreparedModalImpact::new(&system, contact(-1e-12, 1.0), &gate).unwrap();
    assert!(matches!(
        system.apply_normal_impact(&penetrating, &gate),
        Err(ModalImpactError::Contact(_))
    ));
    let active = PreparedModalImpact::new(&system, contact(0.0, 1.0), &gate).unwrap();
    // The first candidate is legal; the second would exceed its ORIGINAL
    // velocity ceiling. Neither accepted component may publish this reset.
    assert!(matches!(
        system.apply_normal_impact(&active, &gate),
        Err(ModalImpactError::Coupling(ModalCouplingError::Component {
            component: 1,
            ..
        }))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.apply_normal_impact(&active, &cancelled),
        Err(ModalImpactError::Coupling(ModalCouplingError::Cancelled))
    ));
    assert_eq!(states(&system), before);
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.last_frame().is_none());
    let mut separating = contact(0.0, 1.0);
    separating.first.shapes[0] = -0.5;
    separating.second.shapes[0] = -1.0;
    let separating = PreparedModalImpact::new(&system, separating, &gate).unwrap();
    assert!(
        !system
            .apply_normal_impact(&separating, &gate)
            .unwrap()
            .contact
            .impacted
    );
    assert_eq!(states(&system), before);
    let same = ModalAttachment {
        component: 0,
        shapes: vec![0.5],
    };
    assert!(matches!(
        PreparedModalImpact::new(
            &system,
            ModalImpactContact {
                first: same.clone(),
                second: same,
                reference_gap_m: 0.0,
                activation_tolerance_m: 0.0,
                restitution: 1.0,
            },
            &gate
        ),
        Err(ModalImpactError::Contact(
            NormalImpactError::InvalidInverseEffectiveMass
        ))
    ));
}

fn cavity_with_striker() -> SharedCavity {
    let gate = CancelGate::new_clock_free();
    let source = HelmholtzModalCavity::new(
        HelmholtzCavity {
            volume_m3: 1e-7,
            neck_radius_m: 1e-4,
            neck_length_m: 0.1,
            inner_end: NeckEnd::Flanged,
            outer_end: NeckEnd::Flanged,
        },
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap(),
        CavityNeckLoss::ExteriorMouth,
    )
    .unwrap();
    SharedCavity::new_with_cylindrical_thermal_wall(
        source,
        vec![
            CavitySurface {
                model: mass(
                    48_000,
                    4.0,
                    -0.02,
                    ModalAcousticTimeBudget::audible_reference(),
                ),
                outward_volume_shapes: vec![0.0],
            },
            CavitySurface {
                model: ModalAcousticTimeModel::try_new(
                    48_000,
                    vec![ModalAcousticMode {
                        angular_frequency_rad_s: 500.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    }],
                    ModalAcousticTimeBudget::audible_reference(),
                )
                .unwrap(),
                outward_volume_shapes: vec![1e-4],
            },
        ],
        48_000,
        ModalAcousticState {
            displacement_m_sqrt_kg: 2e-7,
            velocity_m_sqrt_kg_per_s: 0.0,
        },
        ModalAcousticTimeBudget::audible_reference(),
        config(),
        Some(ViscousNeckOptions {
            poles: 4,
            band_rad_s: [10.0, 1000.0],
            maximum_relative_error: 1e-3,
        }),
        CylindricalCavityThermalOptions {
            radius_m: 0.001,
            poles: 4,
            band_rad_s: [10.0, 1000.0],
            maximum_relative_error: 1e-3,
        },
        Some(vec![0.0; 4]),
        &gate,
    )
    .unwrap()
    .with_exterior_radiation(
        &Model {
            ports: 1,
            poles: vec![Pole {
                omega: 700.0,
                zeta: 0.1,
                coupling: vec![2e6],
            }],
        },
        &[2],
        1,
        &gate,
    )
    .unwrap()
}

#[test]
fn finite_striker_impacts_a_loaded_cavity_without_resetting_its_live_histories() {
    let gate = CancelGate::new_clock_free();
    let mut cavity = cavity_with_striker();
    let mut oracle = cavity_with_striker();
    for _ in 0..8 {
        cavity.step(&[0.0, 1e-4], 0.03, &gate).unwrap();
        oracle.step(&[0.0, 1e-4], 0.03, &gate).unwrap();
    }
    let before = states(cavity.network());
    let memory = cavity.network().material_memory().to_vec();
    let radiation = cavity.network().radiation_state().to_vec();
    let frame = cavity.network().last_frame().unwrap().clone();
    let observation = cavity.observe().unwrap();
    assert!(memory.iter().any(|&x| x != 0.0));
    assert!(radiation.iter().any(|&x| x != 0.0));
    assert!(
        before[3..]
            .iter()
            .any(|x| x.velocity_m_sqrt_kg_per_s != 0.0)
    );
    let speed = before[0].velocity_m_sqrt_kg_per_s * 0.5 - before[1].velocity_m_sqrt_kg_per_s;
    assert!(speed < 0.0);
    let mut description = contact(
        before[1].displacement_m_sqrt_kg - 0.5 * before[0].displacement_m_sqrt_kg + 1e-12,
        0.3,
    );
    description.activation_tolerance_m = 2e-12;
    let prepared = PreparedModalImpact::new(cavity.network(), description, &gate).unwrap();
    let receipt = cavity.apply_surface_impact(&prepared, &gate).unwrap();
    let impulse = -1.3 * speed / 1.25;
    close(receipt.contact.normal_impulse, impulse);
    // An independently computed two-body impulse reaches exactly the same
    // continuing cavity states and subsequent pressure/memory trajectory.
    oracle
        .apply_modal_impulse(&[0.5 * impulse, -impulse], 0.0, &gate)
        .unwrap();
    let after = states(cavity.network());
    for (a, b) in after.iter().zip(&before) {
        assert_eq!(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
    }
    assert_eq!(&after[2..], &before[2..]);
    assert_eq!(cavity.network().material_memory(), memory);
    assert_eq!(cavity.network().radiation_state(), radiation);
    assert_eq!(cavity.network().last_frame(), Some(&frame));
    assert_eq!(cavity.network().samples_rendered(), 8);
    let observed = cavity.observe().unwrap();
    assert_eq!(
        observed.pressure_increment_pa,
        observation.pressure_increment_pa
    );
    assert_eq!(observed.volume_change_m3, observation.volume_change_m3);
    assert_eq!(
        observed.compression_energy_j,
        observation.compression_energy_j
    );
    close(
        2.0 * (after[0].velocity_m_sqrt_kg_per_s - before[0].velocity_m_sqrt_kg_per_s)
            + after[1].velocity_m_sqrt_kg_per_s
            - before[1].velocity_m_sqrt_kg_per_s,
        0.0,
    );
    close(receipt.relative_normal_speed_after_m_s, -0.3 * speed);
    close(receipt.energy_residual_j, 0.0);
    for (a, b) in after.iter().zip(states(oracle.network())) {
        close(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
    }
    for _ in 0..16 {
        let a = cavity.step(&[0.0, 1e-4], 0.03, &gate).unwrap();
        let b = oracle.step(&[0.0, 1e-4], 0.03, &gate).unwrap();
        close(a.total_energy_j(), b.total_energy_j());
        close(
            cavity.observe().unwrap().pressure_increment_pa,
            oracle.observe().unwrap().pressure_increment_pa,
        );
    }
    let forbidden = PreparedModalImpact::new(
        cavity.network(),
        ModalImpactContact {
            first: ModalAttachment {
                component: 0,
                shapes: vec![0.5],
            },
            second: ModalAttachment {
                component: 2,
                shapes: vec![1.0],
            },
            reference_gap_m: 0.0,
            activation_tolerance_m: 0.0,
            restitution: 1.0,
        },
        &gate,
    )
    .unwrap();
    let before_refusal = states(cavity.network());
    assert!(cavity.apply_surface_impact(&forbidden, &gate).is_err());
    assert_eq!(states(cavity.network()), before_refusal);
}
