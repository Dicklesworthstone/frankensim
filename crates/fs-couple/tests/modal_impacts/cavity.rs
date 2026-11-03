//! Full shared-cavity integration, kept distinct from the modal-only tests.
use super::*;
use fs_couple::cavity_modal::shared::{
    CavitySurface, CylindricalCavityThermalOptions, SharedCavity, ViscousNeckOptions,
};
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::{HelmholtzCavity, NeckEnd};

fn cavity() -> SharedCavity {
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
    let surfaces = vec![
        CavitySurface {
            model: ModalAcousticTimeModel::try_free_mass(
                48_000,
                4.0,
                0.0,
                0.02,
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
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
        CavitySurface {
            model: ModalAcousticTimeModel::try_free_mass(
                48_000,
                9.0,
                0.0,
                -0.01,
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
            outward_volume_shapes: vec![0.0],
        },
    ];
    SharedCavity::new_with_cylindrical_thermal_wall(
        source,
        surfaces,
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
        &[3],
        1,
        &gate,
    )
    .unwrap()
}

#[test]
fn two_finite_strikers_share_one_loaded_cavity_impact_and_retain_every_history() {
    let gate = CancelGate::new_clock_free();
    let mut actual = cavity();
    let mut oracle = cavity();
    for _ in 0..8 {
        actual.step(&[0.0, 1e-4, 0.0], 0.03, &gate).unwrap();
        oracle.step(&[0.0, 1e-4, 0.0], 0.03, &gate).unwrap();
    }
    let before = states(actual.network());
    let material = actual.network().material_memory().to_vec();
    let radiation = actual.network().radiation_state().to_vec();
    let observation = actual.observe().unwrap();
    let frame = actual.network().last_frame().unwrap().clone();
    assert!(material.iter().any(|&value| value != 0.0));
    assert!(radiation.iter().any(|&value| value != 0.0));
    assert!(
        before[4..]
            .iter()
            .any(|state| state.velocity_m_sqrt_kg_per_s != 0.0)
    );
    let mut prepared = PreparedModalImpacts::new(
        actual.network(),
        chain(actual.network(), &[4.0, 1.0, 9.0], 0.25),
        MultiImpactConfig::default(),
        &gate,
    )
    .unwrap();
    let receipt = actual.apply_surface_impacts(&mut prepared, &gate).unwrap();
    // Independent inverse of [[5/4,-1],[-1,10/9]], applied to the actual
    // incoming physical speeds; neck/storage forces have no zero-time impulse.
    let rhs0 =
        1.25 * (0.5 * before[0].velocity_m_sqrt_kg_per_s - before[1].velocity_m_sqrt_kg_per_s);
    let rhs1 =
        1.25 * (before[1].velocity_m_sqrt_kg_per_s - before[2].velocity_m_sqrt_kg_per_s / 3.0);
    let j0 = ((10.0 / 9.0) * rhs0 + rhs1) / (7.0 / 18.0);
    let j1 = (rhs0 + 1.25 * rhs1) / (7.0 / 18.0);
    assert!(j0 > 0.0 && j1 > 0.0);
    close(receipt.contacts.normal_impulses[0], j0);
    close(receipt.contacts.normal_impulses[1], j1);
    oracle
        .apply_modal_impulse(&[-0.5 * j0, j0 - j1, j1 / 3.0], 0.0, &gate)
        .unwrap();
    let after = states(actual.network());
    for (a, b) in after.iter().zip(&before) {
        assert_eq!(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
    }
    assert_eq!(&after[3..], &before[3..]);
    assert_eq!(actual.network().material_memory(), material);
    assert_eq!(actual.network().radiation_state(), radiation);
    assert_eq!(actual.network().last_frame(), Some(&frame));
    assert_eq!(actual.network().samples_rendered(), 8);
    let after_observation = actual.observe().unwrap();
    assert_eq!(
        after_observation.pressure_increment_pa,
        observation.pressure_increment_pa
    );
    assert_eq!(
        after_observation.volume_change_m3,
        observation.volume_change_m3
    );
    assert_eq!(
        after_observation.compression_energy_j,
        observation.compression_energy_j
    );
    close(
        2.0 * (after[0].velocity_m_sqrt_kg_per_s - before[0].velocity_m_sqrt_kg_per_s)
            + after[1].velocity_m_sqrt_kg_per_s
            - before[1].velocity_m_sqrt_kg_per_s
            + 3.0 * (after[2].velocity_m_sqrt_kg_per_s - before[2].velocity_m_sqrt_kg_per_s),
        0.0,
    );
    for _ in 0..16 {
        let a = actual.step(&[0.0, 1e-4, 0.0], 0.03, &gate).unwrap();
        let b = oracle.step(&[0.0, 1e-4, 0.0], 0.03, &gate).unwrap();
        close(a.total_energy_j(), b.total_energy_j());
        close(
            actual.observe().unwrap().pressure_increment_pa,
            oracle.observe().unwrap().pressure_increment_pa,
        );
    }
    let mut forbidden = contact(actual.network(), 0, 1, &[4.0, 1.0, 9.0], 0.0);
    forbidden.first.component = 3; // Retained acoustic neck, not a physical surface.
    let mut forbidden = PreparedModalImpacts::new(
        actual.network(),
        vec![forbidden],
        MultiImpactConfig::default(),
        &gate,
    )
    .unwrap();
    let before_refusal = states(actual.network());
    assert!(actual.apply_surface_impacts(&mut forbidden, &gate).is_err());
    assert_eq!(states(actual.network()), before_refusal);
}
