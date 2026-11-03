//! G1/G3: physical pressure/momentum inputs on the continuing original cavity.

use fs_couple::cavity_modal::shared::{
    CavitySurface, CylindricalCavityThermalOptions, SharedCavity, ViscousNeckOptions,
};
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig, ModalCouplingError};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};
use fs_scenario::{HelmholtzCavity, NeckEnd};

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 8,
        max_setup_terms: 100_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 1.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-14,
        energy_relative_tolerance: 1e-8,
    }
}

fn source(small: bool, loss: CavityNeckLoss) -> HelmholtzModalCavity {
    HelmholtzModalCavity::new(
        HelmholtzCavity {
            volume_m3: if small { 1e-7 } else { 0.01 },
            neck_radius_m: if small { 1e-4 } else { 0.01 },
            neck_length_m: if small { 0.1 } else { 0.02 },
            inner_end: NeckEnd::Flanged,
            outer_end: NeckEnd::Flanged,
        },
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap(),
        loss,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 5e-13 + 2e-10 * actual.abs().max(expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn states(cavity: &SharedCavity) -> Vec<ModalAcousticState> {
    cavity
        .network()
        .components()
        .iter()
        .flat_map(|component| component.states().iter().copied())
        .collect()
}

#[test]
fn pressure_impulses_match_physical_momentum_work_and_rank_one_cavity_motion() {
    let gate = CancelGate::new_clock_free();
    let masses = [4.0_f64, 1.0];
    let areas = [0.01, -0.01];
    let pressure = [0.05, -0.04];
    let initial_velocity = [0.02, -0.03];
    let neck_pressure = 0.002;
    let neck_initial_velocity = 0.01;
    let mut prior_jump = None;
    for rate in [48_000, 96_000] {
        let gas = source(false, CavityNeckLoss::Resistance(0.0));
        let neck_mass = gas.neck_mass_kg();
        let neck_area = gas.neck_area_m2();
        let neck_shape = gas.neck_volume_shape_m2_per_sqrt_kg();
        let stiffness = gas.compression_stiffness_n_m();
        let surfaces = (0..2)
            .map(|i| CavitySurface {
                model: ModalAcousticTimeModel::try_free_mass(
                    rate,
                    masses[i],
                    0.0,
                    initial_velocity[i],
                    ModalAcousticTimeBudget::audible_reference(),
                )
                .unwrap(),
                outward_volume_shapes: vec![areas[i] / det::sqrt(masses[i])],
            })
            .collect();
        let mut cavity = SharedCavity::new(
            gas,
            surfaces,
            rate,
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: det::sqrt(neck_mass) * neck_initial_velocity,
            },
            ModalAcousticTimeBudget::audible_reference(),
            config(),
            &gate,
        )
        .unwrap();
        let initial_energy = cavity.network().total_energy_j().unwrap();
        let receipt = cavity
            .apply_pressure_impulses(&pressure, neck_pressure, &gate)
            .unwrap();
        let after = states(&cavity);
        assert!(after.iter().all(|s| s.displacement_m_sqrt_kg == 0.0));
        assert_eq!(cavity.network().samples_rendered(), 0);
        assert!(cavity.network().last_frame().is_none());
        assert_eq!(cavity.observe().unwrap().pressure_increment_pa, 0.0);
        let mut pressure_work = 0.0;
        for i in 0..2 {
            let velocity = after[i].velocity_m_sqrt_kg_per_s / det::sqrt(masses[i]);
            close(
                masses[i] * (velocity - initial_velocity[i]),
                -areas[i] * pressure[i],
            );
            pressure_work -= pressure[i] * areas[i] * f64::midpoint(initial_velocity[i], velocity);
        }
        let neck_velocity = after[2].velocity_m_sqrt_kg_per_s / det::sqrt(neck_mass);
        close(
            neck_mass * (neck_velocity - neck_initial_velocity),
            -neck_area * neck_pressure,
        );
        pressure_work -=
            neck_pressure * neck_area * f64::midpoint(neck_initial_velocity, neck_velocity);
        close(receipt.external_work_j, pressure_work);
        close(receipt.kinetic_energy_change_j, pressure_work);
        close(
            cavity.network().total_energy_j().unwrap() - initial_energy,
            pressure_work,
        );
        if let Some(prior) = prior_jump {
            assert_eq!(after, prior, "an impulse cannot depend on sample duration");
        }
        prior_jump = Some(after.clone());

        // Independent analytic solution of q'' + k b b^T q = 0. Only the
        // collective compression direction oscillates; its nullspace drifts.
        let column = [
            areas[0] / (neck_area * det::sqrt(masses[0])),
            areas[1] / (neck_area * det::sqrt(masses[1])),
            neck_shape / neck_area,
        ];
        let norm2 = column.iter().map(|b| b * b).sum::<f64>();
        let omega = det::sqrt(stiffness * norm2);
        let projection = column
            .iter()
            .zip(&after)
            .map(|(b, s)| b * s.velocity_m_sqrt_kg_per_s)
            .sum::<f64>()
            / norm2;
        let energy = cavity.network().total_energy_j().unwrap();
        for tick in 1..=rate / 750 {
            cavity.step(&[0.0; 2], 0.0, &gate).unwrap();
            let t = f64::from(tick) / f64::from(rate);
            for (i, state) in states(&cavity).iter().enumerate() {
                close(
                    state.displacement_m_sqrt_kg,
                    after[i].velocity_m_sqrt_kg_per_s * t
                        + column[i] * projection * (det::sin(omega * t) / omega - t),
                );
                close(
                    state.velocity_m_sqrt_kg_per_s,
                    after[i].velocity_m_sqrt_kg_per_s
                        + column[i] * projection * (det::cos(omega * t) - 1.0),
                );
            }
            close(cavity.network().total_energy_j().unwrap(), energy);
        }
        assert!(cavity.observe().unwrap().pressure_increment_pa.abs() > 1e-4);
    }
}

fn history_cavity() -> SharedCavity {
    let gate = CancelGate::new_clock_free();
    let surface = CavitySurface {
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
    };
    SharedCavity::new_with_cylindrical_thermal_wall(
        source(true, CavityNeckLoss::ExteriorMouth),
        vec![surface],
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
        &[1],
        1,
        &gate,
    )
    .unwrap()
}

#[test]
fn continuing_kicks_preserve_viscous_thermal_radiation_and_interval_histories() {
    let gate = CancelGate::new_clock_free();
    let mut cavity = history_cavity();
    let mut generalized = history_cavity();
    for _ in 0..8 {
        cavity.step(&[1e-4], 0.03, &gate).unwrap();
        generalized.step(&[1e-4], 0.03, &gate).unwrap();
    }
    let before = states(&cavity);
    let thermal = cavity.network().material_memory().to_vec();
    let radiation = cavity.network().radiation_state().to_vec();
    let frame = cavity.network().last_frame().unwrap().clone();
    let observation = cavity.observe().unwrap();
    let exterior = cavity.neck_exterior_pressure_pa().unwrap();
    let energy = cavity.network().total_energy_j().unwrap();
    assert!(
        before[2..]
            .iter()
            .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0)
    );
    assert!(thermal.iter().any(|x| *x != 0.0));
    assert!(radiation.iter().any(|x| *x != 0.0));
    let pressure = 2e-3;
    let receipt = cavity
        .apply_pressure_impulses(&[pressure], 2e-4, &gate)
        .unwrap();
    let modal = generalized
        .apply_modal_impulse(&[-1e-4 * pressure], 2e-4, &gate)
        .unwrap();
    assert_eq!(receipt, modal);
    let after = states(&cavity);
    assert_eq!(after, states(&generalized));
    assert_eq!(after[2..], before[2..]);
    for (a, b) in after.iter().zip(&before) {
        assert_eq!(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
    }
    assert_eq!(cavity.network().material_memory(), thermal);
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
    assert_eq!(cavity.neck_exterior_pressure_pa().unwrap(), exterior);
    close(
        cavity.network().total_energy_j().unwrap() - energy,
        receipt.external_work_j,
    );
    for _ in 0..16 {
        assert_eq!(
            cavity.step(&[1e-4], 0.03, &gate).unwrap(),
            generalized.step(&[1e-4], 0.03, &gate).unwrap()
        );
    }
}

#[test]
fn pressure_projection_respects_surface_groups_signed_modes_and_zero_participation() {
    let gate = CancelGate::new_clock_free();
    let make = |count| {
        let mut model = ModalAcousticTimeModel::try_new(
            48_000,
            vec![
                ModalAcousticMode {
                    angular_frequency_rad_s: 500.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                };
                count
            ],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap();
        model
            .restore_states(&vec![
                ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: 0.5,
                };
                count
            ])
            .unwrap();
        model
    };
    let mut cavity = SharedCavity::new(
        source(false, CavityNeckLoss::Resistance(0.0)),
        vec![
            CavitySurface {
                model: make(3),
                outward_volume_shapes: vec![0.25, -0.125, 0.0],
            },
            CavitySurface {
                model: make(1),
                outward_volume_shapes: vec![0.0625],
            },
        ],
        48_000,
        ModalAcousticState::default(),
        ModalAcousticTimeBudget::audible_reference(),
        config(),
        &gate,
    )
    .unwrap();
    let receipt = cavity
        .apply_pressure_impulses(&[0.5, -0.25], 0.0, &gate)
        .unwrap();
    let jumps = [-0.125, 0.0625, 0.0, 0.015625];
    for (state, jump) in states(&cavity).iter().zip(jumps) {
        assert_eq!(state.velocity_m_sqrt_kg_per_s, 0.5 + jump);
    }
    assert_eq!(states(&cavity)[4], ModalAcousticState::default());
    let work = jumps.iter().map(|j| 0.5 * j + 0.5 * j * j).sum::<f64>();
    assert!(work < 0.0);
    assert_eq!(receipt.external_work_j, work);
    assert_eq!(receipt.kinetic_energy_change_j, work);
}

#[test]
fn pressure_and_native_budget_refusals_preserve_the_entire_continuing_cavity() {
    let gate = CancelGate::new_clock_free();
    let mut cavity = history_cavity();
    cavity.step(&[1e-4], 0.03, &gate).unwrap();
    let before = states(&cavity);
    let thermal = cavity.network().material_memory().to_vec();
    let radiation = cavity.network().radiation_state().to_vec();
    let frame = cavity.network().last_frame().unwrap().clone();
    let observation = cavity.observe().unwrap();
    assert!(cavity.apply_pressure_impulses(&[], 0.0, &gate).is_err());
    assert!(cavity.apply_modal_impulse(&[], 0.0, &gate).is_err());
    assert!(
        cavity
            .apply_pressure_impulses(&[f64::NAN], 0.0, &gate)
            .is_err()
    );
    assert!(
        cavity
            .apply_modal_impulse(&[0.0], f64::INFINITY, &gate)
            .is_err()
    );
    // A nonzero physical pressure cannot silently lower to zero modal momentum.
    assert!(
        cavity
            .apply_pressure_impulses(&[f64::from_bits(1)], 0.0, &gate)
            .is_err()
    );
    // The surface can stage before the neck violates its original energy ceiling.
    assert!(cavity.apply_modal_impulse(&[1e-7], 1e100, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        cavity.apply_pressure_impulses(&[0.001], 0.002, &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    assert_eq!(states(&cavity), before);
    assert_eq!(cavity.network().material_memory(), thermal);
    assert_eq!(cavity.network().radiation_state(), radiation);
    assert_eq!(cavity.network().last_frame(), Some(&frame));
    assert_eq!(cavity.observe().unwrap(), observation);
    let mut reference = history_cavity();
    reference.step(&[1e-4], 0.03, &gate).unwrap();
    assert_eq!(
        cavity.apply_modal_impulse(&[1e-7], 0.002, &gate).unwrap(),
        reference
            .apply_modal_impulse(&[1e-7], 0.002, &gate)
            .unwrap()
    );
    assert_eq!(
        cavity.step(&[1e-4], 0.03, &gate).unwrap(),
        reference.step(&[1e-4], 0.03, &gate).unwrap()
    );
}
