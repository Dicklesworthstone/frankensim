use super::*;
use crate::modal_acoustic_time::ModalAcousticState;
use crate::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingConfig};
use fs_exec::CancelGate;
use fs_material::gas::GasSpec;

fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}
fn geometry() -> HelmholtzCavity {
    HelmholtzCavity { volume_m3: 0.01, neck_radius_m: 0.01, neck_length_m: 0.02 }
}
fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8, max_connections: 4, max_setup_terms: 1024,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1.0e6, maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1.0e-10, energy_absolute_tolerance_j: 1.0e-12,
        energy_relative_tolerance: 1.0e-8,
    }
}
fn body() -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: 700.0, damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::ZERO,
    }], ModalAcousticTimeBudget::audible_reference()).unwrap()
}

#[test]
fn lowering_preserves_original_helmholtz_storage_and_selected_radiation_loss() {
    let geo = geometry();
    for pressure in [50_000.0, 101_325.0, 200_000.0] {
        let gas = gas(pressure);
        let model = HelmholtzModalCavity::new(geo, gas,
            CavityNeckLoss::RadiationAtResonance(fs_phs::MouthFlange::Flanged)).unwrap();
        let original = fs_phs::helmholtz_resonator_radiating(geo.volume_m3, geo.neck_radius_m,
            geo.neck_length_m, gas.density, gas.sound_speed, fs_phs::MouthFlange::Flanged).unwrap();
        let area = model.neck_area_m2();
        for (x, velocity) in [(0.001, 0.1), (-0.002, 0.0), (0.0, -0.02)] {
            let energy = 0.5 * model.mass * velocity * velocity + 0.5 * model.stiffness * x * x;
            let source = original.hamiltonian(&[area * x, model.mass * velocity / area]);
            assert!((source - energy).abs() <= 2.0e-14 * energy);
            let source_rate = original.structure().1[3] * (area * velocity).powi(2);
            assert!((source_rate - model.resistance * velocity * velocity).abs()
                <= 2.0e-14 * source_rate);
        }
    }
}

#[test]
fn pressure_scaling_changes_inertia_stiffness_and_loss_without_fake_retuning() {
    let loss = CavityNeckLoss::RadiationAtResonance(fs_phs::MouthFlange::Flanged);
    let a = HelmholtzModalCavity::new(geometry(), gas(50_000.0), loss).unwrap();
    let b = HelmholtzModalCavity::new(geometry(), gas(100_000.0), loss).unwrap();
    for ratio in [b.mass / a.mass, b.stiffness / a.stiffness, b.resistance / a.resistance] {
        assert!((ratio - 2.0).abs() < 1.0e-12);
    }
    // At fixed temperature this model's density scaling does not change the
    // ISOLATED resonance, although its loading of a finite body does change.
    assert!((b.isolated_omega_rad_s() / a.isolated_omega_rad_s() - 1.0).abs() < 1.0e-12);
}

#[test]
fn collective_volume_port_preserves_cancellation_and_the_original_energy() {
    let cavity = HelmholtzModalCavity::new(geometry(), gas(101_325.0), CavityNeckLoss::Resistance(0.0)).unwrap();
    let mut plate = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: 700.0, damping_ratio: 0.0, pressure_per_modal_velocity: C64::ZERO,
    }; 2], ModalAcousticTimeBudget::audible_reference()).unwrap();
    plate.restore_states(&[ModalAcousticState {
        displacement_m_sqrt_kg: 0.001, velocity_m_sqrt_kg_per_s: 0.0,
    }; 2]).unwrap();
    let neck = cavity.neck_component(48_000, ModalAcousticTimeBudget::audible_reference()).unwrap();
    let links = cavity.connections(0, 1, &[0.01, -0.01]).unwrap();
    assert_eq!(links.len(), 1); // one common pressure, not per-mode springs
    let gate = CancelGate::new_clock_free();
    let mut network = CoupledModalSystem::new(vec![plate, neck], links, config(), &gate).unwrap();
    for _ in 0..64 {
        network.step(&[0.0; 3]).unwrap();
        assert_eq!(network.components()[1].states()[0], ModalAcousticState::default());
        assert_eq!(network.components()[0].states()[0], network.components()[0].states()[1]);
    }
}

#[test]
fn cavity_and_body_exchange_energy_both_ways_in_one_accepted_step() {
    let cavity = HelmholtzModalCavity::new(geometry(), gas(101_325.0), CavityNeckLoss::Resistance(2000.0)).unwrap();
    let gate = CancelGate::new_clock_free();
    for excite_neck in [false, true] {
        let mut plate = body();
        let mut neck = cavity.neck_component(48_000, ModalAcousticTimeBudget::audible_reference()).unwrap();
        let driven = if excite_neck { &mut neck } else { &mut plate };
        driven.restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.0, velocity_m_sqrt_kg_per_s: 0.01,
        }]).unwrap();
        let mut network = CoupledModalSystem::new(vec![plate, neck],
            cavity.connections(0, 1, &[0.02]).unwrap(), config(), &gate).unwrap();
        let initial = network.total_energy_j().unwrap();
        let mut lost = 0.0;
        for _ in 0..256 {
            let frame = network.step(&[0.0; 2]).unwrap();
            lost += frame.connection_dissipation_j + frame.component_dissipation_j;
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        let receiving = if excite_neck { 0 } else { 1 };
        assert!(network.components()[receiving].states()[0].velocity_m_sqrt_kg_per_s.abs() > 1.0e-8);
        assert!((network.total_energy_j().unwrap() + lost - initial).abs() < 1.0e-12);
        assert!(lost > 0.0);
    }
}

#[test]
fn invalid_geometry_gas_loss_and_attachment_mapping_refuse() {
    let mut bad = geometry(); bad.volume_m3 = 0.0;
    assert!(HelmholtzModalCavity::new(bad, gas(101_325.0), CavityNeckLoss::Resistance(0.0)).is_err());
    let mut stale = gas(101_325.0); stale.sound_speed *= 2.0;
    assert!(HelmholtzModalCavity::new(geometry(), stale, CavityNeckLoss::Resistance(0.0)).is_err());
    for r in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(HelmholtzModalCavity::new(geometry(), gas(101_325.0), CavityNeckLoss::Resistance(r)).is_err());
    }
    let model = HelmholtzModalCavity::new(geometry(), gas(101_325.0), CavityNeckLoss::Resistance(0.0)).unwrap();
    assert!(model.connections(0, 0, &[1.0]).is_err());
    assert!(model.connections(0, 1, &[]).is_err());
    assert!(model.connections(0, 1, &[f64::NAN]).is_err());
}
