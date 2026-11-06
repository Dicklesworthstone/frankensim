use super::*;
use crate::phase::water::regional::{
    RegionalJunctionControls, RegionalWater, RegionalWaterControls,
};
use fs_math::{c64::C64, det};
use fs_qty::Pressure;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;

fn material(pressure: f64, tau: f64) -> RelaxingRegionalWater {
    let bounds = if pressure < 1e6 {
        [300.0, 500.0]
    } else {
        [625.0, 645.0]
    };
    let law = RegionalWater::new(
        Pressure::new(pressure),
        RegionalWaterControls {
            minimum_temperature_k: bounds[0],
            maximum_temperature_k: bounds[1],
            enthalpy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            maximum_saturation_entropy_defect_j_kg_k: 0.1,
            pressure: DensePressureSolve {
                pressure_tolerance_pa: 1e-4,
                density_tolerance_kg_m3: 1e-8,
                maximum_iterations: 128,
            },
            junctions: RegionalJunctionControls {
                maximum_enthalpy_jump_j_kg: 200.0,
                maximum_entropy_jump_j_kg_k: 0.25,
                maximum_volume_jump_m3_kg: 1e-5,
                maximum_temperature_gap_k: 1e-8,
                maximum_iterations: 96,
            },
        },
        || false,
    )
    .unwrap();
    RelaxingRegionalWater::try_new(law.state_at_quality(0.3, || false).unwrap(), tau).unwrap()
}

fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-300),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_frozen_jump_relaxes_with_independent_heat_and_fixed_mass_geometry() {
    for pressure in [1e5, 20e6] {
        let model = material(pressure, 0.01);
        let mut history = WaterCompressionHistory::try_new(model, 1e-3).unwrap();
        let rho0 = model.source().bulk_density_kg_m3();
        let v0 = model.source().specific_volume_m3_kg();
        assert_eq!(history.state().density_kg_m3.to_bits(), rho0.to_bits());
        let e = 1e-5;
        let jump = history.advance(e, 0.0).unwrap();
        let kf = model.frozen().state().adiabatic_bulk_modulus_pa();
        let ke = model.equilibrium().state().adiabatic_bulk_modulus_pa();
        let dk = model.arm().delta_bulk_modulus_pa;
        near(jump.state.pressure_perturbation_pa, kf * e, 4e-15);
        assert_eq!(jump.dissipated_energy_j_m3, 0.0);
        near(jump.compression_work_j_m3, 0.5 * kf * e * e, 4e-15);
        near(
            jump.state.density_kg_m3 * jump.state.specific_volume_m3_kg,
            1.0,
            3e-15,
        );
        near(jump.state.specific_volume_m3_kg / v0, 1.0 - e, 3e-15);
        near(
            jump.state.density_perturbation_kg_m3,
            rho0 * e / (1.0 - e),
            3e-15,
        );
        for ratio in [1e-10, 0.1, 1.0, 1000.0] {
            let mut held = history;
            let step = held.advance(e, ratio * 0.01).unwrap();
            let decay = det::exp(-ratio);
            near(
                step.state.pressure_perturbation_pa,
                (ke + dk * decay) * e,
                8e-15,
            );
            near(
                step.dissipated_energy_j_m3,
                0.5 * dk * e * e * (-det::expm1(-2.0 * ratio)),
                8e-15,
            );
            assert_eq!(step.compression_work_j_m3, 0.0);
            assert_eq!(step.state.density_kg_m3, jump.state.density_kg_m3);
            assert_eq!(*held.source().source(), *model.source());
            assert!(step.energy_residual_j_m3.abs() <= 1e-14 * jump.state.stored_energy_j_m3);
        }
    }
}

#[test]
fn g1_exact_ramp_matches_independent_pressure_work_heat_and_subdivision() {
    let model = material(20e6, 0.01);
    let mut whole = WaterCompressionHistory::try_new(model, 0.001).unwrap();
    whole.advance(0.0004, 0.0).unwrap();
    let mut split = whole;
    let initial = whole.state();
    let target = -0.0002;
    let duration = 0.007;
    let step = whole.advance(target, duration).unwrap();
    let ke = model.equilibrium().state().adiabatic_bulk_modulus_pa();
    let dk = model.arm().delta_bulk_modulus_pa;
    let tau = model.arm().relaxation_time_s;
    let rate = (target - initial.compression) / duration;
    // Independent Simpson integration of the original differential solution.
    let mut work = 0.0;
    let mut heat = 0.0;
    let n = 4096;
    for i in 0..=n {
        let t = duration * i as f64 / n as f64;
        let q =
            initial.branch_compression * det::exp(-t / tau) + rate * tau * (-det::expm1(-t / tau));
        let e = initial.compression + rate * t;
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        work += weight * (ke * e + dk * q) * rate;
        heat += weight * dk * q * q / tau;
    }
    work *= duration / (3.0 * n as f64);
    heat *= duration / (3.0 * n as f64);
    near(step.compression_work_j_m3, work, 3e-12);
    near(step.dissipated_energy_j_m3, heat, 3e-12);
    let prepared = split.prepare_ramp(duration / 32.0).unwrap();
    let mut split_work = 0.0;
    let mut split_heat = 0.0;
    for i in 1..=32 {
        let target_i = initial.compression + (target - initial.compression) * i as f64 / 32.0;
        let accepted = split.advance_prepared(target_i, &prepared).unwrap();
        split_work += accepted.compression_work_j_m3;
        split_heat += accepted.dissipated_energy_j_m3;
    }
    near(
        split.state().pressure_perturbation_pa,
        step.state.pressure_perturbation_pa,
        3e-13,
    );
    near(
        split.state().branch_compression,
        step.state.branch_compression,
        3e-13,
    );
    near(split_work, step.compression_work_j_m3, 3e-13);
    near(split_heat, step.dissipated_energy_j_m3, 3e-13);
}

#[test]
fn g1_periodic_compression_loss_matches_the_same_harmonic_material() {
    let tau = 0.01;
    let model = material(1e5, tau);
    let mut history = WaterCompressionHistory::try_new(model, 1e-3).unwrap();
    let amplitude = 1e-5;
    let points = 1024;
    let period = core::f64::consts::TAU * tau;
    let prepared = history.prepare_ramp(period / points as f64).unwrap();
    let mut last_heat = 0.0;
    let mut last_work = 0.0;
    for cycle in 0..8 {
        last_heat = 0.0;
        last_work = 0.0;
        for i in 1..=points {
            let phase = core::f64::consts::TAU * i as f64 / points as f64;
            let target = amplitude * det::sin_cos(phase).0;
            let step = history.advance_prepared(target, &prepared).unwrap();
            last_heat += step.dissipated_energy_j_m3;
            last_work += step.compression_work_j_m3;
        }
        assert!(last_heat > 0.0, "cycle {cycle}");
    }
    let modulus: C64 = model.response(tau.recip()).unwrap().bulk_modulus_pa();
    let expected = -core::f64::consts::PI * modulus.im * amplitude * amplitude;
    near(last_heat, expected, 8e-6);
    near(last_work, last_heat, 2e-12);
}

#[test]
fn g0_invalid_ramps_are_atomic_and_tiny_density_changes_remain_visible() {
    let model = material(1e5, 0.01);
    for bound in [0.0, -0.01, 1.0, f64::NAN, f64::INFINITY] {
        assert!(WaterCompressionHistory::try_new(model, bound).is_err());
    }
    let mut history = WaterCompressionHistory::try_new(model, 0.001).unwrap();
    history.advance(0.0001, 0.002).unwrap();
    for (compression, dt) in [
        (0.0011, 0.002),
        (-0.0011, 0.002),
        (f64::NAN, 0.002),
        (f64::INFINITY, 0.002),
        (0.0, -0.1),
        (0.0, f64::NAN),
        (0.0, f64::INFINITY),
    ] {
        let before = history;
        assert!(history.advance(compression, dt).is_err());
        assert_eq!(history, before);
    }
    let other = WaterCompressionHistory::try_new(material(1e5, 0.02), 0.001).unwrap();
    let wrong_ramp = other.prepare_ramp(0.01).unwrap();
    let before = history;
    assert!(history.advance_prepared(0.0, &wrong_ramp).is_err());
    assert_eq!(history, before);
    // A finite but overflowing dt/tau still has representable small memory.
    // Refuse it instead of silently replacing the ramp by an infinite hold.
    let fast = RelaxingRegionalWater::try_new(*model.source(), 1e-210).unwrap();
    let mut fast_history = WaterCompressionHistory::try_new(fast, 0.001).unwrap();
    let before = fast_history;
    assert!(fast_history.advance(0.0001, 1e100).is_err());
    assert_eq!(fast_history, before);
    let mut tiny = WaterCompressionHistory::try_new(model, 0.001).unwrap();
    let sample = tiny.advance(1e-20, 0.0).unwrap().state;
    assert_eq!(sample.density_kg_m3, model.source().bulk_density_kg_m3());
    assert!(sample.density_perturbation_kg_m3 > 0.0);
    near(
        sample.density_perturbation_kg_m3,
        model.source().bulk_density_kg_m3() * 1e-20,
        3e-15,
    );
}
