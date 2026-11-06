use super::*;
use crate::phase::water::regional::{
    RegionalJunctionControls, RegionalWater, RegionalWaterControls,
};
use crate::visco::GeneralizedMaxwell;
use fs_math::det;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;

fn chart(pressure_pa: f64) -> RegionalWater {
    let (lo, hi) = if pressure_pa < 1e6 {
        (300.0, 500.0)
    } else {
        (625.0, 645.0)
    };
    RegionalWater::new(
        Pressure::new(pressure_pa),
        RegionalWaterControls {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
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
    .unwrap()
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "{actual:e} != {expected:e}"
    );
}

fn complex_near(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn g1_original_regional_density_limits_and_passive_spectrum() {
    for pressure in [1e5, 20e6] {
        let chart = chart(pressure);
        for quality in [0.01, 0.2, 0.8] {
            let source = chart.state_at_quality(quality, || false).unwrap();
            let law = RelaxingRegionalWater::try_new(source, 0.003).unwrap();
            assert_eq!(law.source(), &source);
            assert_eq!(law.equilibrium().source(), &source);
            assert_eq!(law.frozen().source(), &source);
            let equilibrium = law.equilibrium().state();
            let frozen = law.frozen().state();
            let ke = equilibrium.adiabatic_bulk_modulus_pa();
            let kf = frozen.adiabatic_bulk_modulus_pa();
            let rho = source.bulk_density_kg_m3();
            assert_eq!(law.arm().delta_bulk_modulus_pa, kf - ke);
            assert_eq!(law.arm().relaxation_time_s, 0.003);
            let dc = law.response(0.0).unwrap();
            assert_eq!(dc.bulk_modulus_pa(), C64::from_re(ke));
            assert_eq!(dc.wavenumber_rad_m(), C64::ZERO);
            near(
                dc.specific_impedance_pa_s_m().re,
                equilibrium.characteristic_impedance_pa_s_m(),
                3e-15,
            );
            let high = law.response(1e8 / law.arm().relaxation_time_s).unwrap();
            near(high.bulk_modulus_pa().re, kf, 3e-15);
            near(
                (1e8 / law.arm().relaxation_time_s) / high.wavenumber_rad_m().re,
                frozen.sound_speed_m_s(),
                3e-15,
            );
            for scaled_frequency in [1e-8, 1e-4, 0.1, 1.0, 10.0, 1e4, 1e8] {
                let omega = scaled_frequency / law.arm().relaxation_time_s;
                let response = law.response(omega).unwrap();
                let bulk = response.bulk_modulus_pa();
                let compliance = response.compressibility_per_pa();
                let susceptibility = response.density_susceptibility_kg_m3_pa();
                let wave = response.wavenumber_rad_m();
                let impedance = response.specific_impedance_pa_s_m();
                assert_eq!(response.density_kg_m3().to_bits(), rho.to_bits());
                assert!(bulk.re >= ke && bulk.re <= kf * (1.0 + 3e-15));
                assert!(bulk.im < 0.0 && compliance.im > 0.0 && wave.im > 0.0);
                assert!(susceptibility.im > 0.0 && impedance.im < 0.0);
                complex_near(bulk * compliance, C64::ONE, 5e-15);
                complex_near(bulk * susceptibility, C64::from_re(rho), 5e-15);
                complex_near(wave * impedance, C64::from_re(rho * omega), 5e-15);
                complex_near(wave * wave * bulk, C64::from_re(rho * omega * omega), 8e-15);

                // Unit peak pressure: decay of real acoustic energy flux equals
                // the local work dissipated by the complex compression law.
                let flux_loss = wave.im * impedance.recip().re;
                let heat_rate = -0.5 * omega * bulk.im * compliance.norm_sq();
                near(flux_loss, heat_rate, 2e-13);
            }
        }
    }
}

fn harmonic_history(law: &RelaxingRegionalWater, omega: f64, steps: usize) -> (C64, f64) {
    // An independent time-domain owner integrates exact linear strain ramps.
    // Discard transients, then measure pressure quadratures and actual dashpot heat.
    let arm = law.arm();
    let model = GeneralizedMaxwell::new(
        law.equilibrium().state().adiabatic_bulk_modulus_pa(),
        vec![(arm.delta_bulk_modulus_pa, arm.relaxation_time_s)],
    )
    .unwrap();
    let mut state = model.state();
    let amplitude = 1e-6;
    state.strain = amplitude;
    let dt = core::f64::consts::TAU / (omega * steps as f64);
    let mut quadratures = C64::ZERO;
    let mut heat = 0.0;
    for cycle in 0..12 {
        let heat_before = state.dissipated;
        for step in 1..=steps {
            let angle = core::f64::consts::TAU * (step % steps) as f64 / steps as f64;
            let (sin, cos) = det::sin_cos(angle);
            let pressure = model.step(&mut state, amplitude * cos, dt);
            if cycle == 11 {
                quadratures = quadratures + C64::new(pressure * cos, pressure * sin);
            }
        }
        if cycle == 11 {
            heat = state.dissipated - heat_before;
        }
    }
    (
        quadratures.scale(2.0 / (steps as f64 * amplitude)),
        heat / (amplitude * amplitude),
    )
}

#[test]
fn g1_frequency_response_matches_native_ramp_history_and_cycle_heat_with_refinement() {
    let source = chart(20e6).state_at_quality(0.2, || false).unwrap();
    let law = RelaxingRegionalWater::try_new(source, 0.003).unwrap();
    for omega_tau in [0.2, 2.0] {
        let omega = omega_tau / law.arm().relaxation_time_s;
        let exact = law.response(omega).unwrap().bulk_modulus_pa();
        let (coarse, _) = harmonic_history(&law, omega, 512);
        let (fine, heat) = harmonic_history(&law, omega, 1024);
        assert!((fine - exact).abs() < 0.3 * (coarse - exact).abs());
        complex_near(fine, exact, 1e-5);
        near(heat, -core::f64::consts::PI * exact.im, 1e-5);
        assert!(heat > 0.0);
    }
}

#[test]
fn g3_declared_compression_time_rescales_frequency_without_changing_density() {
    let source = chart(20e6).state_at_quality(0.6, || false).unwrap();
    let a = RelaxingRegionalWater::try_new(source, 0.003).unwrap();
    let b = RelaxingRegionalWater::try_new(source, 0.03).unwrap();
    for omega in [1.0, 70.0, 7000.0] {
        let ra = a.response(omega).unwrap();
        let rb = b.response(omega / 10.0).unwrap();
        complex_near(ra.bulk_modulus_pa(), rb.bulk_modulus_pa(), 3e-15);
        complex_near(
            ra.density_susceptibility_kg_m3_pa(),
            rb.density_susceptibility_kg_m3_pa(),
            3e-15,
        );
        complex_near(
            ra.wavenumber_rad_m(),
            rb.wavenumber_rad_m().scale(10.0),
            3e-15,
        );
        assert_eq!(ra.density_kg_m3(), rb.density_kg_m3());
    }
}

#[test]
fn g0_missing_coexistence_and_invalid_rates_or_frequencies_refuse() {
    let chart = chart(1e5);
    let source = chart.state_at_quality(0.2, || false).unwrap();
    for tau in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(RelaxingRegionalWater::try_new(source, tau).is_err());
    }
    for source in [
        chart.state_at_quality(0.0, || false).unwrap(),
        chart.state_at_quality(1.0, || false).unwrap(),
        chart
            .state_at_temperature(Temperature::new(320.0), || false)
            .unwrap(),
    ] {
        assert!(RelaxingRegionalWater::try_new(source, 0.003).is_err());
    }
    let law = RelaxingRegionalWater::try_new(source, 0.003).unwrap();
    for omega in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(law.response(omega).is_err());
    }
    let arm = law.arm();
    assert_eq!(arm.response(0.0).unwrap(), C64::ZERO);
    complex_near(
        arm.response(arm.relaxation_time_s.recip()).unwrap(),
        C64::new(0.5, -0.5).scale(arm.delta_bulk_modulus_pa),
        3e-15,
    );
    for invalid_value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            BulkRelaxationArm {
                delta_bulk_modulus_pa: invalid_value,
                ..arm
            }
            .response(1.0)
            .is_err()
        );
        assert!(
            BulkRelaxationArm {
                relaxation_time_s: invalid_value,
                ..arm
            }
            .response(1.0)
            .is_err()
        );
    }
    let slow = BulkRelaxationArm {
        delta_bulk_modulus_pa: 1.0,
        relaxation_time_s: f64::MAX,
    };
    assert_eq!(slow.response(f64::MAX).unwrap(), C64::ONE);
}
