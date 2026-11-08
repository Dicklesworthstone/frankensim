//! Joint heat/molecular diffusion in actual rectangular passages.

use fs_duct::DuctError;
use fs_duct::fluid::FluidDuctError;
use fs_duct::fluid::rectangular::{RectangularAdmission, RectangularSection, RectangularWave};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;
use std::f64::consts::PI;

// Manufactured reservoir data test the physical equations, not measured air
// relaxation. The original pressure-dependent material rate remains in use.
fn material(count: usize, time_scale: f64, pressure: f64) -> RelaxationState {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    let species = SpeciesId::new("synthetic").unwrap();
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: (0..count)
            .map(|i| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * i as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * (0.002 + 0.01 * i as f64) * time_scale,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap();
    law.evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
        .unwrap()
}

fn admission(tolerance: f64) -> RectangularAdmission {
    RectangularAdmission {
        maximum_relative_tail: tolerance,
        maximum_transverse_phase: 0.8,
    }
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:?} != {expected:?}"
    );
}
fn compliance(response: RectangularWave, omega: f64, area: f64) -> C64 {
    response.wave.wavenumber / response.wave.characteristic_impedance / C64::from_re(omega * area)
}

// Independent +i omega reservoir sum, not caloric_response's operation tree.
fn heat_capacity(state: &RelaxationState, omega: f64) -> C64 {
    state.processes().iter().fold(
        C64::from_re(state.frozen_cv_j_kg_k() + state.gas().specific_gas_constant),
        |cp, process| {
            cp + C64::from_re(process.heat_capacity_j_kg_k)
                / C64::new(1.0, omega * process.relaxation_time_s)
        },
    )
}

// Independently solve each TWO-dimensional thermal eigenmode's energy
// equation: (s cp(s) + conductivity*lambda/rho) T_mn = s p_mn/rho.
// Production analytically eliminates one axis and uses a one-dimensional sum.
fn temperature_response(a: f64, b: f64, omega: f64, state: &RelaxationState) -> C64 {
    let gas = state.gas();
    let conductivity = gas
        .thermoviscous_state()
        .unwrap()
        .transport()
        .thermal_conductivity_w_m_k();
    let s = C64::new(0.0, omega);
    let cp = heat_capacity(state, omega);
    let mut sum = C64::ZERO;
    for i in 0..512 {
        let m = (2 * i + 1) as f64;
        for j in 0..512 {
            let n = (2 * j + 1) as f64;
            let rate = (conductivity / gas.density) * PI * PI * (m * m / (a * a) + n * n / (b * b));
            sum = sum
                + (s / (s * cp + C64::from_re(rate))).scale(64.0 / (PI.powi(4) * m * m * n * n));
        }
    }
    sum
}

#[test]
fn absent_reservoirs_retain_the_original_rectangular_response_exactly() {
    let state = material(0, 1.0, 101325.0);
    let section = RectangularSection::new(0.0004, 0.002, 256).unwrap();
    for omega in [10.0, 500.0, 5000.0] {
        let original = section
            .evaluate(state.gas(), omega, admission(1e-6))
            .unwrap();
        let molecular = section
            .evaluate_molecular(&state, omega, admission(1e-6))
            .unwrap();
        assert_eq!(molecular.wave.wavenumber, original.wave.wavenumber);
        assert_eq!(
            molecular.wave.characteristic_impedance,
            original.wave.characteristic_impedance
        );
        assert_eq!(
            molecular.thermal.unwrap().mean,
            original.thermal.unwrap().mean
        );
        assert_eq!(
            molecular.thermal.unwrap().tail_real,
            original.thermal.unwrap().tail_real
        );
        let reference = state.gas().characteristic_impedance / section.area_m2();
        assert_eq!(
            molecular.two_port(0.1, reference).unwrap().scattering,
            original.two_port(0.1, reference).unwrap().scattering
        );
    }
}

#[test]
fn coupled_energy_equation_converges_with_the_first_quadrant_tail_bound() {
    let state = material(2, 1.0, 101325.0);
    let gas = state.gas();
    let (a, b) = (0.001, 0.0015);
    for omega in [50.0, 500.0, 5000.0] {
        let cp = heat_capacity(&state, omega);
        let temperature = temperature_response(a, b, omega, &state);
        let expected_h = temperature * cp;
        let expected_kappa = (C64::ONE - temperature.scale(gas.specific_gas_constant))
            .scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature))
            .conj();
        let mut previous = f64::INFINITY;
        for terms in [8, 16, 32] {
            let section = RectangularSection::new(a, b, terms).unwrap();
            let wave = section
                .evaluate_molecular(&state, omega, admission(0.1))
                .unwrap();
            let thermal = wave.thermal.unwrap();
            let error = expected_h - thermal.mean;
            assert!(error.re >= -3e-12 && error.re <= thermal.tail_real + 3e-12);
            assert!(error.im >= -3e-12 && error.im <= thermal.tail_imaginary + 3e-12);
            assert!(error.abs() < 0.2 * previous);
            previous = error.abs();
            let conductivity = gas
                .thermoviscous_state()
                .unwrap()
                .transport()
                .thermal_conductivity_w_m_k();
            let imaginary_ratio = omega * gas.density * cp.re / conductivity / (PI * PI / (a * a));
            let old_pure_imaginary_bound =
                imaginary_ratio.powi(2) / (40.0 * PI * PI * (terms as f64).powi(5));
            if terms == 32 {
                assert!(
                    error.re > old_pure_imaginary_bound,
                    "the positive real diffusion argument needs its N^-3 tail: omega={omega}, error={}, old bound={old_pure_imaginary_bound}",
                    error.re,
                );
            }
        }
        let fine = RectangularSection::new(a, b, 256).unwrap();
        near(
            compliance(
                fine.evaluate_molecular(&state, omega, admission(1e-6))
                    .unwrap(),
                omega,
                fine.area_m2(),
            ),
            expected_kappa,
            2e-8,
        );
    }
}

#[test]
fn isothermal_storage_and_equilibrated_frozen_heat_capacities_have_correct_limits() {
    let section = RectangularSection::new(0.0005, 0.01, 256).unwrap();
    let omega = 0.001;
    let fast = material(2, 1e-9, 101325.0);
    let slow = material(2, 1e9, 101325.0);
    let fast_compliance = compliance(
        section
            .evaluate_molecular(&fast, omega, admission(1e-6))
            .unwrap(),
        omega,
        section.area_m2(),
    );
    let slow_compliance = compliance(
        section
            .evaluate_molecular(&slow, omega, admission(1e-6))
            .unwrap(),
        omega,
        section.area_m2(),
    );
    near(
        fast_compliance,
        C64::from_re(1.0 / fast.gas().pressure),
        2e-6,
    );
    near(slow_compliance, fast_compliance, 1e-10);
    // Compare finite-frequency endpoint physics, keeping the original rho,
    // viscosity and conductivity and the independently known frozen Cv/Cp.
    let section = RectangularSection::new(0.001, 0.004, 256).unwrap();
    let gas = slow.gas();
    let transport = gas.thermoviscous_state().unwrap().transport();
    let cv = slow.frozen_cv_j_kg_k();
    let frozen = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(gas.density, slow.frozen_sound_speed_m_s()).unwrap(),
        FluidTransportState::try_new(
            gas.density,
            cv + gas.specific_gas_constant,
            transport.thermal_conductivity_w_m_k(),
            transport.dynamic_viscosity_pa_s(),
        )
        .unwrap(),
        cv,
    )
    .unwrap();
    for omega in [100.0, 3000.0] {
        let equilibrated = section
            .evaluate_molecular(&fast, omega, admission(1e-6))
            .unwrap();
        let expected_eq = section
            .evaluate(fast.gas(), omega, admission(1e-6))
            .unwrap();
        near(
            equilibrated.wave.wavenumber,
            expected_eq.wave.wavenumber,
            2e-8,
        );
        let unrelaxed = section
            .evaluate_molecular(&slow, omega, admission(1e-6))
            .unwrap();
        let expected_frozen = section.evaluate(&frozen, omega, admission(1e-6)).unwrap();
        near(
            unrelaxed.wave.wavenumber,
            expected_frozen.wave.wavenumber,
            2e-8,
        );
        near(
            unrelaxed.wave.characteristic_impedance,
            expected_frozen.wave.characteristic_impedance,
            2e-8,
        );
        assert!(
            (equilibrated.wave.wavenumber - unrelaxed.wave.wavenumber).abs()
                > 1e-3 * unrelaxed.wave.wavenumber.abs()
        );
    }
}

#[test]
fn pressure_dependent_rates_preserve_scaled_physics_and_passive_scattering() {
    let mut reference = None;
    for scale in [1.0_f64, 0.5, 2.0] {
        let state = material(2, 1.0, 101325.0 * scale);
        let section = RectangularSection::new(0.001 / scale, 0.002 / scale, 256).unwrap();
        let omega = 500.0 * scale;
        let wave = section
            .evaluate_molecular(&state, omega, admission(1e-6))
            .unwrap();
        let impedance = state.gas().characteristic_impedance / section.area_m2();
        let block = wave.two_port(0.2 / scale, impedance).unwrap();
        let (k, z) = (
            wave.wave.wavenumber.scale(1.0 / scale),
            wave.wave
                .characteristic_impedance
                .scale(1.0 / scale.powi(3)),
        );
        if let Some((old_k, old_z, old_s)) = reference {
            near(k, old_k, 3e-13);
            near(z, old_z, 3e-13);
            for (actual, expected) in block.scattering.into_iter().zip(old_s) {
                near(actual, expected, 1e-12);
            }
        } else {
            reference = Some((k, z, block.scattering));
        }
        let [r, t, reverse, opposite] = block.scattering;
        assert_eq!(r, opposite);
        assert_eq!(t, reverse);
        for phase in [0.0_f64, 0.7, 1.9, PI] {
            let other = C64::new(phase.cos(), phase.sin()).scale(0.7);
            let outgoing = (r + t * other).norm_sq() + (t + r * other).norm_sq();
            assert!(outgoing < 1.0 + other.norm_sq());
        }
        let replay = section
            .evaluate_molecular(&state, omega, admission(1e-6))
            .unwrap();
        assert_eq!(replay.wave.wavenumber, wave.wave.wavenumber);
    }
}

#[test]
fn unresolved_joint_loss_invalid_envelopes_and_cancellation_refuse() {
    let state = material(2, 1.0, 101325.0);
    let section = RectangularSection::new(0.001, 0.004, 256).unwrap();
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY, 100000.0] {
        assert!(
            section
                .evaluate_molecular(&state, omega, admission(1e-6))
                .is_err()
        );
    }
    for tolerance in [0.0, -1.0, 1.0, f64::NAN] {
        assert!(
            section
                .evaluate_molecular(&state, 500.0, admission(tolerance))
                .is_err()
        );
    }
    assert!(
        RectangularSection::new(0.001, 0.004, 2)
            .unwrap()
            .evaluate_molecular(&state, 500.0, admission(1e-6))
            .is_err()
    );
    let mut polls = 0;
    let result = section.evaluate_molecular_controlled(&state, 500.0, admission(1e-6), || {
        polls += 1;
        polls == 12 // after viscous work, during the actual thermal solve
    });
    assert!(matches!(
        result,
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert!(
        section
            .evaluate_molecular(&state, 500.0, admission(1e-6))
            .is_ok()
    );
}
