//! Actual-section diffusion, independent limits, and native pressure/flow work.
use fs_duct::compose::{ComposeOptions, TwoPort, compose_chain};
use fs_duct::fluid::rectangular::{
    MAX_RECTANGULAR_TERMS, RectangularAdmission, RectangularSection,
};
use fs_duct::fluid::{self, FluidDuctError};
use fs_duct::{DuctError, LossModel, Segment};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use std::f64::consts::PI;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn admission(tail: f64) -> RectangularAdmission {
    RectangularAdmission {
        maximum_relative_tail: tail,
        maximum_transverse_phase: 0.8,
    }
}
fn options() -> ComposeOptions {
    ComposeOptions {
        preserve_well_conditioned_bits: false,
        ..ComposeOptions::default()
    }
}
fn near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs(),
        "actual={a:?}, expected={b:?}"
    );
}
// Numerical liquid fixture, not measured-material qualification.
fn liquid(conductivity: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4200.0, conductivity, 0.001).unwrap(),
        4200.0,
    )
    .unwrap()
}

#[test]
fn narrow_rectangle_retains_exact_sidewall_resistance_and_isothermal_storage() {
    let gas = air();
    let (a, b, length, omega) = (0.0005, 0.01, 0.03, 0.001);
    let section = RectangularSection::new(a, b, 256).unwrap();
    let response = section.evaluate(&gas, omega, admission(1e-6)).unwrap();
    let series: f64 = (0..256)
        .map(|i| {
            let m = (2 * i + 1) as f64;
            (m * PI * b / (2.0 * a)).tanh() / m.powi(5)
        })
        .sum();
    let sidewall = 1.0 - 192.0 * a * series / (PI.powi(5) * b);
    let resistance = 12.0 * gas.dynamic_viscosity * length / (b * a.powi(3) * sidewall);
    let reference = gas.characteristic_impedance / section.area_m2();
    let block = response.two_port(length, reference).unwrap();
    let open = compose_chain(&[block], Some(C64::ZERO), options()).unwrap();
    assert!((open.z_in.re / resistance - 1.0).abs() < 2e-7);
    let closed = compose_chain(&[block], None, options()).unwrap();
    let compliance =
        section.area_m2() * length * gas.gamma / (gas.density * gas.sound_speed.powi(2));
    assert!((closed.z_in.im * omega * compliance - 1.0).abs() < 2e-7);
    let circle_resistance = 8.0 * PI * gas.dynamic_viscosity * length / section.area_m2().powi(2);
    assert!((resistance / circle_resistance - 9.860009583521856).abs() < 1e-12);
    let swapped = RectangularSection::new(b, a, 256)
        .unwrap()
        .evaluate(&gas, omega, admission(1e-6))
        .unwrap();
    assert_eq!(response.wave.wavenumber, swapped.wave.wavenumber);
}

// Independent DOUBLE sine sum; production analytically eliminates one axis.
fn double_sum(a: f64, b: f64, omega: f64, diffusivity: f64) -> C64 {
    let mut sum = C64::ZERO;
    for i in 0..512 {
        let m = (2 * i + 1) as f64;
        for j in 0..512 {
            let n = (2 * j + 1) as f64;
            let rate = diffusivity * PI * PI * (m * m / (a * a) + n * n / (b * b));
            sum = sum
                + (C64::new(0.0, omega) / C64::new(rate, omega))
                    .scale(64.0 / (PI.powi(4) * m * m * n * n));
        }
    }
    sum
}

#[test]
fn three_fourier_rungs_converge_inside_the_componentwise_tail_bounds() {
    let gas = air();
    let transport = gas.thermoviscous_state().unwrap().transport();
    let (a, b) = (0.001, 0.0015);
    for omega in [50.0, 500.0, 5000.0] {
        let references = [
            transport.kinematic_viscosity_m2_s(),
            transport.thermal_diffusivity_m2_s(),
        ]
        .map(|d| double_sum(a, b, omega, d));
        let mut previous = [f64::INFINITY; 2];
        for terms in [8, 16, 32] {
            let section = RectangularSection::new(a, b, terms).unwrap();
            assert_eq!(section.terms(), terms);
            let response = section.evaluate(&gas, omega, admission(0.1)).unwrap();
            for (i, actual) in [response.viscous, response.thermal.unwrap()]
                .into_iter()
                .enumerate()
            {
                let error = references[i] - actual.mean;
                assert!(error.re >= -2e-13 && error.re <= actual.tail_real + 2e-13);
                assert!(error.im >= -2e-13 && error.im <= actual.tail_imaginary + 2e-13);
                assert!(
                    error.abs() < 0.2 * previous[i],
                    "terms={terms}, error={error:?}"
                );
                previous[i] = error.abs();
            }
        }
    }
}

#[test]
fn slit_and_thin_layer_limits_preserve_actual_liquid_transport() {
    let material = liquid(0.6);
    for omega in [10.0, 1000.0] {
        let response = RectangularSection::new(0.0001, 0.1, 512)
            .unwrap()
            .evaluate(&material, omega, admission(1e-5))
            .unwrap();
        let z = C64::new(0.0, omega / 1e-6).sqrt().scale(0.00005);
        let denominator = (2.0 * z.re).cosh() + (2.0 * z.im).cos();
        let tanh = C64::new((2.0 * z.re).sinh(), (2.0 * z.im).sin()).scale(1.0 / denominator);
        near(response.viscous.mean, C64::ONE - tanh / z, 8e-4);
        assert!(response.thermal.is_none());
    }
    let section = RectangularSection::new(0.002, 0.004, 4096).unwrap();
    let omega = 100000.0;
    let response = section.evaluate(&material, omega, admission(1e-4)).unwrap();
    let perimeter = 2.0 * (0.002 + 0.004);
    let attenuation =
        perimeter / (2.0 * section.area_m2() * 1500.0) * (1e-6 * omega / 2.0_f64).sqrt();
    assert!((response.wave.wavenumber.im / attenuation - 1.0).abs() < 0.005);
    let different_conductivity = section
        .evaluate(&liquid(60.0), omega, admission(1e-4))
        .unwrap();
    assert_eq!(
        response.wave.wavenumber,
        different_conductivity.wave.wavenumber
    );
    assert_eq!(
        response.wave.characteristic_impedance,
        different_conductivity.wave.characteristic_impedance
    );
}

#[test]
fn native_scattering_is_reciprocal_passive_and_invariant_under_subdivision() {
    let gas = air();
    let section = RectangularSection::new(0.001, 0.004, 256).unwrap();
    let response = section.evaluate(&gas, 3000.0, admission(1e-5)).unwrap();
    let reference = gas.characteristic_impedance / section.area_m2();
    let whole = response.two_port(0.4, reference).unwrap();
    let piece = response.two_port(0.1, reference).unwrap();
    let split = piece
        .star(piece, 1e12)
        .unwrap()
        .0
        .star(piece, 1e12)
        .unwrap()
        .0
        .star(piece, 1e12)
        .unwrap()
        .0;
    for (a, b) in whole.scattering.into_iter().zip(split.scattering) {
        near(a, b, 2e-13);
    }
    let [r, t, reverse, opposite] = whole.scattering;
    assert_eq!(r, opposite);
    assert_eq!(t, reverse);
    for phase in [0.0, 0.7, 1.9, PI] {
        let other = C64::new(phase.cos(), phase.sin()).scale(0.7);
        let outgoing = (r + t * other).norm_sq() + (t + r * other).norm_sq();
        assert!(outgoing < 1.0 + other.norm_sq());
    }
    let matched = compose_chain(
        &[whole],
        Some(response.wave.characteristic_impedance),
        options(),
    )
    .unwrap();
    let phase = response.wave.wavenumber.scale(0.4);
    near(matched.z_in, response.wave.characteristic_impedance, 2e-13);
    near(
        matched.u_out_over_u_in,
        C64::new(phase.re.cos(), phase.re.sin()).scale((-phase.im).exp()),
        2e-13,
    );
    // The new actual-area seam preserves the existing cylindrical operation tree.
    let radius = 0.003;
    let wave = fluid::segment_wave_with_loss(&gas, radius, 3000.0, LossModel::Bessel).unwrap();
    let cylinder = TwoPort::from_wave(
        &Segment::Cylinder {
            radius,
            length: 0.4,
        },
        &wave,
        reference,
    )
    .unwrap();
    let area = TwoPort::from_uniform_wave(&wave, PI * radius * radius, 0.4, reference).unwrap();
    assert_eq!(cylinder.scattering, area.scattering);
}

#[test]
fn unresolved_tail_model_envelope_material_and_cancellation_refuse() {
    for (a, b, n) in [
        (0.0, 1.0, 1),
        (1.0, f64::NAN, 1),
        (1.0, 1.0, 0),
        (1.0, 1.0, MAX_RECTANGULAR_TERMS + 1),
    ] {
        assert!(RectangularSection::new(a, b, n).is_err());
    }
    let mut gas = air();
    let section = RectangularSection::new(0.001, 0.004, 256).unwrap();
    assert!(
        RectangularSection::new(0.001, 0.004, 1)
            .unwrap()
            .evaluate(&gas, 3000.0, admission(1e-6))
            .is_err()
    );
    assert!(section.evaluate(&gas, 100000.0, admission(1e-5)).is_err());
    let mut invalid = admission(1e-5);
    invalid.maximum_transverse_phase = PI;
    assert!(section.evaluate(&gas, 1000.0, invalid).is_err());
    let mut calls = 0;
    assert!(matches!(
        section.evaluate_controlled(&gas, 3000.0, admission(1e-5), || {
            calls += 1;
            calls == 3
        }),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    gas.density = f64::NAN;
    assert!(matches!(
        section.evaluate(&gas, 1000.0, admission(1e-5)),
        Err(FluidDuctError::Material(_))
    ));
}
