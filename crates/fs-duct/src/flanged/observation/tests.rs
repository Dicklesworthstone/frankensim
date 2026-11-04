//! Physical field controls, independent aperture quadrature and source work.

use super::*;
use crate::flanged::flanged_modal_impedance_from_material;
use crate::levine_schwinger::gauss_legendre;
use crate::modal::{AxisymmetricMouthBasis, ModeCount, PreparedModalDuct};
use crate::network::PortExcitation;
use crate::{Duct, LossModel, Segment};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;
use std::cell::Cell;

const PI: f64 = core::f64::consts::PI;

fn close(actual: C64, expected: C64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * scale.max(expected.abs()).max(1e-30),
        "{actual:?} != {expected:?}"
    );
}

fn liquid() -> FluidAcousticState {
    // Numerical constitutive control, not a measured liquid material card.
    FluidAcousticState::try_new(1000.0, 1500.0).unwrap()
}

fn integrated_power(observer: &FlangedModalObserver, flow: &[C64]) -> f64 {
    2.0 * PI
        * gauss_legendre(64)
            .iter()
            .map(|&(mu, weight)| {
                weight * observer.far_field(flow, mu).unwrap().power_per_steradian_w
            })
            .sum::<f64>()
}

fn matrix_power(matrix: &[C64], flow: &[C64]) -> f64 {
    let n = flow.len();
    (0..n)
        .map(|row| {
            let pressure = (0..n).fold(C64::ZERO, |p, col| p + matrix[row * n + col] * flow[col]);
            0.5 * (pressure * flow[row].conj()).re
        })
        .sum()
}

#[test]
fn piston_normalization_and_outgoing_pressure_phase_are_physical() {
    let exterior = liquid();
    let radius = 0.025;
    let omega = 3.0 * exterior.sound_speed_m_s() / radius;
    let observer = FlangedModalObserver::from_material(&exterior, radius, omega, 1).unwrap();
    let flow = [C64::new(1e-6, -0.3e-6)];
    let monopole = C64::new(0.0, -exterior.density_kg_m3() * omega / (2.0 * PI)) * flow[0];
    for mu in [0.0, 0.3, 0.8, 1.0] {
        let x = 3.0 * (1.0_f64 - mu * mu).sqrt();
        let directivity = if x == 0.0 {
            1.0
        } else {
            2.0 * fs_math::special::bessel_j1(x) / x
        };
        let field = observer.far_field(&flow, mu).unwrap();
        close(
            field.amplitude_pa_m,
            monopole.scale(directivity),
            monopole.abs(),
        );
        close(
            C64::from_re(field.power_per_steradian_w),
            C64::from_re(
                field.amplitude_pa_m.norm_sq() / (2.0 * exterior.characteristic_impedance_pa_s_m()),
            ),
            field.power_per_steradian_w,
        );
    }
    let distance = 100.0 * radius;
    let quarter_wavelength = 0.5 * PI * exterior.sound_speed_m_s() / omega;
    let first = observer.far_field_pressure(&flow, 0.8, distance).unwrap();
    let second = observer
        .far_field_pressure(&flow, 0.8, distance + quarter_wavelength)
        .unwrap();
    close(
        second,
        (C64::new(0.0, 1.0) * first).scale(distance / (distance + quarter_wavelength)),
        first.abs(),
    );
}

#[test]
fn signed_higher_modes_match_independent_radial_fields_and_removable_limits() {
    let exterior = liquid();
    let radius = 0.025;
    let modes = 6;
    let basis = AxisymmetricMouthBasis::new(ModeCount(modes)).unwrap();
    let flow = [
        C64::new(0.8e-6, 0.1e-6),
        C64::new(-0.3e-6, 0.7e-6),
        C64::new(0.2e-6, -0.6e-6),
        C64::new(-0.4e-6, -0.2e-6),
        C64::new(0.1e-6, 0.3e-6),
        C64::new(-0.2e-6, 0.1e-6),
    ];
    let radial = gauss_legendre(96);
    for ka in [0.2, 4.3, 13.0, 23.0] {
        let omega = ka * exterior.sound_speed_m_s() / radius;
        let observer =
            FlangedModalObserver::from_material(&exterior, radius, omega, modes).unwrap();
        for mu in [0.0, 0.3, 0.8, 1.0] {
            let transverse = ka * (1.0_f64 - mu * mu).sqrt();
            let mut transform = C64::ZERO;
            for (mode, &u) in flow.iter().enumerate() {
                let integral: f64 = radial
                    .iter()
                    .map(|&(eta, weight)| {
                        weight
                            * eta
                            * basis.shape(mode, eta).unwrap()
                            * crate::modal::bessel_j0_certified(transverse * eta)
                    })
                    .sum();
                transform = transform + u.scale(integral);
            }
            let expected = C64::new(0.0, -exterior.density_kg_m3() * omega / PI) * transform;
            close(
                observer.far_field(&flow, mu).unwrap().amplitude_pa_m,
                expected,
                exterior.density_kg_m3() * omega * 1e-6,
            );
        }
        let mut higher_only = flow;
        higher_only[0] = C64::ZERO;
        assert_eq!(
            observer
                .far_field(&higher_only, 1.0)
                .unwrap()
                .amplitude_pa_m,
            C64::ZERO,
            "zero-net-flow radial modes have an exact on-axis null"
        );
        assert!(
            observer
                .far_field(&higher_only, 0.3)
                .unwrap()
                .power_per_steradian_w
                > 0.0
        );
    }
    // Choose a receiver direction at each removable Hankel singularity.
    let ka = 23.0;
    let omega = ka * exterior.sound_speed_m_s() / radius;
    let observer = FlangedModalObserver::from_material(&exterior, radius, omega, modes).unwrap();
    for mode in 1..modes {
        let alpha = basis.radial_eigenvalues()[mode];
        let mu = (1.0 - (alpha / ka).powi(2)).sqrt();
        let mut single = [C64::ZERO; 6];
        single[mode] = C64::from_re(1e-6);
        let expected = C64::new(
            0.0,
            -exterior.density_kg_m3() * omega / (2.0 * PI)
                * crate::modal::bessel_j0_certified(alpha).abs()
                * 1e-6,
        );
        close(
            observer.far_field(&single, mu).unwrap().amplitude_pa_m,
            expected,
            expected.abs(),
        );
    }
}

#[test]
fn coherent_hemispherical_field_power_equals_the_complete_radiation_matrix() {
    for exterior in [FluidAcousticState::try_new(1.2, 343.0).unwrap(), liquid()] {
        let radius = 0.025;
        let flow = [
            C64::new(0.8e-6, 0.1e-6),
            C64::new(-0.3e-6, 0.7e-6),
            C64::new(0.2e-6, -0.6e-6),
            C64::new(-0.4e-6, -0.2e-6),
        ];
        for ka in [0.2, 4.3, 8.1, 23.0] {
            let omega = ka * exterior.sound_speed_m_s() / radius;
            let observer =
                FlangedModalObserver::from_material(&exterior, radius, omega, flow.len()).unwrap();
            let matrix =
                flanged_modal_impedance_from_material(&exterior, radius, omega, flow.len())
                    .unwrap();
            let expected = matrix_power(&matrix, &flow);
            assert!(expected > 0.0);
            close(
                C64::from_re(integrated_power(&observer, &flow)),
                C64::from_re(expected),
                expected,
            );
            if ka == 4.3 {
                let diagonal: f64 = (0..flow.len())
                    .map(|i| 0.5 * flow[i].norm_sq() * matrix[i * flow.len() + i].re)
                    .sum();
                assert!((expected - diagonal).abs() > 1e-3 * expected);
            }
        }
    }
}

#[test]
fn driven_different_fluid_mouth_reaches_the_observer_with_the_same_work() {
    let interior = liquid();
    let exterior = FluidAcousticState::try_new(850.0, 1100.0).unwrap();
    let basis = AxisymmetricMouthBasis::new(ModeCount(3)).unwrap();
    let radius = basis.radial_eigenvalues()[1] / 128.0;
    let duct = Duct {
        segments: [radius, 0.045, 0.038]
            .iter()
            .map(|&radius| Segment::Cylinder {
                radius,
                length: 0.009,
            })
            .collect(),
    };
    let geometry = PreparedModalDuct::new(&duct, 3, 1).unwrap();
    let omega = 128.0 * interior.sound_speed_m_s();
    let reference = interior.characteristic_impedance_pa_s_m() / (PI * radius * radius);
    let pressure = |p, z| PortExcitation::Pressure {
        open_circuit_pressure_pa: p,
        impedance_pa_s_m3: z,
    };
    let sources = [
        pressure(C64::new(2.0, 0.3), C64::new(0.6, 0.2).scale(reference)),
        PortExcitation::VolumeFlow(C64::new(0.25, -0.1).scale(1.0 / reference)),
        pressure(C64::new(0.4, 0.8), C64::new(0.3, -0.2).scale(reference)),
    ];
    let observer =
        FlangedModalObserver::from_material(&exterior, geometry.mouth_radius_m(), omega, 3)
            .unwrap();
    let response = geometry
        .drive_sources_from_material_with_exterior(&interior, &exterior, omega, &sources)
        .unwrap();
    let radiated = integrated_power(&observer, &response.mouth_flow);
    assert!(radiated > 0.0);
    close(
        C64::from_re(radiated),
        C64::from_re(response.inlet.inlet_power_w),
        radiated,
    );
    close(
        C64::from_re(response.inlet.generator_power_w),
        C64::from_re(response.inlet.source_dissipation_w + radiated),
        response.inlet.generator_power_w.abs(),
    );
    let transport =
        FluidTransportState::try_new(interior.density_kg_m3(), 4200.0, 0.6, 0.001).unwrap();
    let lossy = ThermoviscousFluidState::try_new(interior, transport, 4100.0).unwrap();
    let response = geometry
        .drive_sources_from_material_with_exterior_and_loss(
            &lossy,
            &exterior,
            omega,
            LossModel::Bessel,
            &sources,
        )
        .unwrap();
    let radiated = integrated_power(&observer, &response.mouth_flow);
    assert!(radiated > 0.0 && response.inlet.inlet_power_w > radiated);
}

#[test]
fn immutable_observer_samples_the_original_material_once_and_retains_its_domain() {
    struct Counted {
        state: FluidAcousticState,
        calls: Cell<usize>,
    }
    impl FluidAcoustics for Counted {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            self.calls.set(self.calls.get() + 1);
            Ok(self.state)
        }
    }
    let material = Counted {
        state: liquid(),
        calls: Cell::new(0),
    };
    let observer = FlangedModalObserver::from_material(&material, 0.03, 40000.0, 2).unwrap();
    for mu in [0.0, 0.4, 1.0] {
        observer
            .far_field(&[C64::from_re(1e-6), C64::ZERO], mu)
            .unwrap();
        observer
            .far_field_pressure(&[C64::from_re(1e-6), C64::ZERO], mu, 10.0)
            .unwrap();
    }
    assert_eq!(material.calls.get(), 1);
    assert!(observer.far_field(&[C64::ZERO], 0.5).is_err());
    assert!(observer.far_field(&[C64::ZERO; 2], -0.1).is_err());
    assert!(
        observer
            .far_field_pressure(&[C64::ZERO; 2], 0.5, 0.0)
            .is_err()
    );
    assert!(matches!(
        FlangedModalObserver::from_material(&material, 0.03, 1e9, 2),
        Err(FluidDuctError::Duct(DuctError::RadiationKaTooLarge { .. }))
    ));
    struct Failed;
    impl FluidAcoustics for Failed {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            Err(MaterialError::Parameters {
                what: "original exterior refused its phase".into(),
            })
        }
    }
    assert!(matches!(
        FlangedModalObserver::from_material(&Failed, 0.03, 40000.0, 2),
        Err(FluidDuctError::Material(MaterialError::Parameters { what }))
            if what == "original exterior refused its phase"
    ));
}
