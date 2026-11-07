//! Physical angular/spectral responses through the native collision and adjoint owners.

use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::volumes::transport::OpticalMixtureBudget;
use crate::volumes::transport::scattering::tests::with_cx;
use crate::volumes::transport::scattering::{
    MassScatteringSpectrum, ScatteringMixtureComponent, SphericalPhaseTable, SphericalPhaseValue,
};
use fs_geom::Point3;
use fs_math::c64::C64;

fn beam(direction: Vec3, axis: Vec3, wavelength: f64, power: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(direction, axis).unwrap(),
        wavelength,
        DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::CLEAR,
            GlassProvenance::Custom,
        ),
        power,
    )
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-280),
        "{a:e} != {b:e}"
    );
}

fn phase(bias: f64, retardance: f64) -> PhaseFunction {
    let knots: Vec<_> = (0..=64)
        .map(|i| {
            let c = -1.0 + i as f64 / 32.0;
            let power = 1.0 + bias * c;
            let angle = retardance * (1.0 - c * c);
            (
                c,
                SphericalPhaseValue::try_new(
                    [power, power * c * c],
                    C64::new(power * c * angle.cos(), power * c * angle.sin()),
                )
                .unwrap(),
            )
        })
        .collect();
    PhaseFunction::SphericalTable(SphericalPhaseTable::try_new(&knots).unwrap())
}

#[test]
fn spectral_mixtures_share_scattering_weights_and_native_wavelength_pdf() {
    with_cx(|cx, _| {
        let left = phase(-0.6, 0.8);
        let right = phase(0.7, -1.2);
        let spectral =
            PhaseFunction::spectral_table(&[(380.0, left.clone()), (780.0, right.clone())])
                .unwrap();
        let components = [
            ScatteringMixtureComponent {
                mass_weight: 0.3,
                scattering: MassScatteringSpectrum::table(&[(380.0, 1.0), (780.0, 5.0)]).unwrap(),
                phase: spectral.clone(),
            },
            ScatteringMixtureComponent {
                mass_weight: 0.7,
                scattering: MassScatteringSpectrum::table(&[(380.0, 6.0), (780.0, 2.0)]).unwrap(),
                phase: PhaseFunction::Rayleigh,
            },
        ];
        let mixture = PhaseFunction::scattering_mixture(
            cx,
            &components,
            OpticalMixtureBudget { max_components: 2 },
        )
        .unwrap();
        mixture.admit_polarized().unwrap();
        for wavelength in [380.0, 487.0, 550.0, 780.0] {
            let t = (wavelength - 380.0) / 400.0;
            let a = 0.3 * (1.0 + 4.0 * t);
            let b = 0.7 * (6.0 - 4.0 * t);
            let weight = a / (a + b);
            for c in [-0.8, 0.0, 0.37, 1.0] {
                let x = left.mueller_at(wavelength, c).unwrap().0;
                let y = right.mueller_at(wavelength, c).unwrap().0;
                let r = PhaseFunction::Rayleigh.mueller_at(wavelength, c).unwrap().0;
                let actual = mixture.mueller_at(wavelength, c).unwrap().0;
                for row in 0..4 {
                    for column in 0..4 {
                        let expected = weight * ((1.0 - t) * x[row][column] + t * y[row][column])
                            + (1.0 - weight) * r[row][column];
                        assert!((actual[row][column] - expected).abs() < 3e-16);
                    }
                }
                close(actual[0][0], mixture.eval_at(wavelength, c).unwrap(), 3e-14);
            }
            let input = beam(
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 0.0),
                wavelength,
                Stokes::try_new(1.0, 0.2, -0.4, 0.7).unwrap(),
            );
            for samples in [[0.01, 0.23], [0.5, 0.71], [0.99, 0.81]] {
                let sample = mixture.sample_polarized(input, samples).unwrap();
                let cosine = input
                    .frame()
                    .direction()
                    .dot(sample.beam.frame().direction());
                close(
                    sample.pdf,
                    mixture.eval_at(wavelength, cosine).unwrap(),
                    4e-14,
                );
            }
        }
        assert!(mixture.mueller_at(379.0, 0.0).is_err());
    });
}

#[test]
fn particle_sampling_conserves_every_input_polarization_and_reference_rotation() {
    let phase = phase(0.7, 1.1);
    for state in [
        Stokes::unpolarized(1.0),
        Stokes::try_new(1.0, 1.0, 0.0, 0.0).unwrap(),
        Stokes::try_new(1.0, 0.0, 0.0, 1.0).unwrap(),
        Stokes::try_new(1.0, 0.3, -0.4, 0.5).unwrap(),
    ] {
        let input = beam(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            550.0,
            state,
        );
        let rotated = input.with_reference_axis(Vec3::new(0.3, 0.7, 0.0)).unwrap();
        let mut sum = 0.0;
        for i in 0..128 {
            for j in 0..32 {
                let samples = [(i as f64 + 0.5) / 128.0, (j as f64 + 0.5) / 32.0];
                let a = phase.sample_polarized(input, samples).unwrap();
                let b = phase.sample_polarized(rotated, samples).unwrap();
                let x = a.beam.power();
                let y = b.beam.power();
                for (a, b) in [x.i, x.q, x.u, x.v].into_iter().zip([y.i, y.q, y.u, y.v]) {
                    assert!((a - b).abs() < 5e-14);
                }
                assert!(x.q.hypot(x.u).hypot(x.v) <= x.i * (1.0 + 8e-15));
                sum += x.i;
            }
        }
        close(sum / (128.0 * 32.0), 1.0, 2e-13);
    }
}

#[test]
fn particle_connections_match_full_adjoint_and_reciprocal_reversal() {
    let phase = phase(-0.4, 1.7);
    let axis = Vec3::new(0.0, 1.0, 0.0);
    let detector = Stokes::try_new(1.0, 0.2, 0.3, -0.7).unwrap();
    let source = Stokes::try_new(1.0, -0.4, 0.6, 0.3).unwrap();
    let camera = beam(Vec3::new(0.0, 0.0, -1.0), axis, 550.0, detector);
    for c in [-0.9_f64, 0.0, 0.3, 0.95] {
        let direction = Vec3::new((1.0 - c * c).sqrt(), 0.0, c);
        let frame = PolarizationFrame::try_new(direction, axis).unwrap();
        let output = phase.mueller_at(550.0, c).unwrap().apply(source);
        let expected =
            output.i + detector.q * output.q - detector.u * output.u + detector.v * output.v;
        let actual = phase
            .polarized_connection_factor(camera, frame, source)
            .unwrap();
        close(actual, expected, 4e-14);
        let reversed = beam(direction, axis, 550.0, source);
        close(
            phase
                .polarized_connection_factor(reversed, camera.frame(), detector)
                .unwrap(),
            actual,
            4e-14,
        );
        let rotated = camera
            .with_reference_axis(Vec3::new(0.7, 0.2, 0.0))
            .unwrap();
        close(
            phase
                .polarized_connection_factor(rotated, frame, source)
                .unwrap(),
            actual,
            4e-14,
        );
    }
}

#[test]
fn particle_weak_channels_and_true_extinction_survive_connections() {
    let camera = beam(
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        550.0,
        Stokes::unpolarized(1.0),
    );
    let source =
        PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap();
    for weak in [0.0, 2e-240] {
        let value = SphericalPhaseValue::try_new([2.0, weak], C64::ZERO).unwrap();
        let phase = PhaseFunction::SphericalTable(
            SphericalPhaseTable::try_new(&[(-1.0, value), (1.0, value)]).unwrap(),
        );
        let expected = phase
            .spherical_value_at(550.0, 0.0)
            .unwrap()
            .powers_per_sr()[1];
        let factor = phase
            .polarized_connection_factor(
                camera,
                source,
                Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap(),
            )
            .unwrap();
        if weak == 0.0 {
            assert_eq!(factor, 0.0);
        } else {
            assert!(factor > 0.0);
            close(factor, expected, 3e-14);
        }
    }
}

#[test]
fn scalar_only_laws_and_zero_scattering_mixtures_do_not_invent_stokes_response() {
    let beam = beam(
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        550.0,
        Stokes::unpolarized(1.0),
    );
    for phase in [
        PhaseFunction::Isotropic,
        PhaseFunction::henyey_greenstein(0.5).unwrap(),
    ] {
        assert_eq!(
            phase.sample_polarized(beam, [0.3, 0.7]),
            Err(PolarizedScatteringError::UnsupportedPhase)
        );
    }
    with_cx(|cx, _| {
        let components = [ScatteringMixtureComponent {
            mass_weight: 0.0,
            scattering: MassScatteringSpectrum::table(&[(380.0, 1.0), (780.0, 1.0)]).unwrap(),
            phase: PhaseFunction::Isotropic,
        }];
        let mixture = PhaseFunction::scattering_mixture(
            cx,
            &components,
            OpticalMixtureBudget { max_components: 1 },
        )
        .unwrap();
        mixture.admit_polarized().unwrap();
        assert_eq!(
            mixture.apply_polarized_at(550.0, 0.0, Stokes::unpolarized(1.0)),
            Err(PolarizedScatteringError::UnsupportedPhase)
        );
    });
    let large = PhaseFunction::Rayleigh
        .apply_polarized_at(550.0, 0.4, Stokes::try_new(1e308, 1e308, 0.0, 0.0).unwrap())
        .unwrap();
    assert!(large.i.is_finite() && large.i > 0.0);
}
