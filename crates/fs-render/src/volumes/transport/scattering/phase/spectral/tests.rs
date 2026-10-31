use super::*;
use core::f64::consts::PI;

fn dispersive() -> PhaseFunction {
    PhaseFunction::spectral_table(&[
        (380.0, PhaseFunction::henyey_greenstein(0.7).unwrap()),
        (780.0, PhaseFunction::henyey_greenstein(-0.6).unwrap()),
    ])
    .unwrap()
}

#[test]
fn g0_broadband_phase_uses_source_support_and_the_same_matched_packet_law() {
    let source = ContentHash([13; 32]);
    let phase = PhaseFunction::source_spectral_table(
        &[
            (250.0, PhaseFunction::Isotropic),
            (10_000.0, PhaseFunction::Rayleigh),
        ],
        source,
    )
    .unwrap();
    let PhaseFunction::SpectralTable(table) = &phase else {
        unreachable!()
    };
    assert_eq!(table.wavelength_bounds_nm(), (250.0, 10_000.0));
    assert_eq!(table.source_identity(), Some(source));
    let expected = 0.5
        * (PhaseFunction::Isotropic.eval_cos(0.3).unwrap()
            + PhaseFunction::Rayleigh.eval_cos(0.3).unwrap());
    assert_eq!(phase.eval_at(5125.0, 0.3).unwrap(), expected);
    for invalid in [0.0, 249.0, 10_001.0, f64::NAN] {
        assert_eq!(
            phase.eval_at(invalid, 0.0),
            Err(VolumeTransportError::Wavelength)
        );
    }
    let wavelengths = [250.0, 5125.0, 10_000.0];
    let incoming = Vec3::new(0.0, 0.0, 1.0);
    for u in [0.01, 0.2, 0.4, 0.6, 0.8, 0.99] {
        let sample = phase
            .sample_packet(&wavelengths, incoming, u, 0.25)
            .unwrap();
        let evaluated = phase.eval_packet(&wavelengths, sample.direction.z).unwrap();
        assert_eq!(sample.pdf, evaluated.pdf);
        for lane in 0..3 {
            assert_eq!(sample.weight[lane], evaluated.value[lane] / evaluated.pdf);
        }
    }
    assert!(PhaseFunction::Rayleigh.eval_at(150_000.0, 0.3).is_ok());
    assert_eq!(
        PhaseFunction::Rayleigh.eval_at(0.0, 0.3),
        Err(VolumeTransportError::Wavelength)
    );
    assert_eq!(
        dispersive().eval_at(1064.0, 0.3),
        Err(VolumeTransportError::Wavelength)
    );
}

#[test]
fn g0_spectral_interpolation_is_normalized_and_retains_endpoint_laws() {
    let phase = dispersive();
    let PhaseFunction::SpectralTable(table) = &phase else {
        unreachable!()
    };
    for cosine in [-1.0, -0.23, 0.0, 0.82, 1.0] {
        for (wavelength, endpoint) in table.knots() {
            assert_eq!(
                phase.eval_at(*wavelength, cosine).unwrap().to_bits(),
                endpoint.eval_cos(cosine).unwrap().to_bits()
            );
        }
    }
    let n = 32_768;
    for wavelength in [380.0, 417.0, 580.0, 669.0, 780.0] {
        let mut mass = 0.0;
        let mut moment = 0.0;
        for i in 0..n {
            let mu = 2.0 * (f64::from(i) + 0.5) / f64::from(n) - 1.0;
            let weight = phase.eval_at(wavelength, mu).unwrap() * 4.0 * PI / f64::from(n);
            mass += weight;
            moment += weight * mu;
        }
        let t = (wavelength - 380.0) / 400.0;
        assert!((mass - 1.0).abs() < 2e-7);
        assert!((moment - ((1.0 - t) * 0.7 - t * 0.6)).abs() < 2e-7);
    }
}

#[test]
fn g1_packet_sampling_preserves_each_wavelengths_energy_and_angular_moments() {
    let phase = dispersive();
    let wavelengths = [380.0, 497.0, 780.0];
    let incoming = Vec3::new(0.0, 0.0, 1.0);
    let mut mass = [0.0; 3];
    let mut first = [0.0; 3];
    let mut second = [0.0; 3];
    let n = 98_304;
    for i in 0..n {
        let sample = phase
            .sample_packet(
                &wavelengths,
                incoming,
                (f64::from(i) + 0.5) / f64::from(n),
                0.37,
            )
            .unwrap();
        let mu = incoming.dot(sample.direction);
        let evaluated = phase.eval_packet(&wavelengths, mu).unwrap();
        assert_eq!(sample.pdf, evaluated.pdf);
        assert!((sample.direction.dot(sample.direction) - 1.0).abs() < 1e-14);
        for lane in 0..3 {
            let contribution = sample.weight[lane] / f64::from(n);
            mass[lane] += contribution;
            first[lane] += contribution * mu;
            second[lane] += contribution * mu * mu;
        }
    }
    for lane in 0..3 {
        let t = (wavelengths[lane] - 380.0) / 400.0;
        let expected_first = (1.0 - t) * 0.7 - t * 0.6;
        let expected_second =
            ((1.0 - t) * (1.0 + 2.0 * 0.7_f64.powi(2)) + t * (1.0 + 2.0 * 0.6_f64.powi(2))) / 3.0;
        assert!((mass[lane] - 1.0).abs() < 1e-4, "{mass:?}");
        assert!((first[lane] - expected_first).abs() < 1e-4, "{first:?}");
        assert!((second[lane] - expected_second).abs() < 1e-4, "{second:?}");
    }
}

#[test]
fn g0_packet_mixture_covers_disjoint_spectral_angular_support() {
    let phase = PhaseFunction::spectral_table(&[
        (
            380.0,
            PhaseFunction::tabulated(&[(-1.0, 4.0), (-0.2, 0.0), (1.0, 0.0)]).unwrap(),
        ),
        (
            780.0,
            PhaseFunction::tabulated(&[(-1.0, 0.0), (0.2, 0.0), (1.0, 3.0)]).unwrap(),
        ),
    ])
    .unwrap();
    let wavelengths = [380.0, 780.0];
    for u in [0.03, 0.2, 0.6, 0.95] {
        let sample = phase
            .sample_packet(&wavelengths, Vec3::new(0.0, 0.0, 1.0), u, 0.6)
            .unwrap();
        assert!(sample.pdf > 0.0);
        let active = usize::from(u >= 0.5);
        assert_eq!(sample.weight[1 - active], 0.0);
        assert!((sample.weight[active] - 2.0).abs() < 1e-14);
    }
}

#[test]
fn g5_achromatic_packet_sampling_preserves_original_directions_and_pdf_bits() {
    let incoming = Vec3::new(0.6, 0.0, 0.8);
    for phase in [
        PhaseFunction::Isotropic,
        PhaseFunction::Rayleigh,
        PhaseFunction::henyey_greenstein(0.83).unwrap(),
    ] {
        for u in [0.0, 0.17, 0.5, 0.96, fs_math::next_down(1.0)] {
            let scalar = phase.sample(incoming, u, 0.37).unwrap();
            let packet = phase
                .sample_packet(&[400.0, 500.0, 700.0], incoming, u, 0.37)
                .unwrap();
            assert_eq!(scalar.direction, packet.direction);
            assert_eq!(scalar.pdf.to_bits(), packet.pdf.to_bits());
            assert_eq!(packet.weight, [1.0; 3]);
        }
    }
}

#[test]
fn g1_depolarized_molecular_rayleigh_has_analytic_phase_and_second_moment() {
    for d in [0.0, 0.03, 0.3, 0.8] {
        let phase = PhaseFunction::rayleigh_depolarized(d).unwrap();
        let n = 32_768;
        let mut moment = 0.0;
        for i in 0..n {
            let sample = phase
                .sample(
                    Vec3::new(0.0, 0.0, 1.0),
                    (f64::from(i) + 0.5) / f64::from(n),
                    0.25,
                )
                .unwrap();
            let mu = sample.direction.z;
            let expected = 3.0 * ((1.0 + d) + (1.0 - d) * mu * mu) / (8.0 * PI * (2.0 + d));
            assert!((sample.pdf / expected - 1.0).abs() < 2e-14);
            moment += mu * mu / f64::from(n);
        }
        assert!((moment - (4.0 + d) / (5.0 * (2.0 + d))).abs() < 5e-5);
    }
    for d in [-0.01, 6.0 / 7.0, 1.0, f64::NAN, f64::INFINITY] {
        assert!(PhaseFunction::rayleigh_depolarized(d).is_err());
    }
}

#[test]
fn g0_spectral_phase_refuses_gray_use_missing_support_and_recursive_tables() {
    let phase = dispersive();
    let direction = Vec3::new(0.0, 0.0, 1.0);
    assert!(phase.eval_cos(0.0).is_err());
    assert!(phase.sample(direction, 0.2, 0.4).is_err());
    assert!(phase.eval_packet(&[], 0.0).is_err());
    for wavelength in [0.0, 379.0, 781.0, f64::NAN] {
        assert!(phase.eval_at(wavelength, 0.0).is_err());
    }
    for knots in [
        vec![
            (400.0, PhaseFunction::Rayleigh),
            (780.0, PhaseFunction::Rayleigh),
        ],
        vec![(380.0, PhaseFunction::Rayleigh), (780.0, phase)],
        vec![
            (380.0, PhaseFunction::Rayleigh),
            (380.0, PhaseFunction::Rayleigh),
        ],
    ] {
        assert!(PhaseFunction::spectral_table(&knots).is_err());
    }
}
