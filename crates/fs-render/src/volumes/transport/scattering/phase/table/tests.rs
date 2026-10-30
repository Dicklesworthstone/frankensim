use super::*;
use crate::volumes::transport::scattering::PhaseFunction;
use fs_geom::Vec3;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn g1_nonuniform_linear_tables_have_unit_power_and_analytic_moments() {
    const COUNT: usize = 16_384;
    for slope in [-1.0, -0.7, 0.0, 0.7, 1.0] {
        let knots = [-1.0, -0.37, 0.11, 0.72, 1.0].map(|mu| (mu, 1.0 + slope * mu));
        let table = TabulatedPhase::try_new(&knots).unwrap();
        close(table.mean_cosine(), slope / 3.0, 3.0e-16);
        let phase = PhaseFunction::Tabulated(table.clone());
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let (mut integral, mut mean, mut second) = (0.0, 0.0, 0.0);
        for index in 0..COUNT {
            let u = (index as f64 + 0.5) / COUNT as f64;
            let cosine = 2.0 * u - 1.0;
            let expected = (1.0 + slope * cosine) / (4.0 * core::f64::consts::PI);
            close(table.eval_cos(cosine).unwrap(), expected, 6.0e-17);
            integral +=
                phase.eval_cos(cosine).unwrap() * 4.0 * core::f64::consts::PI / COUNT as f64;
            let sample = phase
                .sample(axis, u, (u * 1.618_033_988_749_895).fract())
                .unwrap();
            let c = axis.dot(sample.direction);
            // Independent CDF of (1+a*mu)/2, valid across every table knot.
            close(0.5 * (c + 1.0) + 0.25 * slope * (c * c - 1.0), u, 9.0e-16);
            assert_eq!(sample.pdf, phase.eval_cos(c).unwrap());
            close(sample.direction.dot(sample.direction), 1.0, 8.0e-16);
            mean += c / COUNT as f64;
            second += c * c / COUNT as f64;
        }
        close(integral, 1.0, 4.0e-14);
        close(mean, slope / 3.0, 2.0e-7);
        close(second, 1.0 / 3.0, 4.0e-7);
    }
}

#[test]
fn g1_dark_caps_and_gaps_keep_zero_density_and_sample_only_supported_lobes() {
    let phase = PhaseFunction::tabulated(&[
        (-1.0, 0.0),
        (-0.6, 0.0),
        (-0.4, 1.0),
        (-0.2, 0.0),
        (0.2, 0.0),
        (0.4, 2.0),
        (0.6, 0.0),
        (1.0, 0.0),
    ])
    .unwrap();
    for cosine in [-1.0, -0.8, 0.0, 0.8, 1.0] {
        assert_eq!(phase.eval_cos(cosine).unwrap(), 0.0);
    }
    let diagonal = 1.0 / det::sqrt(3.0);
    for axis in [
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(diagonal, diagonal, diagonal),
    ] {
        for u in [0.0, 1.0 / 3.0, 0.5, 1.0_f64.next_down()] {
            let sample = phase.sample(axis, u, 0.173).unwrap();
            assert!(sample.pdf > 0.0 && sample.pdf.is_finite());
            close(
                sample.pdf,
                phase.eval_cos(axis.dot(sample.direction)).unwrap(),
                3.0e-15,
            );
        }
    }
    const COUNT: usize = 12_288;
    let mut backwards = 0;
    let mut mean = 0.0;
    for index in 0..COUNT {
        let sample = phase
            .sample(
                Vec3::new(0.0, 0.0, 1.0),
                (index as f64 + 0.5) / COUNT as f64,
                0.37,
            )
            .unwrap();
        let cosine = sample.direction.z;
        assert!((-0.6..-0.2).contains(&cosine) || (0.2..0.6).contains(&cosine));
        backwards += usize::from(cosine < 0.0);
        mean += cosine / COUNT as f64;
    }
    assert_eq!(backwards, COUNT / 3);
    close(mean, 0.4 / 3.0, 1.0e-12);
}

#[test]
fn g1_hg_table_refinement_converges_to_known_asymmetry() {
    let g = 0.75_f64;
    let reference = |cosine: f64| (1.0 - g * g) / (1.0 + g * g - 2.0 * g * cosine).powf(1.5);
    let mut errors = Vec::new();
    for intervals in [128, 256, 512] {
        let knots: Vec<_> = (0..=intervals)
            .map(|index| {
                let cosine = -1.0 + 2.0 * index as f64 / intervals as f64;
                (cosine, reference(cosine))
            })
            .collect();
        let table = TabulatedPhase::try_new(&knots).unwrap();
        errors.push((table.mean_cosine() - g).abs());
    }
    assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.3 * errors[1], "{errors:?}");
    assert!(errors[2] < 3.0e-4, "{errors:?}");
}

#[test]
fn g3_normalization_handles_large_small_and_common_scaled_values() {
    let ordinary = TabulatedPhase::try_new(&[(-1.0, 2.0), (-0.25, 1.0), (1.0, 2.0)]).unwrap();
    let rescaled = TabulatedPhase::try_new(&[(-1.0, 16.0), (-0.25, 8.0), (1.0, 16.0)]).unwrap();
    let large =
        TabulatedPhase::try_new(&[(-1.0, f64::MAX), (-0.25, f64::MAX * 0.5), (1.0, f64::MAX)])
            .unwrap();
    assert_eq!(ordinary, rescaled);
    assert_eq!(ordinary, large);
    let tiny =
        TabulatedPhase::try_new(&[(-1.0, f64::MIN_POSITIVE), (1.0, f64::MIN_POSITIVE)]).unwrap();
    for cosine in [-1.0, -0.125, 0.0, 0.7, 1.0] {
        close(
            tiny.eval_cos(cosine).unwrap(),
            1.0 / (4.0 * core::f64::consts::PI),
            2.0e-17,
        );
    }
}

#[test]
fn g0_invalid_or_unresolvable_angular_data_refuses() {
    let invalid: &[&[(f64, f64)]] = &[
        &[],
        &[(-1.0, 1.0)],
        &[(-0.9, 1.0), (1.0, 1.0)],
        &[(-1.0, 1.0), (0.9, 1.0)],
        &[(-1.0, 1.0), (0.0, 1.0), (0.0, 2.0), (1.0, 1.0)],
        &[(-1.0, 1.0), (0.2, 1.0), (-0.2, 1.0), (1.0, 1.0)],
        &[(-1.0, 1.0), (f64::NAN, 1.0), (1.0, 1.0)],
        &[(-1.0, 1.0), (1.0, -1.0)],
        &[(-1.0, f64::INFINITY), (1.0, 1.0)],
        &[(-1.0, 0.0), (1.0, 0.0)],
    ];
    for knots in invalid {
        assert_eq!(
            TabulatedPhase::try_new(knots),
            Err(VolumeTransportError::Phase)
        );
    }
    let maximum: Vec<_> = (0..MAX_PHASE_KNOTS)
        .map(|index| {
            (
                -1.0 + 2.0 * index as f64 / (MAX_PHASE_KNOTS - 1) as f64,
                1.0,
            )
        })
        .collect();
    assert_eq!(
        TabulatedPhase::try_new(&maximum)
            .unwrap()
            .knots_per_sr()
            .len(),
        MAX_PHASE_KNOTS
    );
    assert_eq!(
        TabulatedPhase::try_new(&vec![(-1.0, 1.0); MAX_PHASE_KNOTS + 1]),
        Err(VolumeTransportError::Phase)
    );
    assert_eq!(
        TabulatedPhase::try_new(&[(-1.0, 1.0), (0.0, 1.0), (0.1, 0.0), (1.0, 1.0e-40)]),
        Err(VolumeTransportError::Unrepresentable),
        "a positive tail cannot disappear into a rounded cumulative probability"
    );
    assert_eq!(
        TabulatedPhase::try_new(&[(-1.0, f64::MAX), (1.0, f64::from_bits(1))]),
        Err(VolumeTransportError::Unrepresentable)
    );
    let phase = PhaseFunction::tabulated(&[(-1.0, 1.0), (1.0, 1.0)]).unwrap();
    for u in [-0.1, 1.0, f64::NAN] {
        assert_eq!(
            phase.sample(Vec3::new(0.0, 0.0, 1.0), u, 0.5),
            Err(VolumeTransportError::Phase)
        );
    }
    assert_eq!(
        phase.eval_cos(f64::INFINITY),
        Err(VolumeTransportError::Phase)
    );
}

#[test]
fn g5_table_identity_tracks_the_normalized_angular_law() {
    let identity = |phase: &PhaseFunction| {
        let mut hash = DomainHasher::new("org.frankensim.render.phase-test.v1");
        phase.hash_into(&mut hash);
        hash.finalize()
    };
    let first = PhaseFunction::tabulated(&[(-1.0, 1.0), (0.0, 3.0), (1.0, 2.0)]).unwrap();
    let scaled = PhaseFunction::tabulated(&[(-1.0, 8.0), (0.0, 24.0), (1.0, 16.0)]).unwrap();
    let changed = PhaseFunction::tabulated(&[(-1.0, 1.0), (0.0, 2.0), (1.0, 3.0)]).unwrap();
    assert_eq!(identity(&first), identity(&first.clone()));
    assert_eq!(identity(&first), identity(&scaled));
    assert_ne!(identity(&first), identity(&changed));
    assert_ne!(identity(&first), identity(&PhaseFunction::Rayleigh));
}
