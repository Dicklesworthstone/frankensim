//! G1/G3 finite detectors observe the original partial-coherence ensemble.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumAperture, AngularSpectrumApertureLimits, AngularSpectrumApertureReading,
};
use fs_math::dd::Dd;

fn detector_limits() -> AngularSpectrumApertureLimits {
    AngularSpectrumApertureLimits {
        maximum_mode_pairs: 4096,
        maximum_workspace_bytes: 4 << 20,
    }
}
fn bank<'b, 's, 'm>(
    factors: &'b [PlanarAngularSpectrum<'s, 'm>],
) -> AngularSpectrumCoherence<'b, 's, 'm> {
    AngularSpectrumCoherence::try_new(&factors[0], factors, coherence_limits(), || false).unwrap()
}
fn same_reading(a: AngularSpectrumApertureReading, b: AngularSpectrumApertureReading) {
    let flux_scale = a
        .mean_poynting_w_m2
        .iter()
        .chain(&b.mean_poynting_w_m2)
        .fold(1e-30_f64, |scale, v| scale.max(v.abs()));
    for c in 0..3 {
        near(
            a.mean_poynting_w_m2[c],
            b.mean_poynting_w_m2[c],
            2e-11 * flux_scale,
        );
    }
    close(a.normal_power_w, b.normal_power_w);
}

#[test]
fn g1_finite_detector_resolves_coherence_visibility_with_correct_vector_flux() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let g = grid();
    let kx = core::f64::consts::TAU / (g.width as f64 * g.spacing_m[0]);
    let kz = det::sqrt(1.0 - kx * kx);
    for rho in [0.0, 0.6, 1.0] {
        let mut a = Vec::new();
        let mut b = Vec::new();
        for x in 0..g.width {
            let p = phase(-kx * x as f64 * g.spacing_m[0]);
            a.push([C64::ZERO, C64::ONE + p.scale(rho)]);
            b.push([C64::ZERO, p.scale(det::sqrt(1.0 - rho * rho))]);
        }
        let factors = [prepare(&stack, &a), prepare(&stack, &b)];
        let bank = bank(&factors);
        for center in [0.13, core::f64::consts::PI / kx] {
            for width in [0.5, 6.7] {
                let aperture =
                    AngularSpectrumAperture::try_new([center - width / 2.0, -0.2], [width, 1.7])
                        .unwrap();
                let actual = bank
                    .measure_aperture(Side::Transmitted, 0.0, aperture, detector_limits(), || {
                        false
                    })
                    .unwrap();
                // Independent continuous rectangle integral of the two TE
                // Maxwell waves: Hx=-kz*Ey and Hz=kx*Ey for omega*mu=1.
                let half_phase = kx * width / 2.0;
                let fringe = rho * det::cos(kx * center) * det::sin(half_phase) / half_phase;
                close(actual.mean_poynting_w_m2[0], 0.5 * kx * (1.0 + fringe));
                assert_eq!(actual.mean_poynting_w_m2[1], 0.0);
                close(
                    actual.mean_poynting_w_m2[2],
                    0.5 * (1.0 + kz) * (1.0 + fringe),
                );
                close(
                    actual.normal_power_w,
                    0.5 * (1.0 + kz) * (1.0 + fringe) * aperture.area_m2(),
                );
            }
        }
        let cell = AngularSpectrumAperture::try_new(
            [0.0, 0.0],
            [
                g.width as f64 * g.spacing_m[0],
                g.height as f64 * g.spacing_m[1],
            ],
        )
        .unwrap();
        close(
            bank.measure_aperture(Side::Transmitted, 0.0, cell, detector_limits(), || false)
                .unwrap()
                .normal_power_w,
            bank.power().transmitted_w,
        );
    }
}

#[test]
fn g3_finite_detector_preserves_complex_factor_basis_and_signed_reflected_power() {
    let medium = material(E, M);
    let film = material(4.0 * E, M);
    let exit = material(2.25 * E, M);
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&exit),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.17,
        }],
    )
    .unwrap();
    let (mut a, mut b, mut ua, mut ub) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let j = C64::new(0.0, 1.0);
    for x in 0..grid().width {
        let p = phase(-core::f64::consts::TAU * x as f64 / grid().width as f64);
        let av = [C64::new(0.8, 0.2) + p.scale(0.3), j];
        let bv = [p * C64::new(0.2, -0.6), p.scale(0.4) + C64::new(0.2, 0.1)];
        ua.push(std::array::from_fn(|c| {
            (av[c] + j * bv[c]).scale(core::f64::consts::FRAC_1_SQRT_2)
        }));
        ub.push(std::array::from_fn(|c| {
            (j * av[c] + bv[c]).scale(core::f64::consts::FRAC_1_SQRT_2)
        }));
        a.push(av);
        b.push(bv);
    }
    let factors = [prepare(&stack, &a), prepare(&stack, &b)];
    let rotated = [prepare(&stack, &ua), prepare(&stack, &ub)];
    let original = bank(&factors);
    let mixed = bank(&rotated);
    let aperture = AngularSpectrumAperture::try_new([-0.3, -1.2], [2.0, 1.4]).unwrap();
    for side in [Side::Reflected, Side::TotalIncidentSide, Side::Transmitted] {
        let actual = original
            .measure_aperture(side, 0.21, aperture, detector_limits(), || false)
            .unwrap();
        same_reading(
            actual,
            mixed
                .measure_aperture(side, 0.21, aperture, detector_limits(), || false)
                .unwrap(),
        );
        let readings = factors.each_ref().map(|f| {
            f.measure_aperture(side, 0.21, aperture, detector_limits(), || false)
                .unwrap()
        });
        same_reading(
            actual,
            AngularSpectrumApertureReading {
                mean_poynting_w_m2: std::array::from_fn(|c| {
                    readings[0].mean_poynting_w_m2[c] + readings[1].mean_poynting_w_m2[c]
                }),
                normal_power_w: readings[0].normal_power_w + readings[1].normal_power_w,
            },
        );
        if side == Side::Reflected {
            assert!(actual.mean_poynting_w_m2[2] < 0.0);
            assert!(actual.normal_power_w < 0.0);
        }
    }
}

#[test]
fn g0_g4_detector_admits_aggregate_work_peak_storage_and_darkness_before_work() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let input = vec![[C64::ONE, C64::new(0.0, 0.5)]; grid().width];
    let factors = [prepare(&stack, &input), prepare(&stack, &input)];
    let bank = bank(&factors);
    let g = grid();
    let n = g.width * g.height;
    let exact = AngularSpectrumApertureLimits {
        maximum_mode_pairs: factors.len() * n * (g.width + g.height),
        maximum_workspace_bytes: factors[0].workspace_bytes()
            + (g.width * g.width + g.height * g.height) * core::mem::size_of::<[Dd; 2]>()
            + n * core::mem::size_of::<[[Dd; 2]; 3]>(),
    };
    let aperture = AngularSpectrumAperture::try_new([-0.2, 0.4], [0.7, 1.3]).unwrap();
    let before = bank
        .measure_aperture(Side::Transmitted, 0.2, aperture, exact, || false)
        .unwrap();
    let mut polls = 0;
    assert!(
        bank.measure_aperture(
            Side::Transmitted,
            0.2,
            aperture,
            AngularSpectrumApertureLimits {
                maximum_mode_pairs: exact.maximum_mode_pairs - 1,
                ..exact
            },
            || {
                polls += 1;
                false
            },
        )
        .is_err()
    );
    assert_eq!(
        polls, 1,
        "reject aggregate work before preparing any kernels"
    );
    assert!(
        bank.measure_aperture(
            Side::Transmitted,
            0.2,
            aperture,
            AngularSpectrumApertureLimits {
                maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
                ..exact
            },
            || false,
        )
        .is_err()
    );
    let dark = AngularSpectrumCoherence::try_new(
        &factors[0],
        &[],
        AngularSpectrumCoherenceLimits {
            maximum_modes: 0,
            maximum_point_mode_pairs: 0,
            maximum_fft_work: 0,
            ..coherence_limits()
        },
        || false,
    )
    .unwrap();
    let none = AngularSpectrumApertureLimits {
        maximum_mode_pairs: 0,
        maximum_workspace_bytes: 0,
    };
    assert_eq!(
        dark.measure_aperture(Side::Transmitted, 0.2, aperture, none, || false)
            .unwrap(),
        AngularSpectrumApertureReading {
            mean_poynting_w_m2: [0.0; 3],
            normal_power_w: 0.0,
        }
    );
    for (source, limits) in [(&bank, exact), (&dark, none)] {
        for depth in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                source
                    .measure_aperture(Side::Transmitted, depth, aperture, limits, || false)
                    .is_err()
            );
        }
        assert!(matches!(
            source.measure_aperture(Side::Transmitted, 0.2, aperture, limits, || true),
            Err(AngularSpectrumError::Cancelled)
        ));
    }
    let mut polls = 0;
    assert!(matches!(
        bank.measure_aperture(Side::Transmitted, 0.2, aperture, exact, || {
            polls += 1;
            polls > 30
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(
        bank.measure_aperture(Side::Transmitted, 0.2, aperture, exact, || false)
            .unwrap(),
        before
    );
}
