//! Independent two-point interference oracles and original-owner parity.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumCorrelationLimits, AngularSpectrumMutualCoherence, AngularSpectrumObservation,
};

fn observation(side: Side, distance: f64, x: f64, y: f64) -> AngularSpectrumObservation {
    AngularSpectrumObservation {
        side,
        distance_from_face_m: distance,
        position_m: [x, y],
    }
}
fn query_limits(bank: &AngularSpectrumCoherence<'_, '_, '_>) -> AngularSpectrumCorrelationLimits {
    AngularSpectrumCorrelationLimits {
        maximum_mode_evaluations: bank.correlation_mode_evaluations().unwrap(),
        maximum_workspace_bytes: bank.correlation_workspace_bytes(),
    }
}
fn ey() -> [C64; 3] {
    [C64::ZERO, C64::ONE, C64::ZERO]
}
fn tensor_near(a: [[C64; 3]; 3], b: [[C64; 3]; 3], tolerance: f64) {
    for (a, b) in a.into_iter().flatten().zip(b.into_iter().flatten()) {
        complex_near(a, b, tolerance);
    }
}
fn same(a: &AngularSpectrumMutualCoherence, b: &AngularSpectrumMutualCoherence) {
    tensor_near(a.electric_at_a_v2_m2(), b.electric_at_a_v2_m2(), 3e-12);
    tensor_near(a.electric_at_b_v2_m2(), b.electric_at_b_v2_m2(), 3e-12);
    tensor_near(a.electric_cross_v2_m2(), b.electric_cross_v2_m2(), 3e-12);
}

#[test]
fn g1_identical_arm_intensities_have_different_spatial_coherence_and_fringe_visibility() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let length = grid().width as f64 * grid().spacing_m[0];
    let a = observation(Side::Transmitted, 0.0, length / 4.0, 0.0);
    let b = observation(Side::Transmitted, 0.0, 3.0 * length / 4.0, 0.0);
    for rho in [0.0_f64, 0.6, 1.0] {
        let mut first = Vec::new();
        let mut second = Vec::new();
        for x in 0..grid().width {
            let p = phase(-core::f64::consts::TAU * x as f64 / grid().width as f64);
            first.push([C64::ZERO, C64::ONE + p.scale(rho)]);
            second.push([C64::ZERO, p.scale(det::sqrt(1.0 - rho * rho))]);
        }
        let factors = [prepare(&stack, &first), prepare(&stack, &second)];
        let bank =
            AngularSpectrumCoherence::try_new(&factors[0], &factors, coherence_limits(), || false)
                .unwrap();
        let limits = query_limits(&bank);
        assert_eq!(limits.maximum_mode_evaluations, 2 * factors.len() * grid().width);
        let mutual = bank.mutual_coherence(a, b, limits, || false).unwrap();
        let projected = mutual.project([ey(); 2]).unwrap();
        // Independently: B1(a)=1-i*rho, B1(b)=1+i*rho,
        // B2(a)=-i*sqrt(1-rho²), B2(b)=+i*sqrt(1-rho²).
        // Both arm intensities are 1; Gamma_ab=-i*rho.
        near(projected.arm_mean_square_v2_m2[0], 1.0, 3e-13);
        near(projected.arm_mean_square_v2_m2[1], 1.0, 3e-13);
        complex_near(projected.mutual_v2_m2, C64::new(0.0, -rho), 4e-13);
        near(projected.visibility.unwrap(), rho, 4e-13);
        for phi in [-core::f64::consts::FRAC_PI_2, 0.0, core::f64::consts::FRAC_PI_2] {
            let measured = bank
                .recombined_electric_mean_square_v2_m2(a, b, [ey(); 2], phi, limits, || false)
                .unwrap();
            near(measured, 2.0 - 2.0 * rho * det::sin(phi), 2e-12);
        }
    }
}

#[test]
fn g1_direct_factor_recombination_retains_a_weak_signal_in_a_strong_dark_port() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let epsilon = 2.0_f64.powi(-30);
    let strong = vec![[C64::ZERO, C64::ONE]; grid().width];
    let weak: Vec<_> = (0..grid().width)
        .map(|x| {
            [C64::ZERO, phase(-core::f64::consts::TAU * x as f64 / grid().width as f64).scale(epsilon)]
        })
        .collect();
    let factors = [prepare(&stack, &strong), prepare(&stack, &weak)];
    let bank =
        AngularSpectrumCoherence::try_new(&factors[0], &factors, coherence_limits(), || false)
            .unwrap();
    let a = observation(Side::Transmitted, 0.0, 0.0, 0.0);
    let b = observation(
        Side::Transmitted,
        0.0,
        grid().width as f64 * grid().spacing_m[0] / 2.0,
        0.0,
    );
    let limits = query_limits(&bank);
    let analyzers = [ey(), ey().map(|z| -z)];
    let moments = bank.mutual_coherence(a, b, limits, || false).unwrap().project(analyzers).unwrap();
    let rounded_subtraction = moments.arm_mean_square_v2_m2[0]
        + moments.arm_mean_square_v2_m2[1]
        + 2.0 * moments.mutual_v2_m2.re;
    assert_eq!(rounded_subtraction, 0.0);
    let direct = bank
        .recombined_electric_mean_square_v2_m2(a, b, analyzers, 0.0, limits, || false)
        .unwrap();
    // The strong constant factor cancels; the independent weak factor doubles.
    // Its cycle mean is (2*epsilon)^2/2 = 2^-59, despite rounded moments losing it.
    let expected = 2.0 * epsilon * epsilon;
    near(direct, expected, 3e-12 * expected);
    assert!(direct > 0.0);
}

#[test]
fn g1_lossy_airy_film_preserves_reflection_transmission_and_standing_wave_correlation() {
    let medium = material(E, M);
    let film = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let thickness = 0.4;
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&film),
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let factors = [prepare(&stack, &vec![[C64::ZERO, C64::ONE]; grid().width])];
    let bank =
        AngularSpectrumCoherence::try_new(&factors[0], &factors, coherence_limits(), || false)
            .unwrap();
    // Independent Debye law at omega*tau=1: epsilon/E = 2.5 - 0.5i.
    // Sum the usual scalar Airy series, without querying production responses.
    let index = C64::new(2.5, -0.5).sqrt();
    let k = 1e9 * det::sqrt(E * M);
    let r01 = (C64::ONE - index) / (C64::ONE + index);
    let r12 = -r01;
    let travel = phase(-k * index.re * thickness).scale(det::exp(k * index.im * thickness));
    let denominator = C64::ONE + r01 * r12 * travel * travel;
    let reflected = (r01 + r12 * travel * travel) / denominator;
    let transmitted = (C64::ONE + r01) * (C64::ONE + r12) * travel / denominator;
    let incoming = phase(k * 0.17);
    let outgoing = reflected * phase(-k * 0.17);
    let exit = transmitted * phase(-k * 0.23);
    let b = observation(Side::Transmitted, 0.23, 0.0, 0.0);
    for (side, field) in [
        (Side::Incident, incoming),
        (Side::Reflected, outgoing),
        (Side::TotalIncidentSide, incoming + outgoing),
    ] {
        let a = observation(side, 0.17, 0.0, 0.0);
        let limits = query_limits(&bank);
        let statistics = bank.mutual_coherence(a, b, limits, || false).unwrap().project([ey(); 2]).unwrap();
        near(statistics.arm_mean_square_v2_m2[0], 0.5 * field.norm_sq(), 3e-12);
        near(statistics.arm_mean_square_v2_m2[1], 0.5 * exit.norm_sq(), 3e-12);
        complex_near(statistics.mutual_v2_m2, (field * exit.conj()).scale(0.5), 3e-12);
        let phi = 0.37;
        let direct = bank
            .recombined_electric_mean_square_v2_m2(a, b, [ey(); 2], phi, limits, || false)
            .unwrap();
        near(direct, 0.5 * (field + phase(phi) * exit).norm_sq(), 5e-12);
    }
}

#[test]
fn g1_evanescent_transmission_retains_two_point_phase_and_decay() {
    let medium = material(1.0, 1.0);
    let exit = material(0.01, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&exit), 1.0, &[])
            .unwrap();
    let input: Vec<_> = (0..grid().width)
        .map(|x| [C64::ZERO, phase(-2.0 * core::f64::consts::TAU * x as f64 / grid().width as f64)])
        .collect();
    let factors = [prepare(&stack, &input)];
    let bank =
        AngularSpectrumCoherence::try_new(&factors[0], &factors, coherence_limits(), || false)
            .unwrap();
    let kx = 2.0 * core::f64::consts::TAU / (grid().width as f64 * grid().spacing_m[0]);
    let q0 = det::sqrt(1.0 - kx * kx);
    let alpha = det::sqrt(kx * kx - 0.01);
    let t = C64::from_re(2.0 * q0) / C64::new(q0, -alpha);
    let a = observation(Side::Transmitted, 0.2, 1.3, 0.0);
    let b = observation(Side::Transmitted, 2.1, 9.7, 0.0);
    let result = bank.mutual_coherence(a, b, query_limits(&bank), || false).unwrap();
    let expected = phase(-kx * (a.position_m[0] - b.position_m[0]))
        .scale(0.5 * t.norm_sq() * det::exp(-alpha * (a.distance_from_face_m + b.distance_from_face_m)));
    complex_near(result.electric_cross_v2_m2()[1][1], expected, 3e-12);
    assert!(result.electric_at_b_v2_m2()[1][1].re < result.electric_at_a_v2_m2()[1][1].re);
}

#[test]
fn g3_two_dimensional_tensor_matches_fft_planes_reverses_and_survives_unitary_factor_rotation() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let g = AngularSpectrumGrid { width: 4, height: 4, spacing_m: [8.0, 9.0] };
    let j = C64::new(0.0, 1.0);
    let mut first = Vec::new();
    let mut second = Vec::new();
    let mut rotated_first = Vec::new();
    let mut rotated_second = Vec::new();
    for y in 0..g.height {
        for x in 0..g.width {
            let p = phase(-core::f64::consts::TAU * (x as f64 / g.width as f64 + y as f64 / g.height as f64));
            let q = phase(-core::f64::consts::TAU * (x as f64 / g.width as f64 - y as f64 / g.height as f64));
            let a = [C64::ONE + j * p.scale(0.4), C64::from_re(0.2) + q.scale(0.3)];
            let b = [p.scale(0.2), -j.scale(0.5) + q.scale(0.1)];
            first.push(a);
            second.push(b);
            rotated_first.push(std::array::from_fn(|c| (a[c] + j * b[c]).scale(core::f64::consts::FRAC_1_SQRT_2)));
            rotated_second.push(std::array::from_fn(|c| (j * a[c] + b[c]).scale(core::f64::consts::FRAC_1_SQRT_2)));
        }
    }
    let prepare_grid = |input: &[[C64; 2]]| PlanarAngularSpectrum::prepare(&stack, g, input, limits(), || false).unwrap();
    let factors = [prepare_grid(&first), prepare_grid(&second)];
    let rotated = [prepare_grid(&rotated_first), prepare_grid(&rotated_second)];
    let bank = AngularSpectrumCoherence::try_new(&factors[0], &factors, coherence_limits(), || false).unwrap();
    let rotated_bank = AngularSpectrumCoherence::try_new(&rotated[0], &rotated, coherence_limits(), || false).unwrap();
    let a = observation(Side::Transmitted, 0.3, g.spacing_m[0], 2.0 * g.spacing_m[1]);
    let b = observation(Side::Incident, 1.1, 3.0 * g.spacing_m[0], g.spacing_m[1]);
    let result = bank.mutual_coherence(a, b, query_limits(&bank), || false).unwrap();
    let reversed = bank.mutual_coherence(b, a, query_limits(&bank), || false).unwrap();
    let a_plane = bank.field(a.side, a.distance_from_face_m, || false).unwrap();
    let b_plane = bank.field(b.side, b.distance_from_face_m, || false).unwrap();
    tensor_near(result.electric_at_a_v2_m2(), a_plane.points()[2 * g.width + 1].electric_coherency_v2_m2, 3e-12);
    tensor_near(result.electric_at_b_v2_m2(), b_plane.points()[g.width + 3].electric_coherency_v2_m2, 3e-12);
    assert!(result.electric_at_a_v2_m2()[2][2].re > 1e-5);
    for i in 0..3 {
        for j in 0..3 {
            complex_near(result.electric_cross_v2_m2()[i][j], reversed.electric_cross_v2_m2()[j][i].conj(), 2e-13);
        }
    }
    same(&result, &rotated_bank.mutual_coherence(a, b, query_limits(&rotated_bank), || false).unwrap());
    let translated = AngularSpectrumObservation {
        position_m: [a.position_m[0] + g.width as f64 * g.spacing_m[0], a.position_m[1] - g.height as f64 * g.spacing_m[1]],
        ..a
    };
    same(&result, &bank.mutual_coherence(translated, b, query_limits(&bank), || false).unwrap());
}

#[test]
fn g1_complex_analyzers_use_hermitian_projection_and_independent_polarizations_do_not_cross() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let j = C64::new(0.0, 1.0);
    let circular = [prepare(&stack, &vec![[C64::ONE, j]; grid().width])];
    let bank = AngularSpectrumCoherence::try_new(&circular[0], &circular, coherence_limits(), || false).unwrap();
    let a = observation(Side::Incident, 0.0, 0.0, 0.0);
    let result = bank.mutual_coherence(a, a, query_limits(&bank), || false).unwrap();
    let plus = [C64::ONE, j, C64::ZERO].map(|z| z.scale(core::f64::consts::FRAC_1_SQRT_2));
    let minus = [C64::ONE, -j, C64::ZERO].map(|z| z.scale(core::f64::consts::FRAC_1_SQRT_2));
    let passed = result.project([plus; 2]).unwrap();
    near(passed.arm_mean_square_v2_m2[0], 1.0, 4e-15);
    near(passed.visibility.unwrap(), 1.0, 4e-15);
    let rejected = result.project([plus, minus]).unwrap();
    assert_eq!(rejected.arm_mean_square_v2_m2[1], 0.0);
    assert_eq!(rejected.mutual_v2_m2, C64::ZERO);
    assert_eq!(rejected.visibility, Some(0.0));
    let independent = [
        prepare(&stack, &vec![[C64::ONE, C64::ZERO]; grid().width]),
        prepare(&stack, &vec![[C64::ZERO, C64::ONE]; grid().width]),
    ];
    let bank = AngularSpectrumCoherence::try_new(&independent[0], &independent, coherence_limits(), || false).unwrap();
    let projected = bank.mutual_coherence(a, a, query_limits(&bank), || false).unwrap()
        .project([[C64::ONE, C64::ZERO, C64::ZERO], ey()]).unwrap();
    assert_eq!(projected.arm_mean_square_v2_m2, [0.5; 2]);
    assert_eq!(projected.mutual_v2_m2, C64::ZERO);
    assert_eq!(projected.visibility, Some(0.0));
}

#[test]
fn g0_two_point_queries_admit_work_geometry_darkness_and_cancellation_before_publication() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let factors = [prepare(&stack, &vec![[C64::ZERO, C64::ONE]; grid().width])];
    let bank = AngularSpectrumCoherence::try_new(&factors[0], &factors, coherence_limits(), || false).unwrap();
    let dark = AngularSpectrumCoherence::try_new(&factors[0], &[], coherence_limits(), || false).unwrap();
    let a = observation(Side::Transmitted, 0.0, 0.0, 0.0);
    let zero = dark.mutual_coherence(a, a, query_limits(&dark), || false).unwrap();
    assert_eq!(zero.electric_cross_v2_m2(), [[C64::ZERO; 3]; 3]);
    assert_eq!(zero.project([ey(); 2]).unwrap().visibility, None);
    assert_eq!(dark.correlation_mode_evaluations().unwrap(), 0);
    assert_eq!(dark.recombined_electric_mean_square_v2_m2(a, a, [ey(); 2], 0.0, query_limits(&dark), || false).unwrap(), 0.0);
    for invalid in [
        AngularSpectrumObservation { distance_from_face_m: -1.0, ..a },
        AngularSpectrumObservation { distance_from_face_m: f64::INFINITY, ..a },
        AngularSpectrumObservation { position_m: [f64::NAN, 0.0], ..a },
        AngularSpectrumObservation { position_m: [0.0, f64::INFINITY], ..a },
    ] {
        assert!(dark.mutual_coherence(invalid, a, query_limits(&dark), || false).is_err());
        assert!(dark.recombined_electric_mean_square_v2_m2(a, invalid, [ey(); 2], 0.0, query_limits(&dark), || false).is_err());
    }
    let invalid_analyzers = [[C64::new(f64::NAN, 0.0); 3], ey()];
    assert!(zero.project(invalid_analyzers).is_err());
    assert!(dark.recombined_electric_mean_square_v2_m2(a, a, invalid_analyzers, 0.0, query_limits(&dark), || false).is_err());
    assert!(dark.recombined_electric_mean_square_v2_m2(a, a, [ey(); 2], f64::NAN, query_limits(&dark), || false).is_err());
    for limits in [
        AngularSpectrumCorrelationLimits { maximum_mode_evaluations: bank.correlation_mode_evaluations().unwrap() - 1, ..query_limits(&bank) },
        AngularSpectrumCorrelationLimits { maximum_workspace_bytes: bank.correlation_workspace_bytes() - 1, ..query_limits(&bank) },
    ] {
        assert!(bank.mutual_coherence(a, a, limits, || false).is_err());
        assert!(bank.recombined_electric_mean_square_v2_m2(a, a, [ey(); 2], 0.0, limits, || false).is_err());
    }
    assert!(matches!(dark.mutual_coherence(a, a, query_limits(&dark), || true), Err(AngularSpectrumError::Cancelled)));
    let mut polls = 0;
    assert!(matches!(bank.recombined_electric_mean_square_v2_m2(a, a, [ey(); 2], 0.0, query_limits(&bank), || {
        polls += 1;
        polls > 3
    }), Err(AngularSpectrumError::Cancelled)));
    near(bank.recombined_electric_mean_square_v2_m2(a, a, [ey(); 2], 0.0, query_limits(&bank), || false).unwrap(), 2.0, 1e-14);
}
