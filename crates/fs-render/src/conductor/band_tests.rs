use super::*;
use crate::dielectric::{
    CauchyIor, DielectricCoatingStack, DielectricError, DielectricIor, DielectricIorSample,
    TabulatedIor,
};
use core::f64::consts::PI;

const BAND: (f64, f64) = (8000.0, 14_000.0);

fn source(tag: u8) -> ConductorSource {
    ConductorSource::try_new(
        ContentHash([tag; 32]),
        ConductorDataStatus::CallerAssertedMeasured,
    )
    .unwrap()
}
fn wavelengths(bounds: (f64, f64)) -> [f64; 9] {
    core::array::from_fn(|index| bounds.0 + (bounds.1 - bounds.0) * index as f64 / 8.0)
}
fn metal(bounds: (f64, f64), eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new_in_band(
        wavelengths(bounds)
            .map(|lambda| ConductorIorSample::try_new_in_band(lambda, eta, k).unwrap()),
        (bounds.1 - bounds.0) / 8.0,
        source(1),
    )
    .unwrap()
}
fn film(bounds: (f64, f64), eta: f64, k: f64, thickness: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        wavelengths(bounds)
            .map(|lambda| CoatingIorSample::try_new_in_band(lambda, eta, k).unwrap()),
        thickness,
        (bounds.1 - bounds.0) / 8.0,
        source(2),
    )
    .unwrap()
}
fn real(bounds: (f64, f64), index: f64) -> DielectricIor {
    TabulatedIor::try_new_in_band(
        wavelengths(bounds)
            .map(|lambda| DielectricIorSample::try_new_in_band(lambda, index).unwrap()),
        (bounds.1 - bounds.0) / 8.0,
    )
    .unwrap()
    .into()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4.0e-11 * expected.abs().max(f64::MIN_POSITIVE),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn original_tables_admit_explicit_ir_support_without_changing_visible_identity_or_samples() {
    let samples = core::array::from_fn(|i| {
        ConductorIorSample::try_new_in_band(
            8000.0 + 750.0 * i as f64,
            1.0 + 0.125 * i as f64,
            if i == 8 { 1.0e-200 } else { 2.0 },
        )
        .unwrap()
    });
    let optics = ConductorOptics::try_new_in_band(samples, 750.0, source(4)).unwrap();
    assert_eq!(optics.wavelength_bounds_nm(), BAND);
    assert_eq!(optics.max_knot_span_nm(), 750.0);
    assert_eq!(optics.source(), source(4));
    assert_eq!(optics.samples(), samples);
    for sample in samples {
        let evaluated = optics.complex_ior(sample.wavelength_nm()).unwrap();
        assert_eq!(evaluated.eta().to_bits(), sample.eta().to_bits());
        assert_eq!(evaluated.k().to_bits(), sample.k().to_bits());
    }
    close(optics.complex_ior(8375.0).unwrap().eta(), 1.0625);
    for lambda in [380.0, 7999.0, 14_001.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            optics.complex_ior(lambda),
            Err(ConductorError::InvalidWavelength)
        );
    }
    assert!(ConductorIorSample::try_new(10_000.0, 1.5, 0.2).is_err());
    assert!(CoatingIorSample::try_new(10_000.0, 1.5, 0.0).is_err());
    assert!(DielectricIorSample::try_new(10_000.0, 1.5).is_err());
    let visible = wavelengths((380.0, 780.0));
    let samples = visible.map(|lambda| ConductorIorSample::try_new(lambda, 0.8, 3.2).unwrap());
    let original = ConductorOptics::try_new(samples, source(4)).unwrap();
    let explicit = ConductorOptics::try_new_in_band(samples, 100.0, source(4)).unwrap();
    assert_eq!(original, explicit);
    assert_eq!(original.content_identity(), explicit.content_identity());
    let samples = visible.map(|lambda| CoatingIorSample::try_new(lambda, 1.4, 0.2).unwrap());
    let original = CoherentAbsorbingCoating::try_new(samples, 100.0e-9, source(4)).unwrap();
    let explicit =
        CoherentAbsorbingCoating::try_new_in_band(samples, 100.0e-9, 100.0, source(4)).unwrap();
    assert_eq!(original, explicit);
    let identity = |film| {
        CoherentCoatingStack::try_new(&[CoatingLayer::Absorbing(film)])
            .unwrap()
            .content_identity()
    };
    assert_eq!(identity(&original), identity(&explicit));
    let changed_gap =
        CoherentAbsorbingCoating::try_new_in_band(samples, 100.0e-9, 150.0, source(4)).unwrap();
    assert_ne!(identity(&original), identity(&changed_gap));
    let samples =
        wavelengths(BAND).map(|lambda| DielectricIorSample::try_new_in_band(lambda, 1.5).unwrap());
    let real = TabulatedIor::try_new_in_band(samples, 750.0).unwrap();
    assert_eq!(real.wavelength_bounds_nm(), BAND);
    assert_eq!(real.max_knot_span_nm(), 750.0);
    assert_eq!(real.eval(10_000.0), Ok(1.5));
    assert_eq!(real.group_index(10_000.0), Ok(1.5));
    assert_eq!(real.eval(780.0), Err(DielectricError::InvalidWavelength));
    assert_eq!(
        CauchyIor::try_constant(1.5).unwrap().eval(10_000.0),
        Err(DielectricError::InvalidWavelength)
    );
}

#[test]
fn ir_matched_films_keep_actual_nm_phase_and_analytic_depth_absorption() {
    let eta = 1.5;
    let k = 0.3;
    let substrate = metal(BAND, eta, k);
    let thickness = [2.0e-6, 3.0e-6];
    let films = [
        film(BAND, eta, k, thickness[0]),
        film(BAND, eta, k, thickness[1]),
    ];
    let layers = [
        CoatingLayer::Absorbing(&films[0]),
        CoatingLayer::Absorbing(&films[1]),
    ];
    let stack = CoherentCoatingStack::try_new(&layers).unwrap();
    let reflection = ((1.0 - eta) * (1.0 - eta) + k * k) / ((1.0 + eta) * (1.0 + eta) + k * k);
    let entrance = 4.0 * eta / ((1.0 + eta) * (1.0 + eta) + k * k);
    for lambda in [8000.0, 10_000.0, 14_000.0] {
        close(substrate.fresnel(lambda, 1.0, 1.0).unwrap(), reflection);
        close(
            substrate.fresnel_mueller(lambda, 1.0, 1.0).unwrap().0[0][0],
            reflection,
        );
        let alpha = 4.0 * PI * k / (lambda * 1.0e-9);
        let expected = [
            entrance * -(-alpha * thickness[0]).exp_m1(),
            entrance * (-alpha * thickness[0]).exp() * -(-alpha * thickness[1]).exp_m1(),
        ];
        let work = stack
            .layer_power_balance(&substrate, lambda, 1.0, 1.0)
            .unwrap();
        for pol in 0..2 {
            close(work.total.reflected[pol], reflection);
            close(work.total.absorbed[pol], entrance);
            close(work.layer_absorbed()[0][pol], expected[0]);
            close(work.layer_absorbed()[1][pol], expected[1]);
            close(
                work.substrate_absorbed[pol],
                entrance * (-alpha * (thickness[0] + thickness[1])).exp(),
            );
        }
        let attached = substrate.with_absorbing_coating(films[0]);
        close(attached.fresnel(lambda, 1.0, 1.0).unwrap(), reflection);
        close(
            attached.fresnel_mueller(lambda, 1.0, 1.0).unwrap().0[0][0],
            reflection,
        );
        close(
            attached
                .layer_power_balance(lambda, 1.0, 1.0)
                .unwrap()
                .layer_absorbed()[0][0],
            expected[0],
        );
    }
    // Independent normal-incidence quarter-wave impedance transformation.
    let lambda = 10_000.0;
    let film_eta = 1.5;
    let transparent = CoherentDielectricCoating::try_new(
        real(BAND, film_eta),
        lambda * 1.0e-9 / (4.0 * film_eta),
        source(3),
    )
    .unwrap();
    let (eta, k) = (0.8, 3.2);
    let coated = metal(BAND, eta, k).with_coating(transparent);
    let effective_real = film_eta * film_eta * eta / (eta * eta + k * k);
    let effective_imag = -film_eta * film_eta * k / (eta * eta + k * k);
    let expected = ((1.0 - effective_real).powi(2) + effective_imag.powi(2))
        / ((1.0 + effective_real).powi(2) + effective_imag.powi(2));
    close(coated.fresnel(lambda, 1.0, 1.0).unwrap(), expected);
    close(
        coated.fresnel_mueller(lambda, 1.0, 1.0).unwrap().0[0][0],
        expected,
    );
    assert_eq!(transparent.wavelength_bounds_nm(), BAND);
}

#[test]
fn ir_transmission_preserves_weak_loss_tir_and_reverse_authored_partition() {
    let k = 1.0e-20;
    let depth = 1.0e-6;
    let weak = film(BAND, 1.0, k, depth);
    let layers = [CoatingLayer::Absorbing(&weak)];
    let stack = DielectricCoatingStack::try_mixed_in_band(&layers, BAND).unwrap();
    let work = stack
        .layer_power_balance(10_000.0, false, 1.0, 1.0, 1.0)
        .unwrap();
    close(work.absorbed[0], 4.0 * PI * k * depth / 10_000.0e-9);
    assert!(work.absorbed[0] > 0.0);
    let outgoing = work.fresnel.reflection.0[0][0] + work.fresnel.transmission.0[0][0];
    // Native coherent amplitudes incur ordinary flux roundoff near unity.
    // The physical loss is much smaller: adding it cannot change R+T, yet the
    // positive field-work owner above still returns its analytic value.
    assert!((outgoing - 1.0).abs() <= 4.0 * f64::EPSILON);
    assert_eq!(
        outgoing + 0.5 * (work.absorbed[0] + work.absorbed[1]),
        outgoing
    );

    let films = [film(BAND, 1.2, 0.08, 0.7e-6), film(BAND, 2.0, 0.2, 1.1e-6)];
    let layers = [
        CoatingLayer::Absorbing(&films[0]),
        CoatingLayer::Absorbing(&films[1]),
    ];
    let reversed_layers = [layers[1], layers[0]];
    let stack = DielectricCoatingStack::try_mixed_in_band(&layers, BAND).unwrap();
    let reversed = DielectricCoatingStack::try_mixed_in_band(&reversed_layers, BAND).unwrap();
    for lambda in [8000.0, 11_000.0, 14_000.0] {
        let inside = real(BAND, 1.5).eval(lambda).unwrap();
        let outside = real(BAND, 1.0).eval(lambda).unwrap();
        let tir = stack
            .layer_power_balance(lambda, false, 0.5, inside, outside)
            .unwrap();
        assert!(tir.fresnel.scalar.total_internal_reflection);
        assert_eq!(tir.fresnel.transmission.0[0][0], 0.0);
        assert!(tir.absorbed.iter().all(|fraction| *fraction > 0.0));
        close(
            tir.fresnel.reflection.0[0][0] + 0.5 * (tir.absorbed[0] + tir.absorbed[1]),
            1.0,
        );
        let back = stack
            .layer_power_balance(lambda, true, 0.7, outside, inside)
            .unwrap();
        let reordered = reversed
            .layer_power_balance(lambda, false, 0.7, outside, inside)
            .unwrap();
        for pol in 0..2 {
            close(back.absorbed[pol], reordered.absorbed[pol]);
            close(
                back.layer_absorbed()[0][pol],
                reordered.layer_absorbed()[1][pol],
            );
            close(
                back.layer_absorbed()[1][pol],
                reordered.layer_absorbed()[0][pol],
            );
        }
        close(
            back.fresnel.reflection.0[0][0]
                + back.fresnel.transmission.0[0][0]
                + 0.5 * (back.absorbed[0] + back.absorbed[1]),
            1.0,
        );
    }
}

#[test]
fn every_consumed_law_keeps_support_even_through_zero_depth_grazing_and_bare_shortcuts() {
    let broad = metal((8000.0, 16_000.0), 1.5, 0.3);
    let narrow = (9000.0, 13_000.0);
    let zero_film = film(narrow, 1.5, 0.3, 0.0);
    let transparent =
        CoherentDielectricCoating::try_new(real(narrow, 1.5), 0.0, source(3)).unwrap();
    let legacy =
        CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.5).unwrap(), 0.0, source(3))
            .unwrap();
    let periodic = CoherentPeriodicCoating::try_new(
        &[CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.5).unwrap(), 0.0, 0.0).unwrap()],
        1,
        source(3),
    )
    .unwrap();
    for attached in [
        broad.with_absorbing_coating(zero_film),
        broad.with_coating(transparent),
        broad.with_coating(legacy),
        broad.with_periodic_coating(periodic),
    ] {
        for cosine in [0.0, 1.0] {
            assert_eq!(
                attached.fresnel(8500.0, 1.0, cosine),
                Err(ConductorError::InvalidWavelength)
            );
            assert_eq!(
                attached.fresnel_mueller(8500.0, 1.0, cosine),
                Err(ConductorError::InvalidWavelength)
            );
            assert_eq!(
                attached.power_balance(8500.0, 1.0, cosine),
                Err(ConductorError::InvalidWavelength)
            );
            assert_eq!(
                attached.layer_power_balance(8500.0, 1.0, cosine),
                Err(ConductorError::InvalidWavelength)
            );
        }
    }
    let layers = [CoatingLayer::Absorbing(&zero_film)];
    let stack = CoherentCoatingStack::try_new(&layers).unwrap();
    assert_eq!(
        stack.fresnel(&broad, 8500.0, 1.0, 0.0),
        Err(ConductorError::InvalidWavelength)
    );
    assert_eq!(
        stack.power_balance(&broad, 8500.0, 1.0, 0.0),
        Err(ConductorError::InvalidWavelength)
    );
    assert_eq!(
        stack.layer_power_balance(&broad, 8500.0, 1.0, 0.0),
        Err(ConductorError::InvalidWavelength)
    );
    assert!(matches!(
        DielectricCoatingStack::try_mixed_in_band(&layers, BAND),
        Err(DielectricError::InvalidWavelength)
    ));
    assert!(matches!(
        DielectricCoatingStack::try_mixed(&layers),
        Err(DielectricError::InvalidWavelength)
    ));
    let bare = DielectricCoatingStack::try_mixed_in_band(&[], BAND).unwrap();
    close(
        bare.fresnel(10_000.0, false, 1.0, 1.0, 1.5)
            .unwrap()
            .reflectance,
        0.04,
    );
    for lambda in [780.0, 7999.0, 14_001.0] {
        assert_eq!(
            bare.fresnel(lambda, false, 1.0, 1.0, 1.5),
            Err(DielectricError::InvalidWavelength)
        );
        assert_eq!(
            bare.layer_power_balance(lambda, false, 1.0, 1.0, 1.5),
            Err(DielectricError::InvalidWavelength)
        );
    }
    let point = DielectricCoatingStack::try_mixed_in_band(&layers, (10_000.0, 10_000.0)).unwrap();
    assert_eq!(
        point
            .layer_power_balance(10_000.0, false, 1.0, 1.0, 1.5)
            .unwrap()
            .absorbed,
        [0.0; 2]
    );
    assert!(point.fresnel(10_001.0, false, 1.0, 1.0, 1.5).is_err());
    for bounds in [(0.0, 10_000.0), (14_000.0, 8000.0), (8000.0, f64::INFINITY)] {
        assert!(DielectricCoatingStack::try_mixed_in_band(&[], bounds).is_err());
    }
    let samples = wavelengths(BAND)
        .map(|lambda| ConductorIorSample::try_new_in_band(lambda, 1.5, 0.3).unwrap());
    for gap in [0.0, -1.0, 749.0, f64::NAN, f64::INFINITY] {
        assert!(ConductorOptics::try_new_in_band(samples, gap, source(1)).is_err());
    }
}
