use super::*;
use crate::conductor::{
    CoatingIorSample, CoherentAbsorbingCoating, ConductorDataStatus, ConductorSource,
};
use crate::dielectric::CauchyIor;
use fs_blake3::ContentHash;

fn source() -> ConductorSource {
    ConductorSource::try_new(ContentHash([92; 32]), ConductorDataStatus::Representative).unwrap()
}

fn film(eta: f64, k: f64, thickness_nm: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        thickness_nm * 1e-9,
        source(),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_weak_transmitting_film_absorption_survives_unit_transmittance() {
    let coating = film(1.5, 1e-20, 100.0);
    let layers = [CoatingLayer::Absorbing(&coating)];
    let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
    let balance = stack
        .layer_power_balance(550.0, false, 1.0, 1.5, 1.5)
        .unwrap();
    // Matched real indices: reflection is O(k²), so the first-order heat is
    // exactly Beer-Lambert 4*pi*k*d/lambda, long before R+T can resolve it.
    let expected = 4.0 * core::f64::consts::PI * 1e-20 * 100.0 / 550.0;
    for pol in 0..2 {
        assert!(balance.absorbed[pol] > 0.0);
        close(balance.absorbed[pol], expected, 3e-14 * expected);
        assert_eq!(balance.layer_absorbed()[0][pol], balance.absorbed[pol]);
    }
    assert!(balance.absorbed[0] < f64::EPSILON);
    assert_eq!(
        balance.fresnel,
        stack
            .fresnel_polarized(550.0, false, 1.0, 1.5, 1.5)
            .unwrap()
    );
    close(balance.fresnel.scalar.transmittance, 1.0, 2e-14);
}

#[test]
fn g1_tir_retains_positive_film_work_without_transmitted_heat() {
    for k in [0.02, 1e-20] {
        let coating = film(1.5, k, 80.0);
        let layers = [CoatingLayer::Absorbing(&coating)];
        let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
        let balance = stack
            .layer_power_balance(550.0, false, 0.5, 1.5, 1.0)
            .unwrap();
        assert!(balance.fresnel.scalar.total_internal_reflection);
        assert_eq!(balance.fresnel.scalar.transmittance, 0.0);
        assert!(balance.absorbed.iter().all(|loss| *loss > 0.0));
        assert_eq!(
            balance.fresnel,
            stack
                .fresnel_polarized(550.0, false, 0.5, 1.5, 1.0)
                .unwrap()
        );
        let [s, p] = [
            balance.fresnel.reflection.0[0][0] + balance.fresnel.reflection.0[0][1],
            balance.fresnel.reflection.0[0][0] - balance.fresnel.reflection.0[0][1],
        ];
        close(s + balance.absorbed[0], 1.0, 5e-13);
        close(p + balance.absorbed[1], 1.0, 5e-13);
    }
}

#[test]
fn g3_lossless_and_zero_films_have_exact_zero_heat_and_native_refusals() {
    let transparent =
        CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.0).unwrap(), 90e-9, source())
            .unwrap();
    let zero = film(2.1, 0.3, 0.0);
    let transparent_layers = [transparent];
    let mixed_layers = [
        CoatingLayer::Absorbing(&zero),
        CoatingLayer::Dielectric(&transparent),
    ];
    for stack in [
        DielectricCoatingStack::try_new(&transparent_layers).unwrap(),
        DielectricCoatingStack::try_mixed(&mixed_layers).unwrap(),
    ] {
        for cosine in [0.5, (1.0_f64 - 1.0 / 2.25).sqrt(), 1.0] {
            let balance = stack
                .layer_power_balance(550.0, false, cosine, 1.5, 1.0)
                .unwrap();
            assert_eq!(balance.absorbed, [0.0; 2]);
            assert!(
                balance
                    .layer_absorbed()
                    .iter()
                    .all(|loss| *loss == [0.0; 2])
            );
            assert_eq!(
                balance.fresnel,
                stack
                    .fresnel_polarized(550.0, false, cosine, 1.5, 1.0)
                    .unwrap()
            );
        }
        assert!(
            stack
                .layer_power_balance(550.0, false, 0.0, 1.5, 1.0)
                .is_err()
        );
        assert!(
            stack
                .layer_power_balance(379.0, false, 1.0, 1.5, 1.0)
                .is_err()
        );
    }
    let empty = DielectricCoatingStack::try_mixed(&[])
        .unwrap()
        .layer_power_balance(550.0, false, 0.0, 1.5, 1.5)
        .unwrap();
    assert_eq!(empty.absorbed, [0.0; 2]);
    assert!(empty.layer_absorbed().is_empty());
    assert_eq!(empty.fresnel.scalar.transmittance, 1.0);
}

#[test]
fn g1_opaque_transmitting_film_absorbs_before_the_exit_medium() {
    let opaque = film(2.1, 3.0, 100_000.0);
    let hidden = film(1.5, 0.02, 80.0);
    let layers = [
        CoatingLayer::Absorbing(&opaque),
        CoatingLayer::Absorbing(&hidden),
    ];
    let balance = DielectricCoatingStack::try_mixed(&layers)
        .unwrap()
        .layer_power_balance(550.0, false, 1.0, 1.0, 1.5)
        .unwrap();
    let expected = 4.0 * 2.1 / ((1.0_f64 + 2.1).powi(2) + 3.0_f64.powi(2));
    assert_eq!(balance.fresnel.scalar.transmittance, 0.0);
    assert_eq!(balance.layer_absorbed()[1], [0.0; 2]);
    for pol in 0..2 {
        close(balance.layer_absorbed()[0][pol], expected, 2e-13);
        assert_eq!(balance.layer_absorbed()[0][pol], balance.absorbed[pol]);
    }
}

#[test]
fn g1_oblique_and_reverse_film_work_matches_independent_boundary_fluxes() {
    let front = film(2.1, 0.2, 45.0);
    let spacer = film(1.45, 0.0, 90.0);
    let rear = film(2.5, 0.3, 20.0);
    let layers = [
        CoatingLayer::Absorbing(&front),
        CoatingLayer::Absorbing(&spacer),
        CoatingLayer::Absorbing(&rear),
    ];
    let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
    // Independent Python cmath transfer matrices followed by the normal
    // Poynting flux at each boundary; finite losses are well resolved here.
    // Expected layer entries retain original authored order on both sides.
    for (reverse, eta_i, eta_t, cosine, expected) in [
        (
            false,
            1.0,
            1.7,
            0.83,
            [
                [0.23191197137483371, 0.23969853908524097],
                [0.0, 0.0],
                [0.08103645270686266, 0.09595651235653158],
            ],
        ),
        (
            true,
            1.7,
            1.0,
            0.9,
            [
                [0.18061847959868738, 0.16238361072990382],
                [0.0, 0.0],
                [0.3471661865529634, 0.2304280795950746],
            ],
        ),
    ] {
        let balance = stack
            .layer_power_balance(550.0, reverse, cosine, eta_i, eta_t)
            .unwrap();
        assert_eq!(
            balance.fresnel,
            stack
                .fresnel_polarized(550.0, reverse, cosine, eta_i, eta_t)
                .unwrap()
        );
        assert_eq!(balance.layer_absorbed()[1], [0.0; 2]);
        for (actual, expected) in balance.layer_absorbed().iter().zip(expected) {
            for pol in 0..2 {
                close(actual[pol], expected[pol], 2e-12);
            }
        }
        for pol in 0..2 {
            let sum: f64 = balance.layer_absorbed().iter().map(|loss| loss[pol]).sum();
            assert_eq!(sum, balance.absorbed[pol]);
        }
        if reverse {
            let reversed_layers = [layers[2], layers[1], layers[0]];
            let reversed = DielectricCoatingStack::try_mixed(&reversed_layers)
                .unwrap()
                .layer_power_balance(550.0, false, cosine, eta_i, eta_t)
                .unwrap();
            assert_eq!(balance.fresnel, reversed.fresnel);
            assert_eq!(balance.absorbed, reversed.absorbed);
            for (original, reversed) in balance
                .layer_absorbed()
                .iter()
                .zip(reversed.layer_absorbed().iter().rev())
            {
                assert_eq!(original, reversed);
            }
        }
    }
}
