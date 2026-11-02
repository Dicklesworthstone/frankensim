use super::*;
use crate::conductor::{
    CauchyCoatingLayer, CoatingIorSample, CoherentAbsorbingCoating, CoherentDielectricCoating,
    CoherentPeriodicCoating, ConductorDataStatus, ConductorIorSample, ConductorSource,
};
use crate::dielectric::CauchyIor;
use fs_blake3::ContentHash;

fn source() -> ConductorSource {
    ConductorSource::try_new(ContentHash([91; 32]), ConductorDataStatus::Representative).unwrap()
}
fn substrate(eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        source(),
    )
    .unwrap()
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
fn g1_normal_flux_matches_independent_absorptance_and_survives_unit_reflectance() {
    for (incident, eta, k) in [(1.0, 0.2, 3.0), (1.33, 2.7, 4.1), (1.0, 1e-20, 3.0)] {
        let optics = substrate(eta, k);
        let balance = optics.power_balance(550.0, incident, 1.0).unwrap();
        // Independent normal-incidence flux form, without deriving absorption
        // from the production reflection output.
        let expected = 4.0 * incident * eta / ((eta + incident).powi(2) + k * k);
        for lane in 0..2 {
            close(balance.absorbed[lane], expected, 2e-14 * expected);
            close(balance.absorbed[lane] + balance.reflected[lane], 1.0, 2e-14);
        }
        assert_eq!(balance.absorbed[0], balance.absorbed[1]);
        if eta == 1e-20 {
            assert_eq!(balance.reflected, [1.0; 2]);
            assert!(balance.absorbed[0] > 0.0);
        }
    }
}

#[test]
fn g1_coherent_layers_match_independent_airy_power_and_original_polarizations() {
    let bare = substrate(0.17, 3.1);
    let transparent =
        CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.45).unwrap(), 90e-9, source())
            .unwrap();
    let absorbing = film(2.1, 0.2, 45.0);
    // Independent complex Airy recursion at normal incidence, n0=1,
    // substrate=.17+3.1i, lambda=550 nm. It does not propagate F/G vectors.
    for (optics, expected) in [
        (bare.with_coating(transparent), 0.101_844_608_534_183_31),
        (
            bare.with_absorbing_coating(absorbing),
            0.787_885_410_571_495_1,
        ),
    ] {
        let power = optics.power_balance(550.0, 1.0, 1.0).unwrap();
        close(power.absorbed[0], expected, 5e-13);
        close(power.absorbed[1], expected, 5e-13);
        close(
            power.reflected[0],
            optics.fresnel(550.0, 1.0, 1.0).unwrap(),
            5e-13,
        );
    }
    let layers = [
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.1).unwrap(), 0.05, 45e-9).unwrap(),
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.45).unwrap(), 0.0, 90e-9).unwrap(),
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.5).unwrap(), 0.3, 20e-9).unwrap(),
    ];
    let optics =
        bare.with_periodic_coating(CoherentPeriodicCoating::try_new(&layers, 1, source()).unwrap());
    let power = optics.power_balance(550.0, 1.0, 0.83).unwrap();
    close(power.absorbed[0], 0.226_264_094_373_882_8, 5e-13);
    close(power.absorbed[1], 0.327_327_775_028_380_6, 5e-13);
    let matrix = optics.fresnel_mueller(550.0, 1.0, 0.83).unwrap();
    close(power.reflected[0], matrix.0[0][0] + matrix.0[0][1], 5e-13);
    close(power.reflected[1], matrix.0[0][0] - matrix.0[0][1], 5e-13);
    for pol in 0..2 {
        close(power.reflected[pol] + power.absorbed[pol], 1.0, 5e-13);
    }
}

#[test]
fn g3_zero_thickness_matched_layers_and_explicit_stack_preserve_opaque_balance() {
    let bare = substrate(0.17, 3.1);
    let zero = film(2.1, 0.4, 0.0);
    let matched = film(0.17, 3.1, 20.0);
    let incident =
        CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.0).unwrap(), 30e-9, source())
            .unwrap();
    for cosine in [0.0, 0.17, 0.83, 1.0] {
        let expected = bare.power_balance(550.0, 1.0, cosine).unwrap();
        for optics in [
            bare.with_absorbing_coating(zero),
            bare.with_absorbing_coating(matched),
            bare.with_coating(incident),
        ] {
            assert_eq!(optics.power_balance(550.0, 1.0, cosine).unwrap(), expected);
        }
        let layers = [
            CoatingLayer::Absorbing(&zero),
            CoatingLayer::Dielectric(&incident),
            CoatingLayer::Absorbing(&matched),
        ];
        assert_eq!(
            CoherentCoatingStack::try_new(&layers)
                .unwrap()
                .power_balance(&bare, 550.0, 1.0, cosine)
                .unwrap(),
            expected
        );
        if cosine == 0.0 {
            assert_eq!(
                expected,
                ConductorPowerBalance {
                    reflected: [1.0; 2],
                    absorbed: [0.0; 2]
                }
            );
        }
    }
}

#[test]
fn g0_invalid_surface_queries_refuse_before_any_power_account() {
    let optics = substrate(0.17, 3.1);
    for (wavelength, index, cosine) in [
        (379.0, 1.0, 0.5),
        (781.0, 1.0, 0.5),
        (550.0, 0.0, 0.5),
        (550.0, f64::INFINITY, 0.5),
        (550.0, 1.0, -0.1),
        (550.0, 1.0, f64::NAN),
    ] {
        assert!(optics.power_balance(wavelength, index, cosine).is_err());
    }
    assert!(ConductorIorSample::try_new(550.0, 1.0, 0.0).is_err());
}
