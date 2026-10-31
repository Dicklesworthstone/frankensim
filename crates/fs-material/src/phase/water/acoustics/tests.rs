use super::*;
use crate::fluid_acoustics::FrozenHomogeneousMixture;
use fs_qty::Temperature;

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-20),
        "actual={actual:e}, expected={expected:e}"
    );
}

// An independent finite pressure/volume experiment: impose a nearby saturation
// temperature and recover its phase inventory from trapezoidal reversible work.
// This uses only original u, v and p values, no analytic response derivatives.
fn nearby_adiabatic_state(pair: If97SaturationPair, x: f64, delta_t: f64) -> (f64, f64, f64) {
    let source =
        If97SaturationPair::at_temperature(Temperature::new(pair.temperature().value() + delta_t))
            .unwrap();
    let u0 = (1.0 - x) * pair.liquid().specific_internal_energy().value()
        + x * pair.vapor().specific_internal_energy().value();
    let v0 = (1.0 - x) * pair.liquid().specific_volume().value()
        + x * pair.vapor().specific_volume().value();
    let p_bar = 0.5 * (source.pressure().value() + pair.pressure().value());
    let l = source.liquid();
    let v = source.vapor();
    let du = v.specific_internal_energy().value() - l.specific_internal_energy().value();
    let dv = v.specific_volume().value() - l.specific_volume().value();
    let next_x =
        (u0 - l.specific_internal_energy().value() - p_bar * (l.specific_volume().value() - v0))
            / (du + p_bar * dv);
    assert!((0.0..1.0).contains(&next_x));
    (
        source.pressure().value(),
        l.specific_volume().value() + next_x * dv,
        next_x,
    )
}

#[test]
fn analytic_tangent_matches_independent_reversible_work_compression() {
    for temperature in [300.0, 400.0, 600.0] {
        for quality in [0.01, 0.2, 0.8] {
            let pair = If97SaturationPair::at_temperature(Temperature::new(temperature)).unwrap();
            let response = HomogeneousEquilibriumWater::try_new(pair, quality).unwrap();
            let lower = nearby_adiabatic_state(pair, quality, -0.002);
            let upper = nearby_adiabatic_state(pair, quality, 0.002);
            let v = 1.0 / response.state().density_kg_m3();
            let bulk = -v * (upper.0 - lower.0) / (upper.1 - lower.1);
            near(response.state().adiabatic_bulk_modulus_pa(), bulk, 2e-7);
            near(
                response.tangent().vapor_mass_fraction_derivative_per_pa,
                (upper.2 - lower.2) / (upper.0 - lower.0),
                2e-7,
            );
        }
    }
}

#[test]
fn true_phase_equilibration_changes_compliance_and_preserves_original_source() {
    let pair = If97SaturationPair::at_temperature(Temperature::new(400.0)).unwrap();
    let quality = 1e-3;
    let response = HomogeneousEquilibriumWater::try_new(pair, quality).unwrap();
    let frozen = FrozenHomogeneousMixture::from_mass_fraction(
        pair.liquid().acoustic_state().unwrap(),
        pair.vapor().acoustic_state().unwrap(),
        quality,
    )
    .unwrap();
    assert_eq!(*response.source(), pair);
    assert_eq!(response.vapor_mass_fraction(), quality);
    near(
        response.vapor_volume_fraction(),
        frozen.void_fraction(),
        1e-13,
    );
    near(
        response.state().density_kg_m3(),
        frozen.state().density_kg_m3(),
        1e-13,
    );
    assert!(response.state().sound_speed_m_s() < frozen.state().sound_speed_m_s());
    assert!(response.state().sound_speed_m_s() > 0.0);
}

#[test]
fn original_source_seams_are_visible_without_energy_or_entropy_adjustment() {
    let pair = If97SaturationPair::at_temperature(Temperature::new(500.0)).unwrap();
    let response = HomogeneousEquilibriumWater::try_new(pair, 0.2).unwrap();
    let tangent = response.tangent();
    assert_eq!(
        tangent.source_gibbs_gap_j_kg,
        pair.gibbs_energy_gap().value()
    );
    let energy = tangent.specific_internal_energy_derivative_m3_kg;
    let work = pair.pressure().value() * tangent.specific_volume_derivative_m3_kg_pa;
    assert!((energy + work).abs() < 1e-13 * energy.abs().max(work.abs()));
    // Inexact coexistence matching is the only entropy defect; do not silently
    // replace latent entropy by h_fg/T to manufacture a zero diagnostic.
    let entropy_from_seam = -tangent.vapor_mass_fraction_derivative_per_pa
        * pair.gibbs_energy_gap().value()
        / pair.temperature().value();
    near(
        tangent.source_entropy_derivative_j_kg_k_pa,
        entropy_from_seam,
        5e-7,
    );
    assert!(tangent.source_relative_clapeyron_defect.is_finite());
}

#[test]
fn exact_phase_boundaries_and_invalid_qualities_refuse_ambiguous_acoustics() {
    let pair = If97SaturationPair::at_temperature(Temperature::new(400.0)).unwrap();
    for quality in [0.0, 1.0, -0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(HomogeneousEquilibriumWater::try_new(pair, quality).is_err());
    }
    for quality in [1e-6, 1.0 - 1e-6] {
        assert!(HomogeneousEquilibriumWater::try_new(pair, quality).is_ok());
    }
}
