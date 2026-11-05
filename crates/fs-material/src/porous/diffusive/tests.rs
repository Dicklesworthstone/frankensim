use super::*;
use crate::fluid_acoustics::FluidAcousticState;
use crate::fluid_transport::FluidTransportState;

fn fluid(gamma: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::from_adiabatic_bulk_modulus(1.204, 1.4 * 101_325.0).unwrap(),
        FluidTransportState::try_new(1.204, 1005.0, 1.81e-5 * 1005.0 / 0.71, 1.81e-5).unwrap(),
        gamma,
    )
    .unwrap()
}
fn law(card: super::super::PorousCard) -> PorousFluidLaw {
    PorousFluidLaw::try_new(card, 1.81e-5).unwrap()
}
fn band(material: &PorousFluidLaw) -> [f64; 2] {
    let (lo, hi) = material.reference_card().valid_hz;
    [core::f64::consts::TAU * lo, core::f64::consts::TAU * hi]
}

#[test]
fn measured_jcal_cards_converge_against_original_complex_law_off_grid() {
    for card in [
        super::super::PorousCard::melamine_uf_nguyen2024(),
        super::super::PorousCard::basotect_foam02(),
        super::super::PorousCard::pinta_fleece_foam02(),
    ] {
        let material = law(card);
        let state = fluid(1.4);
        let model = PorousDiffusive::new(&material, &state, 32, band(&material), 1e-6).unwrap();
        assert_eq!(model.material().reference_card(), material.reference_card());
        assert_eq!(model.fluid(), state);
        assert_eq!(model.series_terms().len(), 32);
        assert_eq!(model.shunt_terms().len(), 33);
        assert!(
            model
                .series_terms()
                .iter()
                .chain(model.shunt_terms())
                .all(|t| t.gain_per_s > 0.0 && t.rate_per_s > 0.0)
        );
        let [lo, hi] = model.band_rad_s();
        for i in 0..137 {
            // Deliberately distinct from the producer's logarithmic sample grid.
            let omega = lo + (hi - lo) * ((i as f64 + 0.37) / 137.0);
            let exact = material.effective_fluid(&state, omega).unwrap();
            let (rho, bulk) = model.effective_density_bulk(omega).unwrap();
            assert!(
                (rho - exact.dynamic_density_kg_m3()).abs() / exact.dynamic_density_kg_m3().abs()
                    < 1e-6
            );
            assert!(
                (bulk - exact.dynamic_bulk_modulus_pa()).abs()
                    / exact.dynamic_bulk_modulus_pa().abs()
                    < 1e-6
            );
            // Real losses must converge independently of reversible storage.
            let s = C64::new(0.0, omega);
            let (z, y) = model.normalized_response(omega).unwrap();
            let exact_z = s * exact
                .dynamic_density_kg_m3()
                .scale(1.0 / model.density_infinite_kg_m3());
            let exact_y = s / exact
                .dynamic_bulk_modulus_pa()
                .scale(1.0 / model.bulk_infinite_pa());
            assert!((z.re / exact_z.re - 1.0).abs() < 1e-6);
            assert!((y.re / exact_y.re - 1.0).abs() < 1e-6);
        }
    }
}

#[test]
fn champoux_allard_is_the_original_lafarge_degeneracy() {
    let mut card = super::super::PorousCard::melamine_uf_nguyen2024();
    card.thermal_permeability_m2 = None;
    let a = law(card.clone());
    card.thermal_permeability_m2 =
        Some(card.porosity * card.thermal_length_m * card.thermal_length_m / 8.0);
    let b = law(card);
    let ma = PorousDiffusive::new(&a, &fluid(1.4), 24, band(&a), 1e-8).unwrap();
    let mb = PorousDiffusive::new(&b, &fluid(1.4), 24, band(&b), 1e-8).unwrap();
    assert_eq!(ma.series_terms(), mb.series_terms());
    assert_eq!(ma.shunt_terms(), mb.shunt_terms());
    assert_eq!(
        ma.normalized_response(5000.0).unwrap(),
        mb.normalized_response(5000.0).unwrap()
    );
}

#[test]
fn isolated_thermal_pole_and_weak_thermal_loss_are_not_discarded() {
    let material = law(super::super::PorousCard::melamine_uf_nguyen2024());
    let model =
        PorousDiffusive::new(&material, &fluid(1.0 + 1e-9), 24, band(&material), 1e-8).unwrap();
    let isolated = model.shunt_terms()[0];
    assert!(
        model.shunt_terms()[1..]
            .iter()
            .all(|term| isolated.rate_per_s < term.rate_per_s)
    );
    let omega = isolated
        .rate_per_s
        .clamp(model.band_rad_s()[0], model.band_rad_s()[1]);
    let full = model.normalized_response(omega).unwrap().1.re;
    let mut without_pole = model.clone();
    without_pole.shunt.remove(0);
    let missing = without_pole.normalized_response(omega).unwrap().1.re;
    assert!(full > 0.0);
    assert!(
        (full - missing) / full > 0.5,
        "isolated thermal pole is physical loss"
    );
    assert!(model.sampled_max_real_loss_error() < 1e-8);
}

#[test]
fn cp_equals_cv_removes_only_thermal_memory() {
    let material = law(super::super::PorousCard::melamine_uf_nguyen2024());
    let model = PorousDiffusive::new(&material, &fluid(1.0), 16, band(&material), 1e-8).unwrap();
    assert!(model.shunt_terms().is_empty());
    assert_eq!(model.series_terms().len(), 16);
    assert!(model.series_drag_per_s() > 0.0);
    let (z, y) = model.normalized_response(5000.0).unwrap();
    assert!(z.re > model.series_drag_per_s());
    assert_eq!(y, C64::new(0.0, 5000.0));
    let (_, bulk) = model.effective_density_bulk(5000.0).unwrap();
    assert_eq!(bulk, C64::from_re(model.bulk_infinite_pa()));
}

#[test]
fn original_permeability_and_actual_fluid_modulus_survive_substitution() {
    let material = law(super::super::PorousCard::melamine_uf_nguyen2024());
    let liquid = ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::from_adiabatic_bulk_modulus(1000.0, 2.2e9).unwrap(),
        FluidTransportState::try_new(1000.0, 4180.0, 0.6, 0.001).unwrap(),
        1.0,
    )
    .unwrap();
    let model = PorousDiffusive::new(&material, &liquid, 64, band(&material), 1e-6).unwrap();
    let card = material.reference_card();
    assert_eq!(model.fluid(), liquid);
    assert_eq!(model.bulk_infinite_pa(), 2.2e9 / card.porosity);
    let expected = 0.001 / material.viscous_permeability_m2() / model.density_infinite_kg_m3();
    assert!((model.series_drag_per_s() / expected - 1.0).abs() < 1e-14);
    assert_eq!(model.material().reference_card(), card);
    assert!(model.shunt_terms().is_empty());
}

#[test]
fn source_band_and_unresolved_physical_loss_are_refused() {
    let material = law(super::super::PorousCard::melamine_uf_nguyen2024());
    assert!(PorousDiffusive::new(&material, &fluid(1.4), 2, band(&material), 1e-3).is_err());
    let [lo, hi] = band(&material);
    assert!(PorousDiffusive::new(&material, &fluid(1.4), 32, [lo * 0.99, hi], 1e-4).is_err());
    let model = PorousDiffusive::new(&material, &fluid(1.4), 24, [lo, hi], 1e-8).unwrap();
    assert!(model.normalized_response(lo * 0.99).is_err());
    assert!(model.normalized_response(hi * 1.01).is_err());
    assert!(model.normalized_response(f64::NAN).is_err());
}
