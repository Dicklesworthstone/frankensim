//! G0/G3: finite phase complementarity, temperature direction and basis scaling.

use super::*;

fn law(enthalpy: f64) -> IdealSolubilityV1 {
    IdealSolubilityV1::new(
        SpeciesId::new("explicit-solvent").unwrap(),
        SpeciesId::new("neutral-solute").unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(2.0),
        enthalpy,
        [Temperature::new(250.0), Temperature::new(400.0)],
    )
    .unwrap()
}

#[test]
fn finite_inventory_respects_phase_complementarity() {
    let model = law(0.0);
    for (total, dissolved, solid) in [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (2.0, 2.0, 0.0),
        (5.0, 2.0, 3.0),
    ] {
        let state = model
            .equilibrate(
                Temperature::new(300.0),
                Volume::new(1.0),
                Amount::new(total),
            )
            .unwrap();
        assert_eq!(state.dissolved.value(), dissolved);
        assert_eq!(state.solid.value(), solid);
        assert_eq!(state.balance_residual.value(), 0.0);
        assert!(state.dissolved.value() <= state.capacity.value());
        assert_eq!(
            state.solid.value() * (state.capacity.value() - state.dissolved.value()),
            0.0
        );
    }
}

#[test]
fn solvent_loss_precipitates_and_dilution_redissolves_the_same_stock() {
    let model = law(0.0);
    let total = Amount::new(3.0);
    let initial = model
        .equilibrate(Temperature::new(300.0), Volume::new(2.0), total)
        .unwrap();
    let dry = model
        .equilibrate(
            Temperature::new(300.0),
            Volume::new(0.5),
            initial.total_solute,
        )
        .unwrap();
    let diluted = model
        .equilibrate(Temperature::new(300.0), Volume::new(2.0), dry.total_solute)
        .unwrap();
    assert_eq!(initial.dissolved.value(), 3.0);
    assert_eq!(dry.dissolved.value(), 1.0);
    assert_eq!(dry.solid.value(), 2.0);
    assert_eq!(diluted, initial);
}

#[test]
fn vant_hoff_sign_and_log_temperature_slope() {
    let model = law(24_000.0);
    let cold = model
        .saturation_concentration(Temperature::new(280.0))
        .unwrap()
        .value();
    let reference = model
        .saturation_concentration(Temperature::new(300.0))
        .unwrap()
        .value();
    let hot = model
        .saturation_concentration(Temperature::new(340.0))
        .unwrap()
        .value();
    assert!(cold < reference && reference < hot);
    let independently_expected = 2.0 * (-24_000.0 / R * (1.0 / 340.0 - 1.0 / 300.0)).exp();
    assert!((hot / independently_expected - 1.0).abs() < 2e-14);
    let slope = (hot / cold).ln() / (1.0 / 340.0 - 1.0 / 280.0);
    assert!((slope / (-24_000.0 / R) - 1.0).abs() < 2e-14);
    let exothermic = law(-24_000.0);
    assert!(
        exothermic
            .saturation_concentration(Temperature::new(340.0))
            .unwrap()
            .value()
            < reference
    );
}

#[test]
fn heating_changes_finite_partition_without_creating_solute() {
    let model = law(24_000.0);
    let cold = model
        .equilibrate(Temperature::new(280.0), Volume::new(1.0), Amount::new(4.0))
        .unwrap();
    let hot = model
        .equilibrate(Temperature::new(340.0), Volume::new(1.0), cold.total_solute)
        .unwrap();
    assert!(cold.solid.value() > 0.0);
    assert_eq!(hot.solid.value(), 0.0);
    assert_eq!(hot.dissolved.value(), 4.0);
    assert!(cold.balance_residual.value().abs() < 1e-14);
}

#[test]
fn reference_cached_partition_and_extensive_rescaling_agree() {
    let model = law(0.0);
    let temperature = Temperature::new(350.0);
    assert_eq!(
        model.saturation_concentration(temperature).unwrap(),
        model.reference_solubility()
    );
    for scale in [1e-18, 1.0, 1e18] {
        let volume = Volume::new(0.75 * scale);
        let total = Amount::new(5.0 * scale);
        let direct = model.equilibrate(temperature, volume, total).unwrap();
        assert_eq!(
            direct,
            partition_at_saturation(model.reference_solubility(), volume, total).unwrap()
        );
        assert!((direct.dissolved.value() / scale - 1.5).abs() < 1e-14);
        assert!((direct.solid.value() / scale - 3.5).abs() < 1e-14);
    }
}

#[test]
fn positive_dissolved_trace_survives_dominant_solid() {
    let result = partition_at_saturation(
        AmountConcentration::new(1e-30),
        Volume::new(1.0),
        Amount::new(1.0),
    )
    .unwrap();
    assert_eq!(result.dissolved.value(), 1e-30);
    assert_eq!(result.solid.value(), 1.0);
    assert_eq!(result.balance_residual.value(), 1e-30);
}

#[test]
fn temperature_domain_and_unrepresentable_capacity_refuse() {
    let model = law(24_000.0);
    for t in [0.0, 249.0, 401.0, f64::INFINITY, f64::NAN] {
        assert_eq!(
            model.saturation_concentration(Temperature::new(t)),
            Err(SolubilityErrorV1::OutsideTemperatureDomain)
        );
    }
    for (c, v, n) in [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, -1.0),
        (1e-250, 1e-250, 1.0),
        (1e250, 1e250, 1.0),
    ] {
        assert!(
            partition_at_saturation(AmountConcentration::new(c), Volume::new(v), Amount::new(n))
                .is_err()
        );
    }
    assert!(
        IdealSolubilityV1::new(
            model.solvent().clone(),
            model.solvent().clone(),
            Temperature::new(300.0),
            AmountConcentration::new(1.0),
            0.0,
            model.temperature_domain()
        )
        .is_err()
    );
}

#[test]
fn logarithmic_temperature_evaluation_preserves_representable_result() {
    let model = IdealSolubilityV1::new(
        SpeciesId::new("solvent").unwrap(),
        SpeciesId::new("solute").unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(1e-300),
        800.0 * 600.0 * R,
        [Temperature::new(300.0), Temperature::new(600.0)],
    )
    .unwrap();
    // exp(800) overflows, but 1e-300 * exp(800) is a finite concentration.
    let concentration = model
        .saturation_concentration(Temperature::new(600.0))
        .unwrap()
        .value();
    assert!(concentration.is_finite() && concentration > 0.0);
    assert!((concentration.ln() - (1e-300_f64.ln() + 800.0)).abs() < 2e-13);
}
