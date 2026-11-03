use super::*;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use std::collections::BTreeMap;

fn source(maximum_concentration: f64) -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(AQUEOUS_LANGMUIR_LAW.into()),
        law_version: 1,
        parameters: BTreeMap::from([
            (
                CAPACITY_PARAMETER.into(),
                LawParameter {
                    value: 2.0,
                    dims: AdsorbentLoading::DIMS,
                },
            ),
            (
                HALF_SATURATION_PARAMETER.into(),
                LawParameter {
                    value: 1.0,
                    dims: AmountConcentration::DIMS,
                },
            ),
        ]),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                CONCENTRATION_AXIS,
                concentration_spec(),
                0.0,
                maximum_concentration,
            )
            .with_quantity(TEMPERATURE_AXIS, temperature_spec(), 290.0, 310.0)
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0, 12.0),
        sources: vec![hash_domain(
            "authored-langmuir-fixture",
            b"qmax=2mol/kg;half=1mol/m3",
        )],
        provenance: Provenance {
            source: "Authored arithmetic fixture; not a measured gold/carbon calibration".into(),
            license: "project test fixture".into(),
            artifact: None,
        },
    }
}
fn model(maximum: f64) -> AqueousLangmuir {
    AqueousLangmuir::from_card(ElementId::new("Au").unwrap(), source(maximum)).unwrap()
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 11.0)
        .unwrap()
}
fn inventory(dissolved: f64, adsorbed: f64) -> AdsorptionInventory {
    AdsorptionInventory {
        dissolved: Amount::new(dissolved),
        adsorbed: Amount::new(adsorbed),
    }
}
fn contact(
    model: &AqueousLangmuir,
    initial: AdsorptionInventory,
    mass: f64,
) -> Result<AdsorptionContact> {
    model.equilibrate(
        initial,
        Volume::new(2.0),
        Mass::new(mass),
        Temperature::new(298.15),
        &conditions(),
    )
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1e-15),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_finite_bath_loading_and_preloaded_carbon_desorption() {
    let card = model(10.0);
    let initial = inventory(8.0, 0.0);
    let answer = contact(&card, initial, 3.0).unwrap();
    // c=2, V=2 and C=6: dissolved=4, adsorbed=6*2/(1+2)=4 mol.
    close(answer.final_inventory.dissolved.value(), 4.0);
    close(answer.final_inventory.adsorbed.value(), 4.0);
    close(answer.equilibrium_concentration.value(), 2.0);
    close(answer.loading.unwrap().value(), 4.0 / 3.0);
    close(answer.vacant_capacity.value(), 2.0);
    close(answer.transferred_to_adsorbent.value(), 4.0);
    close(answer.conservation_residual.value(), 0.0);
    assert_eq!(answer.initial, initial);
    assert_eq!(answer.model(), &card);
    assert_eq!(answer.model().source(), &source(10.0));
    let unchanged = contact(&card, answer.final_inventory, 3.0).unwrap();
    close(unchanged.transferred_to_adsorbent.value(), 0.0);

    // Initially loaded adsorbent releases one mole into initially clean water.
    // c=0.5, dissolved=1, adsorbed=6*0.5/1.5=2, total=3 mol.
    let wash = contact(&card, inventory(0.0, 3.0), 3.0).unwrap();
    close(wash.final_inventory.dissolved.value(), 1.0);
    close(wash.final_inventory.adsorbed.value(), 2.0);
    close(wash.transferred_to_adsorbent.value(), -1.0);
    close(wash.vacant_capacity.value(), 4.0);
}

#[test]
fn g1_g3_independent_bisection_oracle_scaling_and_small_vacant_capacity() {
    let card = model(1e20);
    for total in [1e-6, 1.0, 1e6] {
        for capacity in [0.001, 1.0, 100.0] {
            let answer = contact(&card, inventory(total, 0.0), capacity / 2.0).unwrap();
            // Independent monotone material-balance solve, without a quadratic.
            let (mut lo, mut hi) = (0.0, total);
            for _ in 0..200 {
                let x = 0.5 * (lo + hi);
                if x + capacity * (x / 2.0) / (1.0 + x / 2.0) > total {
                    hi = x;
                } else {
                    lo = x;
                }
            }
            close(answer.final_inventory.dissolved.value(), 0.5 * (lo + hi));
            close(
                answer.final_inventory.dissolved.value() + answer.final_inventory.adsorbed.value(),
                total,
            );
        }
    }
    let base = contact(&card, inventory(8.0, 0.0), 3.0).unwrap();
    for factor in [1e-6, 1e6] {
        let answer = card
            .equilibrate(
                inventory(8.0 * factor, 0.0),
                Volume::new(2.0 * factor),
                Mass::new(3.0 * factor),
                Temperature::new(298.15),
                &conditions(),
            )
            .unwrap();
        close(
            answer.equilibrium_concentration.value(),
            base.equilibrium_concentration.value(),
        );
        close(
            answer.final_inventory.adsorbed.value(),
            base.final_inventory.adsorbed.value() * factor,
        );
    }
    let nearly_full = contact(&card, inventory(1e18, 0.0), 3.0).unwrap();
    assert!(nearly_full.vacant_capacity.value() > 0.0);
    close(nearly_full.vacant_capacity.value(), 1.2e-17);
    assert!(nearly_full.final_inventory.adsorbed.value() <= 6.0);
}

#[test]
fn g3_source_support_applies_to_solved_equilibrium_and_cannot_be_overridden() {
    // Initial c=4 is outside support, but the actual equilibrium c=2 is inside.
    let card = model(3.0);
    assert!(contact(&card, inventory(8.0, 0.0), 3.0).is_ok());
    assert!(
        matches!(contact(&model(1.9), inventory(8.0, 0.0), 3.0), Err(AqueousAdsorptionError::OutsideValidity { axis, .. }) if axis == CONCENTRATION_AXIS)
    );
    for query in [
        QueryPoint::new(),
        QueryPoint::new().with("pH", 11.0).unwrap(),
        conditions().with("pH", 9.0).unwrap(),
        conditions()
            .with_quantity(CONCENTRATION_AXIS, concentration_spec(), 1.0)
            .unwrap(),
        conditions()
            .with_quantity(TEMPERATURE_AXIS, temperature_spec(), 298.15)
            .unwrap(),
    ] {
        assert!(
            card.equilibrate(
                inventory(8.0, 0.0),
                Volume::new(2.0),
                Mass::new(3.0),
                Temperature::new(298.15),
                &query
            )
            .is_err()
        );
    }
    assert!(
        card.equilibrate(
            inventory(8.0, 0.0),
            Volume::new(2.0),
            Mass::new(3.0),
            Temperature::new(320.0),
            &conditions()
        )
        .is_err()
    );
}

#[test]
fn g0_empty_batches_zero_adsorbent_and_invalid_stocks() {
    let card = model(10.0);
    let empty = contact(&card, inventory(0.0, 0.0), 3.0).unwrap();
    assert_eq!(empty.final_inventory, inventory(0.0, 0.0));
    close(empty.vacant_capacity.value(), 6.0);
    let bypass = contact(&card, inventory(8.0, 0.0), 0.0).unwrap();
    assert_eq!(bypass.final_inventory, inventory(8.0, 0.0));
    assert!(bypass.loading.is_none());
    for (initial, mass) in [
        (inventory(0.0, 6.1), 3.0),
        (inventory(0.0, 1.0), 0.0),
        (inventory(-1.0, 0.0), 3.0),
        (inventory(1.0, -1.0), 3.0),
        (inventory(1.0, 0.0), -1.0),
    ] {
        assert!(contact(&card, initial, mass).is_err());
    }
    for volume in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            card.equilibrate(
                inventory(1.0, 0.0),
                Volume::new(volume),
                Mass::new(3.0),
                Temperature::new(298.15),
                &conditions()
            )
            .is_err()
        );
    }
}

#[test]
fn g0_card_admission_rejects_gas_units_hidden_state_and_unqualified_domains() {
    let id = ElementId::new("Au").unwrap();
    let mut bad = source(10.0);
    bad.parameters
        .get_mut(HALF_SATURATION_PARAMETER)
        .unwrap()
        .dims = fs_qty::Pressure::DIMS;
    assert!(AqueousLangmuir::from_card(id.clone(), bad).is_err());
    let mut bad = source(10.0);
    bad.sources.clear();
    assert!(AqueousLangmuir::from_card(id.clone(), bad).is_err());
    let mut bad = source(10.0);
    bad.state_schema_version = 1;
    assert!(AqueousLangmuir::from_card(id.clone(), bad).is_err());
    let mut bad = source(10.0);
    bad.law = LawId("gas-langmuir".into());
    assert!(AqueousLangmuir::from_card(id.clone(), bad).is_err());
    let mut bad = source(10.0);
    bad.validity = bad.validity.with_quantity(
        "elapsed-time",
        QuantitySpec::dimensional(fs_qty::Time::DIMS),
        0.0,
        100.0,
    );
    assert!(AqueousLangmuir::from_card(id.clone(), bad).is_err());
    let mut bad = source(10.0);
    bad.validity = bad.validity.with_quantity(
        TEMPERATURE_AXIS,
        QuantitySpec::dimensional(Temperature::DIMS),
        290.0,
        310.0,
    );
    assert!(AqueousLangmuir::from_card(id.clone(), bad).is_err());
    let mut bad = source(10.0);
    bad.validity = bad.validity.with("unknown-untyped-condition", 0.0, 1.0);
    assert!(AqueousLangmuir::from_card(id.clone(), bad).is_err());
    let mut bad = source(10.0);
    bad.parameters.get_mut(CAPACITY_PARAMETER).unwrap().value = 0.0;
    assert!(AqueousLangmuir::from_card(id, bad).is_err());
}

#[test]
fn g0_numerical_loss_refuses_without_a_fabricated_zero_inventory() {
    let card = model(1e20);
    assert!(contact(&card, inventory(f64::INFINITY, 0.0), 3.0).is_err());
    assert!(contact(&card, inventory(1.0, 0.0), f64::MAX).is_err());
    // Scaled positive half-saturation is unrepresentable at this dynamic range.
    // Refuse rather than return complete removal and erase a positive tail.
    let mut tiny = source(1e20);
    tiny.parameters
        .get_mut(HALF_SATURATION_PARAMETER)
        .unwrap()
        .value = 1e-300;
    let tiny = AqueousLangmuir::from_card(ElementId::new("Au").unwrap(), tiny).unwrap();
    assert!(
        tiny.equilibrate(
            inventory(1e300, 0.0),
            Volume::new(1.0),
            Mass::new(1.0),
            Temperature::new(298.15),
            &conditions()
        )
        .is_err()
    );
}
