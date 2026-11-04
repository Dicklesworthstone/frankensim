//! General integer salts under the original source and acid/base contracts.

use super::*;

const FERROUS: &str = "Fe[2+](aq)";
const FERROUS_ARSENATE: &str = "authored-Fe3(AsO4)2(s)";

fn ferrous_model(kdiss: f64) -> SaltPrecipitation {
    let mut fixture = Fixture::new(kdiss);
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * T;
    fixture.records.push(Record::new(
        FERROUS,
        [0, 0, 0, 1, 0, 0],
        2,
        0.0,
        "source-phase-aqueous",
    ));
    fixture.records.push(Record::new(
        FERROUS_ARSENATE,
        [0, 8, 2, 3, 0, 0],
        0,
        -6.0 * rt * KA.ln() + rt * kdiss.ln(),
        "source-phase-solid",
    ));
    let reaction = fixture.reaction(
        "authored-ferrous-arsenate-precipitation",
        &[(FERROUS, -3), (FORMS[3], -2), (FERROUS_ARSENATE, 1)],
    );
    SaltPrecipitation::from_sources(
        fixture.aqueous(false),
        reaction.clone(),
        id(FERROUS),
        id(FERROUS_ARSENATE),
        FAMILY,
        &fixture.sources(&reaction),
    )
    .unwrap()
}

#[test]
fn three_to_two_salt_closes_both_limiting_orientations_and_extensive_scaling() {
    // At h=Ka each of the four acid forms has fraction 1/4. A solid formula
    // contains 3 Fe(II) and 2 arsenate components; the independently prescribed
    // equilibrium has Q=a_Fe^3*a_AsO4^2. The spectator amount fixes its exact
    // electroneutral solution at pH 7, not a frozen-pH approximation in the model.
    for (metal, family, solid) in [(0.003_f64, 0.004_f64, 0.002), (0.009, 0.002, 0.001)] {
        let kdiss = metal.powi(3) * (family / 4.0).powi(2);
        let model = ferrous_model(kdiss);
        assert_eq!(model.cation_stoichiometry(), 3);
        assert_eq!(model.anion_stoichiometry(), 2);
        for scale in [1e-6, 1.0, 1e6] {
            let result = model
                .speciate(
                    Amount::new(scale * (metal + 3.0 * solid)),
                    &[total(FAMILY, scale * (family + 2.0 * solid))],
                    &counterions(scale * (2.0 * metal - 1.5 * family)),
                    Amount::new(scale),
                    options(),
                    never_cancelled,
                )
                .unwrap();
            assert_eq!(result.status, PrecipitationStatus::SolidPresent);
            near(result.dissolved_cation.value() / scale, metal, 4e-9, 0.0);
            near(result.dissolved_family.value() / scale, family, 4e-9, 0.0);
            near(result.solid_amount.value() / scale, solid, 4e-9, 0.0);
            near(result.aqueous.ph.value(), 7.0, 0.0, 3e-9);
            near(
                result.dissolved_cation.value() + 3.0 * result.solid_amount.value(),
                result.cation_total.value(),
                4e-14,
                0.0,
            );
            near(
                result.dissolved_family.value() + 2.0 * result.solid_amount.value(),
                result.family_total.value(),
                4e-14,
                0.0,
            );
            closed(&result);
        }
    }
}

#[test]
fn common_ion_addition_precipitates_more_and_dilution_can_exhaust_the_solid() {
    let model = ferrous_model(2.7e-14);
    let original = solve(&model, 0.009, 0.008, &[]);
    // Add 0.01 mol FeCl2 with its actual 2 monovalent counterions. All acid
    // forms and pH respond in the same solve; no concentration multiplier or
    // prescribed pH substitutes for the original 3:2 activity product.
    let added = solve(&model, 0.019, 0.008, &counterions(0.02));
    assert_eq!(added.status, PrecipitationStatus::SolidPresent);
    assert!(added.solid_amount.value() > original.solid_amount.value());
    assert!(added.dissolved_cation.value() > original.dissolved_cation.value());
    let original_anion = original.aqueous.families[0].forms.last().unwrap();
    let added_anion = added.aqueous.families[0].forms.last().unwrap();
    assert!(added_anion.activity.value() < original_anion.activity.value());
    near(
        added.dissolved_cation.value().powi(3) * added_anion.activity.value().powi(2),
        2.7e-14,
        2e-9,
        0.0,
    );
    assert!((added.aqueous.ph.value() - original.aqueous.ph.value()).abs() > 0.01);
    closed(&added);

    let diluted = model
        .speciate(
            added.cation_total,
            &[total(FAMILY, added.family_total.value())],
            &counterions(0.02),
            Amount::new(1000.0),
            options(),
            never_cancelled,
        )
        .unwrap();
    assert_eq!(diluted.status, PrecipitationStatus::SolidAbsent);
    assert_eq!(diluted.dissolved_cation, added.cation_total);
    assert_eq!(diluted.dissolved_family, added.family_total);
    assert_eq!(diluted.solid_amount.value(), 0.0);
    assert!(diluted.log_saturation_ratio.unwrap().value() < 0.0);
    closed(&diluted);
}

#[test]
fn integer_stoichiometry_preserves_dissolved_traces_independently_of_large_solid() {
    for (metal, family) in [(1e-40_f64, 1.0_f64), (1.0, 1e-40)] {
        let model = ferrous_model(metal.powi(3) * (family / 4.0).powi(2));
        let result = solve(
            &model,
            metal + 3.0,
            family + 2.0,
            &counterions(2.0 * metal - 1.5 * family),
        );
        assert_eq!(result.status, PrecipitationStatus::SolidPresent);
        near(result.dissolved_cation.value(), metal, 4e-9, 0.0);
        near(result.dissolved_family.value(), family, 4e-9, 0.0);
        assert_eq!(result.solid_amount.value(), 1.0);
        let trace_balance = if metal < family {
            result.cation_balance_residual
        } else {
            result.family_balance_residual
        };
        near(trace_balance.value(), 1e-40, 4e-9, 0.0);
        assert!(result.dissolved_limiting_bracket[0].value() > 0.0);
        closed(&result);
    }
}

#[test]
fn stoichiometric_excess_uses_original_products_and_never_a_spurious_trace_floor() {
    let metal = 5.546_775_611_681_588e-7;
    let family = 8.320_163_417_522_381e-7;
    assert_eq!(metal / 2.0, family / 3.0);
    assert!(family - 3.0 * (metal / 2.0) < 0.0);
    let (limit, excess, cation_limits) = partition_basis(metal, family, 2, 3).unwrap();
    assert!(!cation_limits);
    assert_eq!(limit, family / 3.0);
    assert!(excess > 0.0);
    near(excess, 3.529_303_946_892_918e-23, 1e-15, 0.0);

    // Exact complete ratios need no invented epsilon excess, even when the
    // products used to compare them would overflow without the binary scale.
    assert_eq!(partition_basis(6.0, 4.0, 3, 2).unwrap(), (2.0, 0.0, true));
    let (limit, excess, cation_limits) = partition_basis(1e308, 1e308, 3, 2).unwrap();
    assert!(cation_limits);
    assert!(limit.is_finite() && limit > 0.0);
    near(excess, 1e308 / 3.0, 2e-15, 0.0);
    assert!(partition_basis(f64::from_bits(1), 1.0, 3, 2).is_err());
    // The limiting capacity can be positive while its nonlimiting excess
    // falls below representation. Refuse instead of declaring an exact ratio.
    let tiny = f64::from_bits(1);
    assert!(partition_basis(4.0 * tiny, 3.0 * tiny, 3, 2).is_err());

    let model = ferrous_model(2.7e-14);
    for (metal, family) in [(0.0, 0.0), (0.0, 0.1), (0.1, 0.0)] {
        let result = solve(
            &model,
            metal,
            family,
            &counterions(2.0 * metal - 1.5 * family),
        );
        assert_eq!(result.status, PrecipitationStatus::EmptyComponent);
        assert_eq!(result.dissolved_cation.value(), metal);
        assert_eq!(result.dissolved_family.value(), family);
        assert_eq!(result.solid_amount.value(), 0.0);
        closed(&result);
    }
}
