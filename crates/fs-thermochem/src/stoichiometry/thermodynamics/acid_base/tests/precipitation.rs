//! Analytical charge/solubility fixtures; no measured Fe constants are claimed.
use super::*;
use crate::stoichiometry::{
    minerals::MineralReaction, thermodynamics::acid_base::precipitation::*,
};

const METAL: &str = "Fe[3+](aq)";
const SOLID: &str = "FeO(OH)(s)";

fn model_with(
    fixture: &Fixture,
    aqueous: IdealAcidBaseSystem,
    k: f64,
    solid_phase: &str,
) -> Result<HydrolysisPrecipitation> {
    let reaction = MineralReaction::FerricHydrolysisGoethite.build().unwrap();
    let metal_point = point("source-phase-aqueous", true, fixture.reference);
    let solid_point = point(solid_phase, false, fixture.reference);
    let metal = pack(METAL, 0.0, &metal_point, None);
    let gibbs = if k == 1.0 {
        0.0
    } else {
        -crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE * k.ln()
    };
    let solid = pack(SOLID, gibbs, &solid_point, None);
    let sources = reaction
        .terms()
        .iter()
        .map(|term| {
            let (pack, at) = match term.species.as_str() {
                METAL => (&metal, &metal_point),
                SOLID => (&solid, &solid_point),
                other => {
                    let i = fixture.index(other);
                    (&fixture.packs[i], &fixture.points[i])
                }
            };
            FormationGibbsSource {
                species: &term.species,
                pack,
                at,
                reference_axis: fixture.reference,
            }
        })
        .collect::<Vec<_>>();
    HydrolysisPrecipitation::from_sources(aqueous, reaction.clone(), id(METAL), id(SOLID), &sources)
}

fn model(k: f64, buffered: bool) -> HydrolysisPrecipitation {
    let fixture = Fixture::new(1e-14, 1e-4, [1e-4, 1e-8]);
    let families = if buffered {
        vec![fixture.mono()]
    } else {
        vec![]
    };
    model_with(&fixture, fixture.system(families), k, "source-phase-solid").unwrap()
}

fn numerics() -> PrecipitationOptions {
    PrecipitationOptions {
        aqueous: options(),
        ..PrecipitationOptions::default()
    }
}

fn solve(
    model: &HydrolysisPrecipitation,
    total: f64,
    families: &[(String, Amount)],
    spectators: &[SpectatorIon],
) -> PrecipitationSpeciation {
    model
        .speciate(
            Amount::new(total),
            families,
            spectators,
            Amount::new(1.0),
            numerics(),
            || false,
        )
        .unwrap()
}

#[test]
fn coupled_charge_and_mineral_partition_recover_an_independent_known_state() {
    // Choose h=1e-4, Kw=1e-14, K_precip=1e-6: m=h^3/K=1e-6.
    // Electroneutrality then requires Cl=h+3*m-Kw/h, independently of the code.
    let model = model(1e-6, false);
    let chloride = 1e-4 + 3e-6 - 1e-10;
    let result = solve(&model, 2e-5, &[], &[ion("Cl[-](aq)", -1, chloride)]);
    assert_eq!(result.status, PrecipitationStatus::SolidPresent);
    near(result.aqueous.ph.value(), 4.0, 0.0, 2e-9);
    near(result.dissolved_amount.value(), 1e-6, 2e-9, 0.0);
    near(result.solid_amount.value(), 1.9e-5, 2e-9, 0.0);
    assert!(result.log_saturation_ratio.unwrap().value().abs() <= 1e-9);
    assert!(
        result.log_saturation_ratio.unwrap().value().abs()
            + result.log_saturation_numerical_spread.unwrap().value()
            <= 1e-9
    );
    assert!(result.component_balance_residual.value().abs() < 1e-20);
    assert!(result.aqueous.charge_balance_residual.value().abs() < 1e-15);
    assert!(result.dissolved_bracket[0].value() <= result.dissolved_amount.value());
    assert!(result.dissolved_bracket[1].value() >= result.dissolved_amount.value());
    assert!(result.aqueous_solves >= result.iterations + 2);
    assert!(result.aqueous_iterations > result.aqueous_solves);
    assert_eq!(result.model().cation().as_str(), METAL);
    assert_eq!(result.model().solid().as_str(), SOLID);
    assert_eq!(result.model().reaction().terms(), model.reaction().terms());
    assert_eq!(result.model().standard().molar(), model.standard().molar());
}

#[test]
fn buffering_is_solved_with_the_mineral_at_one_common_proton_activity() {
    // At h=Ka=1e-4 a 1e-3 mol HA/A- pool splits equally. Positive spectator
    // charge is 5e-4+Kw/h-h-3m, giving the SAME known m=1e-6 and pH=4.
    let model = model(1e-6, true);
    let families = [total("monoprotic", 1e-3)];
    let result = solve(
        &model,
        2e-5,
        &families,
        &[ion("Na+(aq)", 1, 5e-4 + 1e-10 - 1e-4 - 3e-6)],
    );
    assert_eq!(result.status, PrecipitationStatus::SolidPresent);
    near(result.aqueous.ph.value(), 4.0, 0.0, 2e-9);
    near(result.dissolved_amount.value(), 1e-6, 2e-9, 0.0);
    for form in &result.aqueous.families[0].forms {
        near(form.amount.value(), 5e-4, 2e-9, 0.0);
    }
    assert!(result.aqueous.max_log_mass_action_residual.value() < 1e-10);
}

#[test]
fn absent_solid_obeys_undersaturation_and_empty_component_keeps_its_distinct_status() {
    let model = model(1e-6, false);
    let result = solve(&model, 1e-10, &[], &[ion("Cl[-](aq)", -1, 1e-4)]);
    assert_eq!(result.status, PrecipitationStatus::SolidAbsent);
    assert_eq!(result.dissolved_amount.value(), 1e-10);
    assert_eq!(result.solid_amount.value(), 0.0);
    assert!(result.log_saturation_ratio.unwrap().value() < 0.0);
    assert_eq!(result.iterations, 0);
    let empty = solve(&model, 0.0, &[], &[]);
    assert_eq!(empty.status, PrecipitationStatus::EmptyComponent);
    assert_eq!(empty.log_saturation_ratio, None);
    assert_eq!(empty.dissolved_amount.value(), 0.0);
    assert_eq!(empty.solid_amount.value(), 0.0);
    near(empty.aqueous.ph.value(), 7.0, 0.0, 2e-11);
}

#[test]
fn trace_dissolved_metal_survives_when_extent_subtraction_would_erase_it() {
    let model = model(1.0, false);
    let result = solve(&model, 1e-3, &[], &[]);
    assert_eq!(result.status, PrecipitationStatus::SolidPresent);
    // Charge: h+3m=Kw/h, hence h is almost 1e-7 and m almost 1e-21.
    near(result.dissolved_amount.value(), 1e-21, 2e-9, 0.0);
    assert_eq!(result.solid_amount.value(), result.component_total.value());
    assert!(result.dissolved_amount.value() > 0.0);
    assert!(result.component_balance_residual.value() > 0.0);
    assert!(result.log_saturation_ratio.unwrap().value().abs() < 1e-9);
}

#[test]
fn extensive_scaling_and_spectator_order_preserve_the_partition_physics() {
    let model = model(1e-6, false);
    let ions = [
        ion("Cl[-](aq)", -1, 1e-4 + 3e-6 - 1e-10 + 2e-4),
        ion("Na+(aq)", 1, 2e-4),
    ];
    let first = solve(&model, 2e-5, &[], &ions);
    let reordered = solve(&model, 2e-5, &[], &[ions[1].clone(), ions[0].clone()]);
    assert_eq!(first.dissolved_amount, reordered.dissolved_amount);
    assert_eq!(first.aqueous.ph, reordered.aqueous.ph);
    let scaled_ions = ions.map(|mut ion| {
        ion.amount = Amount::new(7.0 * ion.amount.value());
        ion
    });
    let scaled = model
        .speciate(
            Amount::new(7.0 * 2e-5),
            &[],
            &scaled_ions,
            Amount::new(7.0),
            numerics(),
            || false,
        )
        .unwrap();
    near(
        scaled.aqueous.ph.value(),
        first.aqueous.ph.value(),
        0.0,
        2e-9,
    );
    near(
        scaled.dissolved_amount.value(),
        7.0 * first.dissolved_amount.value(),
        3e-9,
        0.0,
    );
    near(
        scaled.solid_amount.value(),
        7.0 * first.solid_amount.value(),
        3e-9,
        0.0,
    );
}

#[test]
fn original_shared_sources_and_pure_solid_phase_are_required() {
    let fixture = Fixture::new(1e-14, 1e-4, [1e-4, 1e-8]);
    let other = Fixture::with_reference(
        1e-14,
        1e-4,
        [1e-4, 1e-8],
        "source-elemental-reference-other-fixture",
    );
    assert!(model_with(&fixture, other.system(vec![]), 1e-6, "source-phase-solid").is_err());
    assert!(
        model_with(
            &fixture,
            fixture.system(vec![]),
            1e-6,
            "source-phase-aqueous"
        )
        .is_err()
    );
    assert!(
        model_with(
            &fixture,
            fixture.system(vec![]),
            1e-6,
            "source-phase-liquid"
        )
        .is_err()
    );
}

#[test]
fn bad_inputs_bounds_budget_and_cancellation_cannot_publish_partial_results() {
    let model = model(1e-6, false);
    let ions = [ion("Cl[-](aq)", -1, 1e-4 + 3e-6 - 1e-10)];
    let before = ions.clone();
    for total in [f64::NAN, -1.0, f64::INFINITY] {
        assert!(
            model
                .speciate(
                    Amount::new(total),
                    &[],
                    &ions,
                    Amount::new(1.0),
                    numerics(),
                    || false
                )
                .is_err()
        );
    }
    assert!(
        model
            .speciate(
                Amount::new(2e-5),
                &[],
                &[ion(METAL, 3, 1e-4)],
                Amount::new(1.0),
                numerics(),
                || false
            )
            .is_err()
    );
    assert!(matches!(
        model.speciate(
            Amount::new(2e-5),
            &[],
            &ions,
            Amount::new(1.0),
            PrecipitationOptions {
                minimum_dissolved_amount: Amount::new(2e-6),
                ..numerics()
            },
            || false
        ),
        Err(AcidBaseError::OutsideBracket)
    ));
    assert!(matches!(
        model.speciate(
            Amount::new(2e-5),
            &[],
            &ions,
            Amount::new(1.0),
            PrecipitationOptions {
                max_iterations: 1,
                ..numerics()
            },
            || false
        ),
        Err(AcidBaseError::Budget)
    ));
    let mut calls = 0;
    assert!(matches!(
        model.speciate(
            Amount::new(2e-5),
            &[],
            &ions,
            Amount::new(1.0),
            numerics(),
            || {
                calls += 1;
                calls >= 10
            }
        ),
        Err(AcidBaseError::Cancelled)
    ));
    assert_eq!(ions, before);
}

#[test]
fn loose_inner_ph_cannot_masquerade_as_tightly_resolved_mineral_equilibrium() {
    let model = model(1e-6, false);
    let loose = PrecipitationOptions {
        aqueous: SpeciationOptions {
            ph_tolerance: Dimensionless::new(0.5),
            absolute_charge_tolerance: Amount::new(1.0),
            relative_charge_tolerance: Dimensionless::new(0.01),
            ..options()
        },
        log_activity_tolerance: Dimensionless::new(1e-12),
        log_saturation_tolerance: Dimensionless::new(1e-12),
        ..numerics()
    };
    assert!(
        model
            .speciate(
                Amount::new(2e-5),
                &[],
                &[ion("Cl[-](aq)", -1, 1e-4 + 3e-6 - 1e-10)],
                Amount::new(1.0),
                loose,
                || false
            )
            .is_err()
    );
}
