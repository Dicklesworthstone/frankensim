//! Independent finite-bath fixtures; all kinetic coefficients are synthetic.
use super::super::{
    AXES, CONCENTRATION_AXIS, DIFFUSION_RESISTANCE, DiffusionResistance, FILM_RESISTANCE, LAW,
    LOADING_AXIS, RADIUS_AXIS, SOLID_MOLAR_DENSITY, SURFACE_RESISTANCE, TransferResistance,
};
use super::*;
use crate::reaction_batch::{ELAPSED_TIME_AXIS, TEMPERATURE_AXIS, VOLUME_AXIS};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::QuantitySpec;
use fs_thermochem::stoichiometry::minerals::MineralReaction;
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix, StoichiometricMatrix};

fn card(resistances: [f64; 3]) -> ConstitutiveModelCard {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims) in AXES {
        let (lo, hi) = match axis {
            TEMPERATURE_AXIS => (298.15, 298.15),
            VOLUME_AXIS => (1.0, 1.0),
            RADIUS_AXIS => (1e-9, 1.0),
            _ => (0.0, 100.0),
        };
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi);
    }
    ConstitutiveModelCard {
        law: LawId(LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [
            (FILM_RESISTANCE, resistances[0], TransferResistance::DIMS),
            (SURFACE_RESISTANCE, resistances[1], TransferResistance::DIMS),
            (
                DIFFUSION_RESISTANCE,
                resistances[2],
                DiffusionResistance::DIMS,
            ),
            (SOLID_MOLAR_DENSITY, 1000.0, AmountConcentration::DIMS),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        validity,
        sources: vec![hash_domain(
            "shrinking-core-population-fixture",
            b"authored analytical coefficients; no measured ore kinetics",
        )],
        provenance: Provenance {
            source: "synthetic finite-bath population fixture, not plant calibration".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

fn model(
    recipe: MineralReaction,
    solid: &str,
    source: ConstitutiveModelCard,
) -> ShrinkingCoreReaction {
    ShrinkingCoreReaction::from_card(recipe.build().unwrap(), id(solid), id("H+(aq)"), source)
        .unwrap()
}

fn population(
    name: &str,
    model: ShrinkingCoreReaction,
    radius: f64,
    amount: f64,
) -> LeachPopulation {
    LeachPopulation::new(
        name,
        model,
        Length::new(radius),
        Amount::new(amount),
        QueryPoint::new(),
    )
    .unwrap()
}

fn inventory(populations: &[LeachPopulation], acid: f64) -> Vec<(SpeciesId, Amount)> {
    let mut amounts = BTreeMap::<SpeciesId, f64>::new();
    for p in populations {
        for term in p.model().reaction().terms() {
            amounts.entry(term.species.clone()).or_default();
        }
        *amounts.entry(p.model().reactants().0.clone()).or_default() += p.original_amount().value();
    }
    amounts.insert(id("H+(aq)"), acid);
    amounts
        .into_iter()
        .map(|(species, n)| (species, Amount::new(n)))
        .collect()
}

fn batch(populations: Vec<LeachPopulation>, acid: f64) -> ShrinkingCoreBatch {
    let initial = inventory(&populations, acid);
    ShrinkingCoreBatch::new(
        populations,
        id("H+(aq)"),
        &initial,
        Volume::new(1.0),
        Temperature::new(298.15),
    )
    .unwrap()
}

fn copper_bins() -> Vec<LeachPopulation> {
    let copper = model(
        MineralReaction::CupricOxideAcidLeach,
        "CuO(s)",
        card([1.5, 0.0, 0.0]),
    );
    vec![
        population("fine", copper.clone(), 0.001, 0.5),
        population("coarse", copper, 0.002, 0.5),
    ]
}

fn advance(batch: &mut ShrinkingCoreBatch, duration: f64) -> PopulationBatchAdvance {
    let settings = PopulationNumerics::default();
    let result = batch
        .advance(Time::new(duration), settings, || false)
        .unwrap();
    assert!(result.work <= settings.max_work);
    assert!(
        result.estimated_time_error.value()
            <= settings.absolute_time_tolerance.value()
                + settings.relative_time_tolerance * duration
    );
    result
}

fn amount(inventory: &[(SpeciesId, Amount)], species: &str) -> f64 {
    inventory
        .iter()
        .find(|(id, _)| id.as_str() == species)
        .unwrap()
        .1
        .value()
}
fn state<'a>(report: &'a PopulationBatchAdvance, name: &str) -> &'a PopulationState {
    report.populations.iter().find(|p| p.id == name).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}
fn snapshot(
    batch: &ShrinkingCoreBatch,
) -> (Vec<(SpeciesId, Amount)>, Vec<PopulationState>, f64, f64) {
    (
        batch.inventory().to_vec(),
        batch.populations().to_vec(),
        batch.elapsed().value(),
        batch.exposure().value(),
    )
}

#[test]
fn two_size_bins_share_acid_cross_completion_and_retain_exposure_during_hold() {
    // alpha=2, rho=1000, rf=1.5: F_fine=X and F_coarse=2X.
    // Before the first completion, C=3-1.5H and H=2(1-exp(-1.5t)).
    let mut bath = batch(copper_bins(), 3.0);
    let first_time = (2.0 / 3.0) * (4.0_f64 / 3.0).ln();
    let first = advance(&mut bath, first_time);
    close(state(&first, "fine").conversion.value(), 0.5, 4e-8);
    close(state(&first, "coarse").conversion.value(), 0.25, 4e-8);
    close(amount(bath.inventory(), "H+(aq)"), 2.25, 6e-8);
    close(first.cumulative_exposure.value(), 0.5, 4e-8);
    assert_eq!(first.cumulative.reactions.len(), 1);
    close(first.cumulative.reactions[0].extent.value(), 0.375, 4e-8);

    let first_completion = (2.0 / 3.0) * core::f64::consts::LN_2;
    let second_time = first_completion + 2.0 * (6.0_f64 / 5.0).ln();
    let second = advance(&mut bath, second_time - first_time);
    assert_eq!(state(&second, "fine").remaining.value(), 0.0);
    close(state(&second, "coarse").conversion.value(), 0.75, 5e-8);
    close(amount(bath.inventory(), "H+(aq)"), 1.25, 5e-8);
    close(second.cumulative_exposure.value(), 1.5, 5e-8);
    // The report is cumulative from original stocks, not this call's .5 mol.
    close(second.cumulative.reactions[0].extent.value(), 0.875, 5e-8);
    close(amount(&second.cumulative.consumed, "H+(aq)"), 1.75, 5e-8);
    close(amount(&second.cumulative.formed, "Cu[2+](aq)"), 0.875, 5e-8);
    close(second.previous_elapsed.value(), first_time, 1e-15);

    let completion_time = first_completion + 2.0 * (3.0_f64 / 2.0).ln();
    let completed = advance(&mut bath, completion_time + 0.25 - second_time);
    assert!(
        completed
            .populations
            .iter()
            .all(|p| p.remaining.value() == 0.0)
    );
    close(amount(bath.inventory(), "H+(aq)"), 1.0, 1e-14);
    close(completed.cumulative_exposure.value(), 2.25, 8e-8);
    let inventory_before_hold = bath.inventory().to_vec();
    let held = advance(&mut bath, 2.0);
    assert_eq!(held.exposure_increment.value(), 2.0);
    close(held.cumulative_exposure.value(), 4.25, 8e-8);
    assert_eq!(bath.inventory(), inventory_before_hold);
    for (p, point) in held.models().iter().zip(&held.final_points) {
        assert_eq!(point.axes()[RADIUS_AXIS], p.radius().value());
        assert_eq!(point.axes()[LOADING_AXIS], 1.0);
        close(
            point.axes()[ELAPSED_TIME_AXIS],
            completion_time + 2.25,
            1e-14,
        );
    }
}

#[test]
fn mixed_resistances_rebin_permute_and_resume_without_resetting_loading_or_cores() {
    let mut source = card([1.0, 2.0, 3000.0]);
    source.validity = source.validity.clone().with_quantity(
        LOADING_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        1.0,
        1.0,
    );
    let calcite = model(
        MineralReaction::CalciteBicarbonateNeutralization,
        "CaCO3(s)",
        source.clone(),
    );
    let bins = vec![
        population("a", calcite.clone(), 0.001, 0.25),
        population("b", calcite.clone(), 0.001, 0.75),
    ];
    let mut reversed = bins.clone();
    reversed.reverse();
    let mut split = batch(bins.clone(), 1.0);
    let mut permuted = batch(reversed, 1.0);
    let mut resumed = batch(bins, 1.0);
    let mut single = batch(vec![population("whole", calcite, 0.001, 1.0)], 1.0);
    // Exact finite-bath integral at y=.5: ln2 + 3 + 3(1-ln2).
    let t = 6.0 - 2.0 * core::f64::consts::LN_2;
    let result = advance(&mut split, t);
    let permutation = advance(&mut permuted, t);
    let unsplit = advance(&mut single, t);
    advance(&mut resumed, t * 0.5);
    let continuation = advance(&mut resumed, t * 0.5);
    assert_eq!(result.populations, permutation.populations);
    assert_eq!(split.inventory(), permuted.inventory());
    for p in &result.populations {
        close(p.conversion.value(), 0.875, 4e-8);
        close(p.core_radius.value(), 0.0005, 2e-11);
        assert!(p.remaining.value() > 0.0);
    }
    close(amount(split.inventory(), "H+(aq)"), 0.125, 4e-8);
    close(unsplit.populations[0].conversion.value(), 0.875, 4e-8);
    for (a, b) in result.populations.iter().zip(&continuation.populations) {
        close(a.conversion.value(), b.conversion.value(), 6e-8);
        close(a.remaining.value(), b.remaining.value(), 6e-8);
    }
    for ((p, initial), final_point) in result
        .models()
        .iter()
        .zip(&result.initial_points)
        .zip(&result.final_points)
    {
        assert_eq!(p.model().source(), &source);
        assert_eq!(initial.axes()[LOADING_AXIS], 1.0);
        assert_eq!(final_point.axes()[LOADING_AXIS], 1.0);
    }
}

#[test]
fn distinct_minerals_compete_with_their_own_stoichiometry_and_keep_finite_remainders() {
    let copper = model(
        MineralReaction::CupricOxideAcidLeach,
        "CuO(s)",
        card([1.5, 0.0, 0.0]),
    );
    let calcite = model(
        MineralReaction::CalciteBicarbonateNeutralization,
        "CaCO3(s)",
        card([3.0, 0.0, 0.0]),
    );
    let mut bath = batch(
        vec![
            population("copper", copper, 0.001, 1.0),
            population("gangue", calcite, 0.001, 1.0),
        ],
        3.0,
    );
    // Both intrinsic F=X; acid demand is 2X+X, so X=1-exp(-3t).
    let half_life = core::f64::consts::LN_2 / 3.0;
    let first = advance(&mut bath, half_life);
    assert_eq!(first.cumulative.reactions.len(), 2);
    for p in &first.populations {
        close(p.conversion.value(), 0.5, 5e-8);
    }
    close(amount(bath.inventory(), "H+(aq)"), 1.5, 8e-8);
    for product in ["Cu[2+](aq)", "Ca[2+](aq)", "HCO3[-](aq)", "H2O(l)"] {
        close(amount(bath.inventory(), product), 0.5, 5e-8);
    }
    let fourth = advance(&mut bath, 3.0 * half_life);
    for p in &fourth.populations {
        close(p.conversion.value(), 15.0 / 16.0, 6e-8);
        assert!(p.remaining.value() > 0.0);
    }
    close(amount(bath.inventory(), "H+(aq)"), 3.0 / 16.0, 1e-7);
    close(fourth.cumulative_exposure.value(), 15.0 / 16.0, 6e-8);
}

#[test]
fn pure_diffusion_and_near_stoichiometric_early_time_have_independent_answers() {
    let diffusion = model(
        MineralReaction::CalciteBicarbonateNeutralization,
        "CaCO3(s)",
        card([0.0, 0.0, 3000.0]),
    );
    let mut bath = batch(vec![population("diffusion", diffusion, 0.001, 1.0)], 1.0);
    let result = advance(&mut bath, 3.0 * (1.0 - core::f64::consts::LN_2));
    close(result.populations[0].conversion.value(), 0.875, 4e-8);
    close(amount(bath.inventory(), "H+(aq)"), 0.125, 4e-8);

    let surface = model(
        MineralReaction::CalciteBicarbonateNeutralization,
        "CaCO3(s)",
        card([0.0, 1.0, 0.0]),
    );
    let mut early = batch(vec![population("surface", surface, 0.001, 0.9999)], 1.0);
    let settings = PopulationNumerics {
        max_work: 200_000,
        ..PopulationNumerics::default()
    };
    let t = 0.01_f64;
    let out = early.advance(Time::new(t), settings, || false).unwrap();
    let stoichiometric = 1.0 - (1.0 + 2.0 * t).powf(-1.5);
    let maintained = 1.0 - (1.0 - t).powi(3);
    let x = out.populations[0].conversion.value();
    assert!(x >= stoichiometric - 4e-8 && x < maintained);
    assert!(out.populations[0].remaining.value() > 0.0);
    assert!(amount(early.inventory(), "H+(aq)") > 0.0);
    assert!(out.work <= settings.max_work);
}

#[test]
fn absolute_source_clock_and_depleted_concentration_refuse_atomically() {
    let mut clock = card([3.0, 0.0, 0.0]);
    clock.validity = clock.validity.clone().with_quantity(
        ELAPSED_TIME_AXIS,
        QuantitySpec::dimensional(Time::DIMS),
        0.0,
        0.3,
    );
    let calcite = model(
        MineralReaction::CalciteBicarbonateNeutralization,
        "CaCO3(s)",
        clock,
    );
    let mut bath = batch(vec![population("aged", calcite, 0.001, 1.0)], 2.0);
    let first = advance(&mut bath, 0.2);
    assert_eq!(first.initial_points[0].axes()[ELAPSED_TIME_AXIS], 0.0);
    assert_eq!(first.final_points[0].axes()[ELAPSED_TIME_AXIS], 0.2);
    let accepted = snapshot(&bath);
    assert!(matches!(
        bath.advance(Time::new(0.2), PopulationNumerics::default(), || false),
        Err(PopulationBatchError::Core(
            ShrinkingCoreError::OutsideValidity
        ))
    ));
    assert_eq!(snapshot(&bath), accepted);

    let mut restricted = card([1.5, 0.0, 0.0]);
    restricted.validity = restricted.validity.clone().with_quantity(
        CONCENTRATION_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        2.99,
        3.0,
    );
    let copper = model(MineralReaction::CupricOxideAcidLeach, "CuO(s)", restricted);
    let mut domain = batch(vec![population("restricted", copper, 0.001, 1.0)], 3.0);
    let before = snapshot(&domain);
    assert!(matches!(
        domain.advance(Time::new(0.1), PopulationNumerics::default(), || false),
        Err(PopulationBatchError::Core(
            ShrinkingCoreError::OutsideValidity
        ))
    ));
    assert_eq!(snapshot(&domain), before);
}

// A conserved synthetic chain solely for the forbidden-feedback constructor
// test. These assigned solid roles/coefficients are not a chemical mechanism.
fn synthetic_water_step(
    name: &str,
    solid: &str,
    product: &str,
    solid_oxygen: i128,
) -> ShrinkingCoreReaction {
    let species = vec![id(solid), id("O2(aq)"), id(product)];
    let elements = ElementalMatrix::new(
        vec![ElementId::new("H").unwrap(), ElementId::new("O").unwrap()],
        species.clone(),
        vec![vec![2, 0, 2], vec![solid_oxygen, 2, solid_oxygen + 1]],
    )
    .unwrap();
    let reaction = ReactionId::new(name).unwrap();
    let stoichiometry = StoichiometricMatrix::new(
        species.clone(),
        vec![reaction.clone()],
        vec![vec![-2], vec![-1], vec![2]],
    )
    .unwrap();
    let charges = ChargeVector::new(species, vec![0, 0, 0]).unwrap();
    ShrinkingCoreReaction::from_card(
        ConservedReaction::new(&elements, &stoichiometry, &charges, reaction).unwrap(),
        id(solid),
        id("O2(aq)"),
        card([1.0, 0.0, 0.0]),
    )
    .unwrap()
}

#[test]
fn population_assignment_identity_and_reaction_feedback_are_admitted_explicitly() {
    let bins = copper_bins();
    let original = inventory(&bins, 3.0);
    let mut wrong_stock = original.clone();
    wrong_stock
        .iter_mut()
        .find(|(id, _)| id.as_str() == "CuO(s)")
        .unwrap()
        .1 = Amount::new(0.9);
    let construct = |populations, inventory: &[(SpeciesId, Amount)]| {
        ShrinkingCoreBatch::new(
            populations,
            id("H+(aq)"),
            inventory,
            Volume::new(1.0),
            Temperature::new(298.15),
        )
    };
    assert!(matches!(
        construct(bins.clone(), &wrong_stock),
        Err(PopulationBatchError::Invalid(_))
    ));
    assert!(matches!(
        construct(vec![bins[0].clone(), bins[0].clone()], &original),
        Err(PopulationBatchError::Invalid(_))
    ));
    let mut duplicated = original.clone();
    duplicated.push(original[0].clone());
    assert!(matches!(
        construct(bins, &duplicated),
        Err(PopulationBatchError::Network(_))
    ));

    let hydration = ShrinkingCoreReaction::from_card(
        MineralReaction::LimeHydration.build().unwrap(),
        id("CaO(s)"),
        id("H2O(l)"),
        card([1.0, 0.0, 0.0]),
    )
    .unwrap();
    assert!(matches!(
        construct(
            vec![population("wrong-reagent", hydration, 0.001, 1.0)],
            &original
        ),
        Err(PopulationBatchError::Invalid(
            "every population must consume the same reagent"
        ))
    ));
    let feedback = vec![
        population(
            "first",
            synthetic_water_step("first", "H2(s)", "H2O(s)", 0),
            0.001,
            1.0,
        ),
        population(
            "second",
            synthetic_water_step("second", "H2O(s)", "H2O2(s)", 1),
            0.001,
            1.0,
        ),
    ];
    assert!(matches!(
        ShrinkingCoreBatch::new(
            feedback,
            id("O2(aq)"),
            &[],
            Volume::new(1.0),
            Temperature::new(298.15)
        ),
        Err(PopulationBatchError::Invalid(
            "reaction products cannot feed modeled reactants"
        ))
    ));
}

#[test]
fn work_cancellation_zero_duration_and_zero_reagent_preserve_accepted_state() {
    let mut bath = batch(copper_bins(), 3.0);
    advance(&mut bath, 0.1);
    let accepted = snapshot(&bath);
    let limited = PopulationNumerics {
        max_work: 1,
        ..PopulationNumerics::default()
    };
    assert!(matches!(
        bath.advance(Time::new(0.2), limited, || false),
        Err(PopulationBatchError::Core(ShrinkingCoreError::Budget))
    ));
    assert_eq!(snapshot(&bath), accepted);
    let mut polls = 0;
    assert!(matches!(
        bath.advance(Time::new(0.2), PopulationNumerics::default(), || {
            polls += 1;
            polls > 7
        }),
        Err(PopulationBatchError::Core(ShrinkingCoreError::Cancelled))
    ));
    assert!(polls > 7);
    assert_eq!(snapshot(&bath), accepted);
    let identity = advance(&mut bath, 0.0);
    assert_eq!(identity.exposure_increment.value(), 0.0);
    assert_eq!(snapshot(&bath), accepted);
    let mut exhausted = batch(copper_bins(), 0.0);
    let stocks = exhausted.inventory().to_vec();
    let idle = advance(&mut exhausted, 1.0);
    assert_eq!(idle.exposure_increment.value(), 0.0);
    assert_eq!(exhausted.inventory(), stocks);
    assert!(idle.populations.iter().all(|p| p.conversion.value() == 0.0));
}
