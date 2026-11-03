//! G0/G1/G3/G4 checks of simultaneous, explicitly charge-balanced CIL.
//! Chemical fixtures are shared with the original cyanidation owner; sorbent
//! cards below are authored arithmetic fixtures, not measured carbon or ore.

use super::*;
use crate::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, AdsorbentLoading, AqueousLangmuir, CAPACITY_PARAMETER,
    CONCENTRATION_AXIS, HALF_SATURATION_PARAMETER,
    flow::{LOCAL_LDF_LAW, LocalLangmuirTransfer, RetainedSorbent},
};
use crate::cyanide_leach::tests::{
    card as leach_card, charge,
    competition::{competing_batch, copper_charge},
    options, system,
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, Frequency, Mass, QuantitySpec, Time, Volume,
    chemistry::ElementId as AnalyticalElementId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

pub(crate) fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}",
    );
}

#[test]
fn actual_dissolved_oxygen_source_exit_is_atomic_after_a_valid_interval() {
    let template = cil(charge(1.0), 1.0);
    let local = template.state().populations[0].transfer.model().clone();
    let mut equilibrium_source = local.equilibrium().source().clone();
    equilibrium_source.validity = equilibrium_source.validity.with_quantity(
        OXYGEN_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        0.24,
        0.25,
    );
    let equilibrium =
        AqueousLangmuir::from_card(AnalyticalElementId::new("Au").unwrap(), equilibrium_source)
            .unwrap();
    let transfer = LocalLangmuirTransfer::from_card(equilibrium, local.source().clone()).unwrap();
    let sorbent = RetainedSorbent::new(
        "oxygen-qualified-carbon",
        transfer,
        Mass::new(1.0),
        Amount::new(0.0),
        QueryPoint::new(),
    )
    .unwrap();
    let mut bath = CarbonInLeachBatch::from_batch(
        original(charge(1.0)),
        IonPairTransfer::sodium().unwrap(),
        vec![CarbonPopulation::new(CarbonRole::ProductCarbon, sorbent).unwrap()],
    )
    .unwrap();
    let first = bath
        .advance(Time::new(0.02), numerics(1.0), || false)
        .unwrap();
    let transfer = &first.state.populations[0].transfer;
    for point in [&transfer.conditions, &transfer.equilibrium.conditions] {
        close(
            point.axes()[ELAPSED_TIME_AXIS],
            first.state.elapsed.value(),
            0.0,
        );
        close(
            point.axes()[OXYGEN_AXIS],
            first.state.solution.dissolved_oxygen.value(),
            0.0,
        );
    }
    assert!(first.state.solution.dissolved_oxygen.value() >= 0.24);
    let before = format!("{:?}", bath.state());
    // The leach card permits this clock and DO; the original carbon isotherm
    // does not permit the lower oxygen reached along the requested trajectory.
    assert!(
        bath.advance(Time::new(0.08), numerics(1.0), || false)
            .is_err()
    );
    assert_eq!(format!("{:?}", bath.state()), before);
}

fn local_source(
    name: &str,
    law: &str,
    parameters: &[(&str, f64, Dims)],
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    let artifact = hash_domain("simultaneous-cil-arithmetic-test-v1", name.as_bytes());
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: parameters
            .iter()
            .map(|&(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: format!("Authored {name} arithmetic; no ore/carbon calibration"),
            license: "project test fixture".into(),
            artifact: Some(artifact),
        },
    }
}

pub(crate) fn population(
    name: &str,
    role: CarbonRole,
    mass: f64,
    preloaded: f64,
    capacity: f64,
    half: f64,
    rate: f64,
    concentration_ceiling: f64,
) -> CarbonPopulation {
    let domain = ValidityDomain::unconstrained()
        .with_quantity(
            CONCENTRATION_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            0.0,
            concentration_ceiling,
        )
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
            298.15,
        )
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), -1.0, 15.0)
        .with_quantity(
            FREE_CYANIDE_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            0.0,
            10.0,
        )
        .with_quantity(
            COUNTERION_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            0.0,
            10.0,
        )
        .with_quantity(
            OXYGEN_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            0.0,
            10.0,
        );
    let equilibrium = AqueousLangmuir::from_card(
        AnalyticalElementId::new("Au").unwrap(),
        local_source(
            name,
            AQUEOUS_LANGMUIR_LAW,
            &[
                (CAPACITY_PARAMETER, capacity, AdsorbentLoading::DIMS),
                (HALF_SATURATION_PARAMETER, half, AmountConcentration::DIMS),
            ],
            domain.clone(),
        ),
    )
    .unwrap();
    let rate_source = local_source(
        name,
        LOCAL_LDF_LAW,
        &[(
            crate::aqueous_adsorption::flow::RATE_PARAMETER,
            rate,
            Frequency::DIMS,
        )],
        domain
            .with_quantity(
                crate::aqueous_adsorption::flow::LOADING_AXIS,
                QuantitySpec::dimensional(AdsorbentLoading::DIMS),
                0.0,
                capacity,
            )
            .with_quantity(
                crate::aqueous_adsorption::flow::VOLUME_AXIS,
                QuantitySpec::dimensional(Volume::DIMS),
                0.001,
                10.0,
            ),
    );
    CarbonPopulation::new(
        role,
        RetainedSorbent::new(
            name,
            LocalLangmuirTransfer::from_card(equilibrium, rate_source).unwrap(),
            Mass::new(mass),
            Amount::new(preloaded),
            QueryPoint::new(),
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn dynamic_chemistry_reaches_both_source_points_and_owned_extras_are_refused() {
    let initial = charge(1.0);
    let mut bath = cil(initial.clone(), 1.0);
    let initial_free_cn = bath.state().solution.free_cyanide.value();
    let report = bath
        .advance(Time::new(0.3), numerics(1.0), || false)
        .unwrap();
    let state = &report.state;
    assert!(state.solution.free_cyanide.value() < initial_free_cn);
    assert!(state.dissolved_counterion.value() < 1.0);
    for population in &state.populations {
        for point in [
            &population.transfer.conditions,
            &population.transfer.equilibrium.conditions,
        ] {
            for (axis, expected, dims) in [
                (PH_AXIS, state.solution.speciation.ph.value(), Dims::NONE),
                (ELAPSED_TIME_AXIS, state.elapsed.value(), Time::DIMS),
                (
                    OXYGEN_AXIS,
                    state.solution.dissolved_oxygen.value() / initial.solution_volume.value(),
                    AmountConcentration::DIMS,
                ),
                (
                    FREE_CYANIDE_AXIS,
                    state.solution.free_cyanide.value() / initial.solution_volume.value(),
                    AmountConcentration::DIMS,
                ),
                (
                    COUNTERION_AXIS,
                    state.dissolved_counterion.value() / initial.solution_volume.value(),
                    AmountConcentration::DIMS,
                ),
                (
                    CONCENTRATION_AXIS,
                    state.solution.gold_complex.value() / initial.solution_volume.value(),
                    AmountConcentration::DIMS,
                ),
            ] {
                close(
                    *point.axes().get(axis).expect("actual source coordinate"),
                    expected,
                    0.0,
                );
                assert_eq!(
                    point.axis_quantities().get(axis),
                    Some(&QuantitySpec::dimensional(dims))
                );
            }
        }
    }
    let model = state.populations[0].transfer.model().clone();
    // Neither original source law admits a batch-duration validity axis, so
    // the retained-source owner rejects this override before CIL attachment.
    assert!(
        RetainedSorbent::new(
            "source-clock-override",
            model.clone(),
            Mass::new(1.0),
            Amount::new(0.0),
            QueryPoint::new()
                .with_quantity(
                    ELAPSED_TIME_AXIS,
                    QuantitySpec::dimensional(Time::DIMS),
                    0.1,
                )
                .unwrap(),
        )
        .is_err()
    );
    for (axis, dims) in [
        (PH_AXIS, Dims::NONE),
        (FREE_CYANIDE_AXIS, AmountConcentration::DIMS),
        (COUNTERION_AXIS, AmountConcentration::DIMS),
        (OXYGEN_AXIS, AmountConcentration::DIMS),
    ] {
        let caller_override = QueryPoint::new()
            .with_quantity(axis, QuantitySpec::dimensional(dims), 0.1)
            .unwrap();
        let sorbent = RetainedSorbent::new(
            "owned-axis-override",
            model.clone(),
            Mass::new(1.0),
            Amount::new(0.0),
            caller_override,
        )
        .unwrap();
        assert!(
            CarbonPopulation::new(CarbonRole::ProductCarbon, sorbent).is_err(),
            "caller must not override {axis}"
        );
    }
}

#[test]
fn fresh_cyanide_and_caustic_doses_survive_attachment_and_enter_all_accounts() {
    use crate::cyanide_leach::dosing::DissolvedSodiumReagent;
    let initial = copper_charge(0.2);
    let mut fresh = competing_batch(initial.clone(), 0.2);
    fresh
        .dose_reagent(
            DissolvedSodiumReagent::SodiumCyanide,
            Amount::new(0.4),
            10_000_000,
            || false,
        )
        .unwrap();
    fresh
        .dose_reagent(
            DissolvedSodiumReagent::SodiumHydroxide,
            Amount::new(0.2),
            10_000_000,
            || false,
        )
        .unwrap();
    let inputs = fresh.reagent_inputs();
    let prepared_solution = fresh.state().clone();
    assert_eq!(
        fresh.initial_charge().uncomplexed_cyanide,
        initial.uncomplexed_cyanide
    );
    let mut bath =
        CarbonInLeachBatch::from_batch(fresh, IonPairTransfer::sodium().unwrap(), populations(1.0))
            .unwrap();
    assert_eq!(bath.reagent_inputs(), inputs);
    close(
        bath.state().solution.free_cyanide.value(),
        prepared_solution.free_cyanide.value(),
        2e-12,
    );
    close(
        bath.state().solution.speciation.ph.value(),
        prepared_solution.speciation.ph.value(),
        2e-11,
    );
    close(bath.state().dissolved_counterion.value(), 1.6, 4e-16);
    let report = bath
        .advance(Time::new(0.3), numerics(1.0), || false)
        .unwrap();
    assert_eq!(bath.reagent_inputs(), inputs);
    close(inputs.sodium_cyanide.value(), 0.4, 0.0);
    close(inputs.sodium_hydroxide.value(), 0.2, 0.0);
    close(inputs.sodium.value(), 0.6, 2e-16);
    // Reconstruct original plus delivered totals independently for the account;
    // neither the runtime's original charge nor its reagent history is rebased.
    let mut account_target = initial.clone();
    account_target.uncomplexed_cyanide =
        Amount::new(initial.uncomplexed_cyanide.value() + inputs.sodium_cyanide.value());
    let sodium = account_target
        .spectators
        .iter_mut()
        .find(|ion| ion.species.as_str() == "Na[+](aq)")
        .unwrap();
    sodium.amount = Amount::new(sodium.amount.value() + inputs.sodium.value());
    accounts(
        &report.state,
        &account_target,
        0.0,
        &[("activated-carbon", 0.6), ("natural-preg-robber", 0.3)],
    );
    report_accounts(&report);
}

pub(crate) fn populations(scale: f64) -> Vec<CarbonPopulation> {
    vec![
        population(
            "activated-carbon",
            CarbonRole::ProductCarbon,
            scale,
            0.0,
            0.6,
            0.1,
            1.4,
            10.0,
        ),
        population(
            "natural-preg-robber",
            CarbonRole::PregRobbing,
            scale,
            0.0,
            0.3,
            0.2,
            0.35,
            10.0,
        ),
    ]
}

pub(crate) fn original(c: CyanidationCharge) -> GoldCyanidationBatch {
    GoldCyanidationBatch::new(
        GoldCyanidation::from_card(system(false), "uncomplexed-cyanide", leach_card(4.0)).unwrap(),
        c,
        QueryPoint::new(),
        options(),
    )
    .unwrap()
}

pub(crate) fn cil(c: CyanidationCharge, scale: f64) -> CarbonInLeachBatch {
    CarbonInLeachBatch::from_batch(
        original(c),
        IonPairTransfer::sodium().unwrap(),
        populations(scale),
    )
    .unwrap()
}

pub(crate) fn numerics(scale: f64) -> CilNumerics {
    CilNumerics {
        maximum_step: Time::new(0.1),
        absolute_amount_error: Amount::new(2e-9 * scale),
        relative_amount_error: Dimensionless::new(1e-8),
        maximum_attempts: 8192,
        maximum_work: 10_000_000,
    }
}

fn sum_adsorbed(state: &CilState) -> f64 {
    state.populations.iter().map(|p| p.adsorbed.value()).sum()
}

fn member(state: &CilState, name: &str) -> f64 {
    state
        .populations
        .iter()
        .find(|p| p.id == name)
        .unwrap()
        .adsorbed
        .value()
}

fn accounts(
    state: &CilState,
    initial: &CyanidationCharge,
    initial_adsorbed: f64,
    capacities: &[(&str, f64)],
) {
    let solution = &state.solution;
    let bound = solution
        .competition
        .as_ref()
        .map_or(0.0, |c| c.bound_ligand.value());
    let adsorbed = sum_adsorbed(state);
    let sodium = initial
        .spectators
        .iter()
        .find(|s| s.species.as_str() == "Na[+](aq)")
        .unwrap();
    let scale =
        initial.accessible_gold.value() + initial.gold_complex.value() + initial_adsorbed + 1.0;
    close(
        solution.accessible_gold.value() + solution.gold_complex.value() + adsorbed,
        initial.accessible_gold.value() + initial.gold_complex.value() + initial_adsorbed,
        4e-13 * scale,
    );
    close(
        solution.uncomplexed_cyanide.value()
            + bound
            + 2.0 * (solution.gold_complex.value() + adsorbed),
        initial.uncomplexed_cyanide.value()
            + 2.0 * (initial.gold_complex.value() + initial_adsorbed),
        4e-9 * scale,
    );
    close(
        state.dissolved_counterion.value() + adsorbed,
        sodium.amount.value() + initial_adsorbed,
        4e-13 * scale,
    );
    close(
        solution.dissolved_oxygen.value() + 0.25 * solution.gold_dissolved.value(),
        initial.dissolved_oxygen.value(),
        4e-13 * scale,
    );
    close(
        solution
            .speciation
            .spectators
            .iter()
            .find(|s| s.species.as_str() == "Na[+](aq)")
            .unwrap()
            .amount
            .value(),
        state.dissolved_counterion.value(),
        4e-13 * scale,
    );
    close(
        solution.speciation.charge_balance_residual.value(),
        0.0,
        4e-10 * scale,
    );
    for population in &state.populations {
        let capacity = capacities
            .iter()
            .find(|(id, _)| *id == population.id)
            .unwrap()
            .1;
        assert!(population.adsorbed.value() >= 0.0);
        assert!(population.vacant_capacity.value() >= 0.0);
        close(
            population.adsorbed.value() + population.vacant_capacity.value(),
            capacity,
            4e-13 * scale,
        );
    }
    assert!(solution.gold_complex.value() >= 0.0);
    assert!(state.dissolved_counterion.value() >= 0.0);
}

fn report_accounts(report: &CilAdvance) {
    for residual in [
        report.conservation.gold,
        report.conservation.counterion,
        report.conservation.dissolved_oxygen,
        report.conservation.maximum_site_residual,
    ] {
        close(residual.value(), 0.0, 2e-12);
    }
    close(report.conservation.cyanide.value(), 0.0, 1e-8);
    assert!(report.accepted_steps > 0);
    assert!(report.work > 0);
}

#[test]
fn simultaneous_dissolution_and_independently_declared_carbon_populations() {
    let initial = charge(1.0);
    let mut bath = cil(initial.clone(), 1.0);
    let mut leach_only = original(initial.clone());
    let report = bath
        .advance(Time::new(0.6), numerics(1.0), || false)
        .unwrap();
    let reference = leach_only
        .advance(
            Time::new(0.6),
            CyanidationNumerics {
                maximum_step: Time::new(0.1),
                absolute_gold_error: Amount::new(2e-10),
                relative_gold_error: Dimensionless::new(1e-9),
                ..CyanidationNumerics::default()
            },
            || false,
        )
        .unwrap();
    assert!(report.state.solution.gold_dissolved.value() > 0.1);
    assert!(report.state.solution.accessible_gold.value() > 0.0);
    assert!(report.state.solution.gold_complex.value() > 0.0);
    let product = member(&report.state, "activated-carbon");
    let preg = member(&report.state, "natural-preg-robber");
    assert!(product > preg && preg > 0.0);
    close(report.interval_product_capture.value(), product, 2e-13);
    close(report.interval_preg_robbing_capture.value(), preg, 2e-13);
    close(
        report.interval_gold_dissolved.value(),
        report.state.solution.gold_dissolved.value(),
        2e-13,
    );
    assert!(report.state.solution.gold_complex.value() < reference.state.gold_complex.value());
    // The declared rate has no product inhibition: neutral pair transfer must
    // not create free cyanide, change pH, or spuriously accelerate dissolution.
    close(
        report.state.solution.gold_dissolved.value(),
        reference.state.gold_dissolved.value(),
        8e-8,
    );
    close(
        report.state.solution.free_cyanide.value(),
        reference.state.free_cyanide.value(),
        2e-7,
    );
    close(
        report.state.solution.speciation.ph.value(),
        reference.state.speciation.ph.value(),
        2e-7,
    );
    accounts(
        &report.state,
        &initial,
        0.0,
        &[("activated-carbon", 0.6), ("natural-preg-robber", 0.3)],
    );
    report_accounts(&report);
}

#[test]
fn preloaded_desorption_has_an_independent_exponential_limit_and_pair_account() {
    let mut initial = charge(1.0);
    initial.accessible_gold = Amount::new(0.0);
    let reference = original(initial.clone()).state().clone();
    let mut bath = CarbonInLeachBatch::from_batch(
        original(initial.clone()),
        IonPairTransfer::sodium().unwrap(),
        vec![population(
            "preloaded-product",
            CarbonRole::ProductCarbon,
            1.0,
            0.2,
            2.0,
            1e12,
            0.7,
            10.0,
        )],
    )
    .unwrap();
    let report = bath
        .advance(Time::new(1.0), numerics(1.0), || false)
        .unwrap();
    let adsorbed = member(&report.state, "preloaded-product");
    // With half concentration 1e12, the total reload correction is <4e-13 mol.
    close(adsorbed, 0.2 * (-0.7_f64).exp(), 8e-9);
    assert!(report.interval_product_capture.value() < 0.0);
    close(
        report.interval_product_capture.value(),
        adsorbed - 0.2,
        3e-14,
    );
    close(
        report.state.solution.gold_complex.value(),
        0.2 - adsorbed,
        3e-14,
    );
    close(
        report.state.dissolved_counterion.value(),
        1.2 - adsorbed,
        3e-14,
    );
    close(report.state.solution.gold_dissolved.value(), 0.0, 0.0);
    close(
        report.state.solution.free_cyanide.value(),
        reference.free_cyanide.value(),
        2e-12,
    );
    close(
        report.state.solution.speciation.ph.value(),
        reference.speciation.ph.value(),
        2e-11,
    );
    accounts(&report.state, &initial, 0.2, &[("preloaded-product", 2.0)]);
    report_accounts(&report);
}

#[test]
fn active_copper_preserves_metal_mass_action_and_all_cil_pools() {
    let initial = copper_charge(0.2);
    let mut control = competing_batch(initial.clone(), 0.2);
    let mut bath = CarbonInLeachBatch::from_batch(
        competing_batch(initial.clone(), 0.2),
        IonPairTransfer::sodium().unwrap(),
        populations(1.0),
    )
    .unwrap();
    let report = bath
        .advance(Time::new(0.4), numerics(1.0), || false)
        .unwrap();
    let reference = control
        .advance(Time::new(0.4), CyanidationNumerics::default(), || false)
        .unwrap();
    let metal = &report.state.solution.competition.as_ref().unwrap().metals[0];
    close(
        metal.members.iter().map(|m| m.amount.value()).sum(),
        0.2,
        2e-12,
    );
    let activity = |name: &str| {
        metal
            .members
            .iter()
            .find(|m| m.species.as_str() == name)
            .unwrap()
            .activity
            .value()
    };
    let free_cn =
        report.state.solution.free_cyanide.value() / initial.amount_at_unit_activity.value();
    close(
        activity("CuCN(aq)") / (activity("Cu[+](aq)") * free_cn),
        10.0,
        3e-8,
    );
    close(
        activity("Cu(CN)2[-](aq)") / (activity("CuCN(aq)") * free_cn),
        5.0,
        3e-8,
    );
    close(
        report.state.solution.gold_dissolved.value(),
        reference.state.gold_dissolved.value(),
        2e-7,
    );
    close(
        report.state.solution.free_cyanide.value(),
        reference.state.free_cyanide.value(),
        4e-7,
    );
    assert!(sum_adsorbed(&report.state) > 0.0);
    accounts(
        &report.state,
        &initial,
        0.0,
        &[("activated-carbon", 0.6), ("natural-preg-robber", 0.3)],
    );
    report_accounts(&report);
}

#[test]
fn continuation_and_extensive_scaling_preserve_the_coupled_trajectory() {
    let initial = charge(1.0);
    let mut whole = cil(initial.clone(), 1.0);
    let mut resumed = cil(initial.clone(), 1.0);
    let whole_report = whole
        .advance(Time::new(0.6), numerics(1.0), || false)
        .unwrap();
    let first = resumed
        .advance(Time::new(0.2), numerics(1.0), || false)
        .unwrap();
    let second = resumed
        .advance(Time::new(0.4), numerics(1.0), || false)
        .unwrap();
    close(second.state.elapsed.value(), 0.6, 2e-16);
    close(
        second.state.solution.gold_dissolved.value(),
        whole_report.state.solution.gold_dissolved.value(),
        2e-7,
    );
    close(
        first.interval_product_capture.value() + second.interval_product_capture.value(),
        member(&second.state, "activated-carbon"),
        2e-13,
    );
    close(
        first.interval_preg_robbing_capture.value() + second.interval_preg_robbing_capture.value(),
        member(&second.state, "natural-preg-robber"),
        2e-13,
    );
    for p in &whole_report.state.populations {
        close(member(&second.state, &p.id), p.adsorbed.value(), 2e-7);
    }
    let factor = 1.5;
    let mut scaled = initial.clone();
    scaled.accessible_gold = Amount::new(initial.accessible_gold.value() * factor);
    scaled.uncomplexed_cyanide = Amount::new(initial.uncomplexed_cyanide.value() * factor);
    scaled.dissolved_oxygen = Amount::new(initial.dissolved_oxygen.value() * factor);
    scaled.gold_complex = Amount::new(initial.gold_complex.value() * factor);
    scaled.amount_at_unit_activity = Amount::new(initial.amount_at_unit_activity.value() * factor);
    scaled.solution_volume = Volume::new(initial.solution_volume.value() * factor);
    for ion in &mut scaled.spectators {
        ion.amount = Amount::new(ion.amount.value() * factor);
    }
    let scaled_report = cil(scaled.clone(), factor)
        .advance(Time::new(0.6), numerics(factor), || false)
        .unwrap();
    close(
        scaled_report.state.solution.gold_dissolved.value() / factor,
        whole_report.state.solution.gold_dissolved.value(),
        2e-7,
    );
    close(
        scaled_report.state.solution.speciation.ph.value(),
        whole_report.state.solution.speciation.ph.value(),
        2e-7,
    );
    for p in &whole_report.state.populations {
        close(
            member(&scaled_report.state, &p.id) / factor,
            p.adsorbed.value(),
            2e-7,
        );
    }
    accounts(
        &scaled_report.state,
        &scaled,
        0.0,
        &[
            ("activated-carbon", 0.6 * factor),
            ("natural-preg-robber", 0.3 * factor),
        ],
    );
}

#[test]
fn source_exit_cancellation_and_work_exhaustion_publish_no_partial_state() {
    let mut source_limited = CarbonInLeachBatch::from_batch(
        original(charge(1.0)),
        IonPairTransfer::sodium().unwrap(),
        vec![population(
            "narrow-source",
            CarbonRole::ProductCarbon,
            1.0,
            0.0,
            0.6,
            0.1,
            1.4,
            1e-8,
        )],
    )
    .unwrap();
    let before = format!("{:?}", source_limited.state());
    assert!(
        source_limited
            .advance(Time::new(0.1), numerics(1.0), || false)
            .is_err()
    );
    assert_eq!(format!("{:?}", source_limited.state()), before);

    let mut cancelled = cil(charge(1.0), 1.0);
    let before = format!("{:?}", cancelled.state());
    let mut polls = 0;
    assert!(
        cancelled
            .advance(Time::new(0.3), numerics(1.0), || {
                polls += 1;
                polls >= 8
            })
            .is_err()
    );
    assert!(polls >= 8);
    assert_eq!(format!("{:?}", cancelled.state()), before);

    let mut budgeted = cil(charge(1.0), 1.0);
    let before = format!("{:?}", budgeted.state());
    let budget = CilNumerics {
        maximum_work: 1,
        ..numerics(1.0)
    };
    assert!(budgeted.advance(Time::new(0.3), budget, || false).is_err());
    assert_eq!(format!("{:?}", budgeted.state()), before);
    // Failure leaves the continuation/error history usable, too.
    let after = budgeted
        .advance(Time::new(0.1), numerics(1.0), || false)
        .unwrap();
    let clean = cil(charge(1.0), 1.0)
        .advance(Time::new(0.1), numerics(1.0), || false)
        .unwrap();
    assert_eq!(format!("{:?}", after.state), format!("{:?}", clean.state));
    assert_eq!(
        after.accumulated_error_estimate,
        clean.accumulated_error_estimate
    );
}

#[test]
fn attachment_rejects_insufficient_counterion_duplicate_ids_and_used_baths() {
    assert!(
        CarbonInLeachBatch::from_batch(
            original(charge(0.5)),
            IonPairTransfer::sodium().unwrap(),
            populations(1.0),
        )
        .is_err()
    );
    assert!(
        CarbonInLeachBatch::from_batch(
            original(charge(1.0)),
            IonPairTransfer::sodium().unwrap(),
            vec![
                population(
                    "duplicate",
                    CarbonRole::ProductCarbon,
                    1.0,
                    0.0,
                    0.6,
                    0.1,
                    1.4,
                    10.0
                ),
                population(
                    "duplicate",
                    CarbonRole::PregRobbing,
                    1.0,
                    0.0,
                    0.3,
                    0.2,
                    0.35,
                    10.0
                ),
            ],
        )
        .is_err()
    );
    let mut used = original(charge(1.0));
    used.advance(Time::new(0.01), CyanidationNumerics::default(), || false)
        .unwrap();
    assert!(
        CarbonInLeachBatch::from_batch(used, IonPairTransfer::sodium().unwrap(), populations(1.0),)
            .is_err()
    );
}


#[test]
fn live_cil_reagent_doses_respeciate_without_moving_carbon_or_progress() {
    use crate::cyanide_leach::dosing::DissolvedSodiumReagent::{SodiumCyanide, SodiumHydroxide};

    let initial = charge(1.0);
    let mut bath = cil(initial, 1.0);
    bath.advance(Time::new(0.2), numerics(1.0), || false)
        .unwrap();
    let before = bath.state().clone();
    let before_inputs = bath.reagent_inputs();

    let cyanide = bath
        .dose_reagent(SodiumCyanide, Amount::new(0.2), 10_000_000, || false)
        .unwrap();
    assert_eq!(cyanide.state.elapsed, before.elapsed);
    assert_eq!(
        cyanide.state.solution.gold_dissolved,
        before.solution.gold_dissolved
    );
    for (after, old) in cyanide.state.populations.iter().zip(&before.populations) {
        assert_eq!(after.adsorbed, old.adsorbed);
        assert_eq!(after.vacant_capacity, old.vacant_capacity);
    }
    assert!(
        cyanide.state.solution.free_cyanide.value()
            > before.solution.free_cyanide.value()
    );
    close(
        cyanide.inputs.sodium_cyanide.value(),
        before_inputs.sodium_cyanide.value() + 0.2,
        2e-16,
    );
    close(
        cyanide.inputs.sodium.value(),
        before_inputs.sodium.value() + 0.2,
        2e-16,
    );
    report_accounts(&CilAdvance {
        layout: cyanide.layout.clone(),
        state: cyanide.state.clone(),
        interval_gold_dissolved: Amount::new(0.0),
        interval_product_capture: Amount::new(0.0),
        interval_preg_robbing_capture: Amount::new(0.0),
        conservation: cyanide.conservation,
        accepted_steps: 0,
        rejected_steps: 0,
        work: cyanide.work,
        maximum_local_error: Amount::new(0.0),
        accumulated_error_estimate: Amount::new(0.0),
    });

    let ph_before = bath.state().solution.speciation.ph.value();
    let hydroxide = bath
        .dose_reagent(SodiumHydroxide, Amount::new(0.1), 10_000_000, || false)
        .unwrap();
    assert!(hydroxide.state.solution.speciation.ph.value() > ph_before);
    close(
        hydroxide.inputs.sodium_hydroxide.value(),
        before_inputs.sodium_hydroxide.value() + 0.1,
        2e-16,
    );
    for residual in [
        hydroxide.conservation.gold,
        hydroxide.conservation.counterion,
        hydroxide.conservation.dissolved_oxygen,
        hydroxide.conservation.maximum_site_residual,
    ] {
        close(residual.value(), 0.0, 2e-12);
    }
    close(hydroxide.conservation.cyanide.value(), 0.0, 1e-8);
}
