//! G0/G1/G3 analytical flow/adsorption limits and complete metal accounts.
use super::*;

#[path = "tests/instantaneous.rs"]
mod instantaneous;
use crate::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, CAPACITY_PARAMETER, HALF_SATURATION_PARAMETER,
    kinetics::{
        AdsorptionBatch, AdsorptionNumerics, ELAPSED_TIME_AXIS, LDF_LAW, LangmuirKinetics,
        SorbentPopulation,
    },
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::{Dims, chemistry::ElementId};

fn source(
    law: &str,
    parameters: &[(&str, f64, Dims)],
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    let artifact = hash_domain(
        "continuous-carbon-test-v1",
        b"synthetic arithmetic, not measured carbon behavior",
    );
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
            source: "Authored ideal-contact arithmetic fixture, not ore calibration".into(),
            license: "project fixture".into(),
            artifact: Some(artifact),
        },
    }
}
fn equilibrium_source(half: f64) -> ConstitutiveModelCard {
    source(
        AQUEOUS_LANGMUIR_LAW,
        &[
            (CAPACITY_PARAMETER, 2.0, AdsorbentLoading::DIMS),
            (HALF_SATURATION_PARAMETER, half, AmountConcentration::DIMS),
        ],
        ValidityDomain::unconstrained()
            .with_quantity(CONCENTRATION_AXIS, concentration_spec(), 0.0, 100.0)
            .with_quantity(TEMPERATURE_AXIS, temperature_spec(), 290.0, 310.0)
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0, 12.0),
    )
}
fn equilibrium(half: f64) -> AqueousLangmuir {
    AqueousLangmuir::from_card(ElementId::new("Au").unwrap(), equilibrium_source(half)).unwrap()
}
fn rate_source(rate: f64) -> ConstitutiveModelCard {
    let validity = equilibrium_source(1.0)
        .validity
        .with_quantity(
            LOADING_AXIS,
            QuantitySpec::dimensional(AdsorbentLoading::DIMS),
            0.0,
            2.0,
        )
        .with_quantity(
            VOLUME_AXIS,
            QuantitySpec::dimensional(Volume::DIMS),
            0.001,
            1000.0,
        );
    source(
        LOCAL_LDF_LAW,
        &[(RATE_PARAMETER, rate, Frequency::DIMS)],
        validity,
    )
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 11.0)
        .unwrap()
}
fn population(id: &str, mass: f64, loaded: f64, rate: f64, half: f64) -> RetainedSorbent {
    RetainedSorbent::new(
        id,
        LocalLangmuirTransfer::from_card(equilibrium(half), rate_source(rate)).unwrap(),
        Mass::new(mass),
        Amount::new(loaded),
        conditions(),
    )
    .unwrap()
}
fn tank(populations: Vec<RetainedSorbent>, dissolved: f64) -> AdsorptionCstr {
    AdsorptionCstr::new(
        populations,
        Amount::new(dissolved),
        Volume::new(2.0),
        Temperature::new(298.15),
    )
    .unwrap()
}
fn advance(
    tank: &mut AdsorptionCstr,
    q: f64,
    inlet: f64,
    t: f64,
    budget: AdsorptionFlowBudget,
) -> AdsorptionFlowAdvance {
    tank.advance(
        VolumetricFlowRate::new(q),
        AmountConcentration::new(inlet),
        Time::new(t),
        budget,
        || false,
    )
    .unwrap()
}
fn fixed_step(step: f64) -> AdsorptionFlowBudget {
    AdsorptionFlowBudget {
        absolute_amount_error: Amount::new(100.0),
        relative_amount_error: Dimensionless::new(0.0),
        maximum_step: Time::new(step),
        maximum_work: 5_000_000,
        ..AdsorptionFlowBudget::default()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
fn balances(result: &AdsorptionFlowAdvance, tolerance: f64) {
    assert!(result.dissolved.value() >= 0.0 && result.outlet_amount.value() >= 0.0);
    close(result.interval_balance_residual.value(), 0.0, tolerance);
    close(result.cumulative_balance_residual.value(), 0.0, tolerance);
    for (state, model) in result.populations.iter().zip(result.models()) {
        assert_eq!(state.id, model.id());
        assert!(state.adsorbed.value() >= 0.0 && state.adsorbed.value() <= model.capacity);
        assert!(state.vacant_capacity.value() >= 0.0);
        close(
            state.adsorbed.value() + state.vacant_capacity.value(),
            model.capacity,
            tolerance,
        );
        close(
            state.loading.value() * model.mass().value(),
            state.adsorbed.value(),
            tolerance,
        );
        close(
            state.cumulative_transfer.value(),
            state.adsorbed.value() - model.initial_adsorbed().value(),
            tolerance,
        );
    }
    close(
        result.dissolved.value()
            + result.outlet_amount.value()
            + result.transfers.iter().map(|(_, x)| x.value()).sum::<f64>(),
        result.initial_dissolved.value() + result.inlet_amount.value(),
        tolerance,
    );
}

#[test]
fn zero_rate_matches_independent_washout_and_outlet_integral_with_first_order_refinement() {
    let original = tank(vec![population("retained", 3.0, 2.0, 0.0, 1.0)], 8.0);
    // V=2, Q=.4, C_in=1: D(t)=2+6 exp(-.2t), O(t)=.4t+6(1-exp(-.2t)).
    let exact_dissolved = 2.0 + 6.0 * (-0.8_f64).exp();
    let exact_outlet = 1.6 + 6.0 * (1.0 - (-0.8_f64).exp());
    let mut errors = Vec::new();
    // Exact binary time partitions isolate convergence order from a final
    // sub-ulp clock remainder, which is tested separately below.
    for step in [0.5, 0.25, 0.125] {
        let mut model = original.clone();
        let result = advance(&mut model, 0.4, 1.0, 4.0, fixed_step(step));
        errors.push((result.dissolved.value() - exact_dissolved).abs());
        close(
            result.outlet_amount.value(),
            exact_outlet,
            errors.last().copied().unwrap() + 1e-11,
        );
        assert_eq!(result.populations[0].adsorbed.value(), 2.0);
        assert_eq!(result.transfers[0].1.value(), 0.0);
        balances(&result, 3e-12);
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.9 && pair[0] / pair[1] < 2.1,
            "{errors:?}"
        );
    }
    // Decimal timesteps must cover the requested endpoint without discarding
    // a small positive tail or trying to halve an unresolvable final duration.
    let mut decimal = original;
    let completed = advance(&mut decimal, 0.4, 1.0, 4.0, fixed_step(0.4));
    assert_eq!(completed.elapsed.value(), 4.0);
    assert_eq!(completed.inlet_amount.value(), 1.6);
    close(completed.dissolved.value(), exact_dissolved, 0.05);
    balances(&completed, 3e-12);
}

#[test]
fn stopped_flow_matches_independent_finite_bath_curve_and_existing_batch_owner() {
    // V=2, capacity=6, c_half=1, total=8: analytical ODE
    // dA/dt=.5*(4-A)*(12-A)/(10-A); its integral gives the time at A=2.
    let time = -(0.75 * 0.5_f64.ln() + 0.25 * (5.0_f64 / 6.0).ln()) / 0.5;
    let mut model = tank(vec![population("carbon", 3.0, 0.0, 0.5, 1.0)], 8.0);
    let budget = AdsorptionFlowBudget {
        absolute_amount_error: Amount::new(1e-9),
        relative_amount_error: Dimensionless::new(1e-8),
        maximum_work: 5_000_000,
        ..AdsorptionFlowBudget::default()
    };
    let result = advance(&mut model, 0.0, 99.0, time, budget);
    close(result.populations[0].adsorbed.value(), 2.0, 4e-4);
    assert_eq!(result.inlet_amount.value(), 0.0);
    assert_eq!(result.outlet_amount.value(), 0.0);
    assert!(result.rejected_steps > 0);
    assert!(result.maximum_local_error.value() <= 1e-9 + 8e-8);
    balances(&result, 5e-10);

    // An independently declared batch card is used ONLY by the batch owner.
    let mut batch_rate = rate_source(0.5);
    batch_rate.law = LawId(LDF_LAW.into());
    batch_rate.validity = batch_rate.validity.with_quantity(
        ELAPSED_TIME_AXIS,
        QuantitySpec::dimensional(Time::DIMS),
        0.0,
        1000.0,
    );
    let sorbent = SorbentPopulation::new(
        "carbon",
        LangmuirKinetics::from_card(equilibrium(1.0), batch_rate).unwrap(),
        Mass::new(3.0),
        Amount::new(0.0),
        conditions(),
    )
    .unwrap();
    let mut batch = AdsorptionBatch::new(
        vec![sorbent],
        Amount::new(8.0),
        Volume::new(2.0),
        Temperature::new(298.15),
    )
    .unwrap();
    batch
        .advance(
            Time::new(time),
            AdsorptionNumerics {
                absolute_amount_error: Amount::new(1e-9),
                relative_amount_error: 1e-8,
                max_work: 5_000_000,
                ..AdsorptionNumerics::default()
            },
            || false,
        )
        .unwrap();
    close(batch.dissolved().value(), result.dissolved.value(), 1e-9);
    close(
        batch.populations()[0].adsorbed.value(),
        result.populations[0].adsorbed.value(),
        1e-9,
    );
}

#[test]
fn stationary_loading_keeps_retained_carbon_and_independently_exports_the_feed() {
    let mut model = tank(
        vec![
            population("fast", 3.0, 4.0, 2.0, 1.0),
            population("slow", 1.0, 4.0 / 3.0, 0.1, 1.0),
        ],
        4.0,
    );
    let before = model.populations().unwrap();
    let result = advance(&mut model, 0.25, 2.0, 10.0, AdsorptionFlowBudget::default());
    assert_eq!(result.dissolved.value(), 4.0);
    assert_eq!(result.populations, before);
    assert_eq!(result.inlet_amount.value(), 5.0);
    assert_eq!(result.outlet_amount.value(), 5.0);
    assert!(result.transfers.iter().all(|(_, n)| n.value() == 0.0));
    assert_eq!(result.models(), model.models());
    balances(&result, 1e-12);
}

#[test]
fn shared_metal_and_loaded_carbon_survive_feed_changes_stops_and_restart() {
    let mut model = tank(
        vec![
            population("strong", 1.0, 0.0, 0.8, 0.05),
            population("weak", 1.0, 1.8, 0.3, 1.0),
        ],
        1.0,
    );
    let a = advance(&mut model, 0.0, 0.0, 1.0, AdsorptionFlowBudget::default());
    assert!(a.transfers[0].1.value() > 0.0);
    assert!(a.transfers[1].1.value() < 0.0);
    balances(&a, 1e-10);
    let b = advance(&mut model, 0.4, 3.0, 3.0, AdsorptionFlowBudget::default());
    assert!(b.outlet_amount.value() > 0.0);
    close(b.inlet_amount.value(), 3.6, 1e-14);
    balances(&b, 3e-10);
    let c = advance(&mut model, 0.0, 0.0, 2.0, AdsorptionFlowBudget::default());
    assert_eq!(c.outlet_amount.value(), 0.0);
    assert_eq!(c.cumulative_outlet, b.cumulative_outlet);
    balances(&c, 3e-10);
    let d = advance(&mut model, 1.0, 0.0, 12.0, AdsorptionFlowBudget::default());
    assert!(d.transfers.iter().all(|(_, n)| n.value() < 0.0));
    assert!(d.cumulative_outlet.value() > b.cumulative_outlet.value());
    assert_eq!(d.cumulative_inlet, b.cumulative_inlet);
    assert_eq!(d.previous_elapsed.value(), 6.0);
    assert_eq!(d.elapsed.value(), 18.0);
    assert_eq!(model.elapsed(), d.elapsed);
    balances(&d, 5e-10);
}

#[test]
fn canonical_order_subdivision_scaling_and_continuation_preserve_the_same_flow_problem() {
    let original = tank(
        vec![
            population("a", 3.0, 0.5, 0.5, 1.0),
            population("b", 2.0, 0.1, 0.2, 0.2),
        ],
        8.0,
    );
    let mut reverse = tank(original.models().iter().rev().cloned().collect(), 8.0);
    let mut once = original.clone();
    let numerics = fixed_step(0.125);
    let whole = advance(&mut once, 0.25, 1.0, 2.0, numerics);
    let reversed = advance(&mut reverse, 0.25, 1.0, 2.0, numerics);
    assert_eq!(once, reverse);
    assert_eq!(whole.outlet_amount, reversed.outlet_amount);
    let mut resumed = original.clone();
    let first = advance(&mut resumed, 0.25, 1.0, 1.0, numerics);
    let second = advance(&mut resumed, 0.25, 1.0, 1.0, numerics);
    assert_eq!(resumed.dissolved(), once.dissolved());
    assert_eq!(resumed.populations().unwrap(), once.populations().unwrap());
    close(
        first.outlet_amount.value() + second.outlet_amount.value(),
        whole.outlet_amount.value(),
        1e-13,
    );
    balances(&second, 1e-11);

    let mut split = tank(
        vec![
            population("a-1", 1.0, 0.5 / 3.0, 0.5, 1.0),
            population("a-2", 2.0, 1.0 / 3.0, 0.5, 1.0),
            population("b", 2.0, 0.1, 0.2, 0.2),
        ],
        8.0,
    );
    let divided = advance(&mut split, 0.25, 1.0, 2.0, numerics);
    close(divided.dissolved.value(), whole.dissolved.value(), 2e-11);
    close(
        divided.outlet_amount.value(),
        whole.outlet_amount.value(),
        2e-11,
    );
    close(
        divided.populations[0].adsorbed.value() + divided.populations[1].adsorbed.value(),
        whole.populations[0].adsorbed.value(),
        2e-11,
    );

    let mut scaled = AdsorptionCstr::new(
        vec![
            population("a", 12.0, 2.0, 0.5, 1.0),
            population("b", 8.0, 0.4, 0.2, 0.2),
        ],
        Amount::new(32.0),
        Volume::new(8.0),
        Temperature::new(298.15),
    )
    .unwrap();
    let larger = advance(&mut scaled, 1.0, 1.0, 2.0, numerics);
    close(
        larger.dissolved.value(),
        4.0 * whole.dissolved.value(),
        5e-11,
    );
    close(
        larger.outlet_amount.value(),
        4.0 * whole.outlet_amount.value(),
        5e-11,
    );
    balances(&larger, 1e-10);
}

#[test]
fn source_identity_dimensions_and_complete_physical_conditions_are_enforced() {
    let mut bad = rate_source(0.5);
    bad.law = LawId(LDF_LAW.into());
    assert!(LocalLangmuirTransfer::from_card(equilibrium(1.0), bad).is_err());
    for axis in [
        "elapsed-time",
        "residence-time",
        "mean-residence-time",
        "batch-volume",
    ] {
        let mut bad = rate_source(0.5);
        bad.validity =
            bad.validity
                .with_quantity(axis, QuantitySpec::dimensional(Time::DIMS), 0.0, 1000.0);
        assert!(LocalLangmuirTransfer::from_card(equilibrium(1.0), bad).is_err());
    }
    let mut bad = rate_source(0.5);
    bad.parameters.get_mut(RATE_PARAMETER).unwrap().dims = Amount::DIMS;
    assert!(LocalLangmuirTransfer::from_card(equilibrium(1.0), bad).is_err());
    let local = LocalLangmuirTransfer::from_card(equilibrium(1.0), rate_source(0.5)).unwrap();
    assert!(
        RetainedSorbent::new(
            "bad",
            local.clone(),
            Mass::new(1.0),
            Amount::new(3.0),
            conditions()
        )
        .is_err()
    );
    let forged = conditions()
        .with_quantity(CONCENTRATION_AXIS, concentration_spec(), 1.0)
        .unwrap();
    assert!(
        RetainedSorbent::new(
            "bad",
            local.clone(),
            Mass::new(1.0),
            Amount::new(0.0),
            forged
        )
        .is_err()
    );
    let incomplete = RetainedSorbent::new(
        "bad",
        local,
        Mass::new(1.0),
        Amount::new(0.0),
        QueryPoint::new(),
    )
    .unwrap();
    assert!(
        AdsorptionCstr::new(
            vec![incomplete],
            Amount::new(1.0),
            Volume::new(2.0),
            Temperature::new(298.15)
        )
        .is_err()
    );
    let mut model = tank(vec![population("carbon", 1.0, 0.0, 0.5, 1.0)], 1.0);
    let before = model.clone();
    for (flow, inlet, time) in [(-1.0, 0.0, 1.0), (0.1, f64::NAN, 1.0), (0.1, 1.0, -1.0)] {
        assert!(
            model
                .advance(
                    VolumetricFlowRate::new(flow),
                    AmountConcentration::new(inlet),
                    Time::new(time),
                    AdsorptionFlowBudget::default(),
                    || false
                )
                .is_err()
        );
        assert_eq!(model, before);
    }
}

#[test]
fn source_departure_work_exhaustion_and_cancellation_do_not_publish_partial_history() {
    let mut model = tank(vec![population("carbon", 3.0, 0.0, 0.5, 1.0)], 8.0);
    let before = model.clone();
    for budget in [
        AdsorptionFlowBudget {
            maximum_attempts: 1,
            maximum_step: Time::new(0.01),
            ..AdsorptionFlowBudget::default()
        },
        AdsorptionFlowBudget {
            maximum_root_iterations: 1,
            ..AdsorptionFlowBudget::default()
        },
        AdsorptionFlowBudget {
            maximum_work: 4,
            ..AdsorptionFlowBudget::default()
        },
    ] {
        assert_eq!(
            model
                .advance(
                    VolumetricFlowRate::new(0.1),
                    AmountConcentration::new(1.0),
                    Time::new(2.0),
                    budget,
                    || false
                )
                .unwrap_err(),
            AdsorptionFlowError::Budget
        );
        assert_eq!(model, before);
    }
    let mut calls = 0;
    assert_eq!(
        model
            .advance(
                VolumetricFlowRate::new(0.1),
                AmountConcentration::new(1.0),
                Time::new(2.0),
                AdsorptionFlowBudget::default(),
                || {
                    calls += 1;
                    calls > 15
                }
            )
            .unwrap_err(),
        AdsorptionFlowError::Cancelled
    );
    assert_eq!(model, before);
    assert!(matches!(
        model.advance(
            VolumetricFlowRate::new(10.0),
            AmountConcentration::new(1000.0),
            Time::new(10.0),
            AdsorptionFlowBudget::default(),
            || false
        ),
        Err(AdsorptionFlowError::OutsideConcentrationSupport)
    ));
    assert_eq!(model, before);
}

#[test]
fn zero_duration_empty_solution_and_positive_trace_stocks_keep_their_meaning() {
    let mut empty = tank(vec![population("carbon", 1.0, 0.0, 0.5, 1.0)], 0.0);
    let before = empty.clone();
    let result = advance(&mut empty, 10.0, 1.0, 0.0, AdsorptionFlowBudget::default());
    assert_eq!(empty, before);
    assert_eq!(result.inlet_amount.value(), 0.0);
    assert_eq!(result.outlet_amount.value(), 0.0);
    let hold = advance(&mut empty, 0.1, 0.0, 2.0, AdsorptionFlowBudget::default());
    assert_eq!(hold.dissolved.value(), 0.0);
    balances(&hold, 0.0);

    let mut trace = tank(vec![population("carbon", 1.0, 0.0, 0.0, 1.0)], 1e-200);
    let result = advance(
        &mut trace,
        0.4,
        0.0,
        1.0,
        AdsorptionFlowBudget {
            absolute_amount_error: Amount::new(1e-208),
            relative_amount_error: Dimensionless::new(1e-6),
            ..AdsorptionFlowBudget::default()
        },
    );
    assert!(result.dissolved.value() > 0.0 && result.outlet_amount.value() > 0.0);
    close(
        result.dissolved.value() + result.outlet_amount.value(),
        1e-200,
        1e-212,
    );
}
