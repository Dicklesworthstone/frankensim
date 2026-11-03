//! Independent analytic finite-bath reference and authored multi-sorbent fixtures.

use super::*;
use crate::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, AdsorptionInventory, CAPACITY_PARAMETER, HALF_SATURATION_PARAMETER,
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
        "aqueous-ldf-kinetics-test-v1",
        b"synthetic arithmetic fixture, not measured carbon kinetics",
    );
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: parameters
            .iter()
            .map(|&(key, value, dims)| (key.into(), LawParameter { value, dims }))
            .collect(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: "Authored arithmetic fixture; not calibrated gold/carbon behavior".into(),
            license: "project test fixture".into(),
            artifact: Some(artifact),
        },
    }
}
fn equilibrium_source() -> ConstitutiveModelCard {
    source(
        AQUEOUS_LANGMUIR_LAW,
        &[
            (CAPACITY_PARAMETER, 2.0, AdsorbentLoading::DIMS),
            (HALF_SATURATION_PARAMETER, 1.0, AmountConcentration::DIMS),
        ],
        ValidityDomain::unconstrained()
            .with_quantity(CONCENTRATION_AXIS, concentration_spec(), 0.0, 100.0)
            .with_quantity(TEMPERATURE_AXIS, temperature_spec(), 290.0, 310.0)
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0, 12.0),
    )
}
fn rate_source(rate: f64) -> ConstitutiveModelCard {
    let mut domain = equilibrium_source().validity;
    for (axis, dims, lo, hi) in [
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 1000.0),
        (LOADING_AXIS, AdsorbentLoading::DIMS, 0.0, 2.0),
        (VOLUME_AXIS, Volume::DIMS, 0.001, 1e6),
    ] {
        domain = domain.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi);
    }
    source(LDF_LAW, &[(RATE_PARAMETER, rate, Frequency::DIMS)], domain)
}
fn equilibrium() -> AqueousLangmuir {
    AqueousLangmuir::from_card(ElementId::new("Au").unwrap(), equilibrium_source()).unwrap()
}
fn model(rate: f64) -> LangmuirKinetics {
    LangmuirKinetics::from_card(equilibrium(), rate_source(rate)).unwrap()
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 11.0)
        .unwrap()
}
fn population(id: &str, mass: f64, adsorbed: f64, rate: f64) -> SorbentPopulation {
    SorbentPopulation::new(
        id,
        model(rate),
        Mass::new(mass),
        Amount::new(adsorbed),
        conditions(),
    )
    .unwrap()
}
fn batch(populations: Vec<SorbentPopulation>, dissolved: f64) -> AdsorptionBatch {
    AdsorptionBatch::new(
        populations,
        Amount::new(dissolved),
        Volume::new(2.0),
        Temperature::new(298.15),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance={tolerance:e}"
    );
}
fn assert_balances(batch: &AdsorptionBatch) {
    let total = batch.dissolved().value()
        + batch
            .populations()
            .iter()
            .map(|p| p.adsorbed.value())
            .sum::<f64>();
    close(
        total,
        batch.total_metal().value(),
        2e-12 * total.max(1e-300),
    );
    for (state, source) in batch.populations().iter().zip(batch.models()) {
        assert!(state.adsorbed.value() >= 0.0 && state.adsorbed.value() <= source.capacity);
        assert!(state.vacant_capacity.value() >= 0.0);
        close(
            state.adsorbed.value() + state.vacant_capacity.value(),
            source.capacity,
            1e-12 * source.capacity,
        );
    }
}

#[test]
fn independent_finite_bath_analytic_curve_and_equilibrium_limit() {
    // V=2, qmax*m=6, c_half=1, total=8. The loading ODE factors into
    // dx/dt = k (4-x)(12-x)/(10-x). Integrating from x=0 to x=2 gives:
    let time_for_two_mol = -(0.75 * 0.5_f64.ln() + 0.25 * (5.0_f64 / 6.0).ln()) / 0.5;
    let mut coarse = batch(vec![population("carbon", 3.0, 0.0, 0.5)], 8.0);
    let mut fine = coarse.clone();
    coarse
        .advance(
            Time::new(time_for_two_mol),
            AdsorptionNumerics {
                absolute_amount_error: Amount::new(1e-3),
                relative_amount_error: 1e-3,
                ..AdsorptionNumerics::default()
            },
            || false,
        )
        .unwrap();
    let answer = fine
        .advance(
            Time::new(time_for_two_mol),
            AdsorptionNumerics {
                absolute_amount_error: Amount::new(1e-9),
                relative_amount_error: 1e-8,
                max_work: 5_000_000,
                ..AdsorptionNumerics::default()
            },
            || false,
        )
        .unwrap();
    let coarse_error = (coarse.populations()[0].adsorbed.value() - 2.0).abs();
    let fine_error = (fine.populations()[0].adsorbed.value() - 2.0).abs();
    assert!(fine_error < coarse_error / 8.0);
    close(fine.populations()[0].adsorbed.value(), 2.0, 4e-4);
    assert!(answer.rejected_steps > 0 && answer.accepted_steps > 1);
    assert!(answer.maximum_local_error.value() <= 1e-9 + 8e-8);
    assert!(answer.accumulated_error_estimate.value() >= answer.maximum_local_error.value());
    assert_balances(&fine);
    let original = equilibrium()
        .equilibrate(
            AdsorptionInventory {
                dissolved: Amount::new(8.0),
                adsorbed: Amount::new(0.0),
            },
            Volume::new(2.0),
            Mass::new(3.0),
            Temperature::new(298.15),
            &conditions(),
        )
        .unwrap();
    let mut limit = batch(vec![population("carbon", 3.0, 0.0, 0.5)], 8.0);
    limit
        .advance(Time::new(80.0), AdsorptionNumerics::default(), || false)
        .unwrap();
    close(
        limit.dissolved().value(),
        original.final_inventory.dissolved.value(),
        2e-7,
    );
    close(
        limit.populations()[0].adsorbed.value(),
        original.final_inventory.adsorbed.value(),
        2e-7,
    );
    assert_balances(&limit);
}

#[test]
fn rates_and_time_scale_together_and_resumed_history_keeps_absolute_source_clock() {
    let mut slow = batch(vec![population("carbon", 3.0, 0.0, 0.5)], 8.0);
    let mut fast = batch(vec![population("carbon", 3.0, 0.0, 1.0)], 8.0);
    let numerics = AdsorptionNumerics {
        maximum_step: Time::new(0.125),
        ..AdsorptionNumerics::default()
    };
    let first = slow.advance(Time::new(1.0), numerics, || false).unwrap();
    fast.advance(
        Time::new(0.5),
        AdsorptionNumerics {
            maximum_step: Time::new(0.0625),
            ..numerics
        },
        || false,
    )
    .unwrap();
    assert_eq!(slow.dissolved(), fast.dissolved());
    assert_eq!(
        slow.populations()[0].adsorbed,
        fast.populations()[0].adsorbed
    );
    let mut uninterrupted = batch(vec![population("carbon", 3.0, 0.0, 0.5)], 8.0);
    uninterrupted
        .advance(Time::new(2.0), numerics, || false)
        .unwrap();
    let second = slow.advance(Time::new(1.0), numerics, || false).unwrap();
    close(
        slow.dissolved().value(),
        uninterrupted.dissolved().value(),
        2e-4,
    );
    assert_eq!(second.previous_elapsed, Time::new(1.0));
    assert_eq!(slow.elapsed(), Time::new(2.0));
    assert_eq!(second.initial_dissolved, first.dissolved);
    assert_eq!(second.initial_points[0].axes()[ELAPSED_TIME_AXIS], 1.0);
    assert_eq!(second.final_points[0].axes()[ELAPSED_TIME_AXIS], 2.0);
    assert_eq!(second.models()[0].model().equilibrium(), &equilibrium());
    assert_eq!(second.models()[0].model().source(), &rate_source(0.5));
    assert_eq!(second.models()[0].initially_adsorbed(), Amount::new(0.0));
    assert!(second.accumulated_error_estimate.value() >= first.accumulated_error_estimate.value());
    assert_balances(&slow);
}

#[test]
fn shared_solution_competition_preloading_and_reversible_desorption() {
    let commercial = population("activated", 1.0, 0.0, 1.0);
    let natural = population("natural", 1.0, 0.0, 0.4);
    let mut alone = batch(vec![commercial.clone()], 2.0);
    let mut competing = batch(vec![natural, commercial], 2.0);
    let single = alone
        .advance(Time::new(2.0), AdsorptionNumerics::default(), || false)
        .unwrap();
    let shared = competing
        .advance(Time::new(2.0), AdsorptionNumerics::default(), || false)
        .unwrap();
    assert!(shared.populations[0].adsorbed.value() < single.populations[0].adsorbed.value());
    assert!(shared.populations[1].adsorbed.value() > 0.0);
    assert_eq!(
        shared.final_points[0].axes()[CONCENTRATION_AXIS],
        shared.final_points[1].axes()[CONCENTRATION_AXIS]
    );
    assert_balances(&competing);
    let mut wash = batch(
        vec![
            population("activated", 1.0, 0.0, 1.0),
            population("preloaded", 1.0, 1.5, 0.4),
        ],
        0.0,
    );
    let washed = wash
        .advance(Time::new(2.0), AdsorptionNumerics::default(), || false)
        .unwrap();
    assert!(washed.dissolved.value() > 0.0);
    assert!(washed.transfers[0].1.value() > 0.0);
    assert!(washed.transfers[1].1.value() < 0.0);
    assert!(washed.populations[1].adsorbed.value() < 1.5);
    assert_balances(&wash);
}

#[test]
fn canonical_order_population_subdivision_and_extensive_scaling_preserve_physics() {
    let a = population("a", 1.0, 0.2, 0.5);
    let z = population("z", 2.0, 0.1, 0.3);
    let mut first = batch(vec![z.clone(), a.clone()], 8.0);
    let mut permuted = batch(vec![a, z], 8.0);
    first
        .advance(Time::new(0.75), AdsorptionNumerics::default(), || false)
        .unwrap();
    permuted
        .advance(Time::new(0.75), AdsorptionNumerics::default(), || false)
        .unwrap();
    assert_eq!(first, permuted);
    let mut whole = batch(vec![population("whole", 3.0, 0.0, 0.5)], 8.0);
    let mut split = batch(
        vec![
            population("half-a", 1.5, 0.0, 0.5),
            population("half-b", 1.5, 0.0, 0.5),
        ],
        8.0,
    );
    let mut scaled = AdsorptionBatch::new(
        vec![population("whole", 30.0, 0.0, 0.5)],
        Amount::new(80.0),
        Volume::new(20.0),
        Temperature::new(298.15),
    )
    .unwrap();
    whole
        .advance(Time::new(0.75), AdsorptionNumerics::default(), || false)
        .unwrap();
    split
        .advance(Time::new(0.75), AdsorptionNumerics::default(), || false)
        .unwrap();
    scaled
        .advance(
            Time::new(0.75),
            AdsorptionNumerics {
                absolute_amount_error: Amount::new(1e-7),
                ..AdsorptionNumerics::default()
            },
            || false,
        )
        .unwrap();
    close(split.dissolved().value(), whole.dissolved().value(), 2e-9);
    close(
        split.populations()[0].adsorbed.value() * 2.0,
        whole.populations()[0].adsorbed.value(),
        2e-9,
    );
    close(
        scaled.dissolved().value() / 10.0,
        whole.dissolved().value(),
        2e-9,
    );
    assert_balances(&split);
    assert_balances(&scaled);
}

#[test]
fn zero_time_zero_metal_and_explicit_zero_rates_do_not_invent_transfer() {
    let mut empty = batch(vec![population("carbon", 3.0, 0.0, 0.5)], 0.0);
    let answer = empty
        .advance(Time::new(4.0), AdsorptionNumerics::default(), || false)
        .unwrap();
    assert_eq!(answer.dissolved.value(), 0.0);
    assert_eq!(answer.populations[0].adsorbed.value(), 0.0);
    assert_eq!(answer.populations[0].vacant_capacity.value(), 6.0);
    let mut disabled = batch(vec![population("carbon", 3.0, 2.0, 0.0)], 8.0);
    let original = disabled.clone();
    let none = disabled
        .advance(Time::new(0.0), AdsorptionNumerics::default(), || false)
        .unwrap();
    assert_eq!(disabled, original);
    assert_eq!(none.accepted_steps, 0);
    disabled
        .advance(Time::new(4.0), AdsorptionNumerics::default(), || false)
        .unwrap();
    assert_eq!(disabled.dissolved(), original.dissolved());
    assert_eq!(disabled.populations(), original.populations());
    assert_eq!(disabled.elapsed(), Time::new(4.0));
    assert_balances(&disabled);
}

#[test]
fn sources_layouts_capacities_owned_coordinates_and_numerical_inputs_are_checked() {
    for rate in [-0.1, f64::INFINITY, f64::NAN] {
        assert!(LangmuirKinetics::from_card(equilibrium(), rate_source(rate)).is_err());
    }
    let mut wrong = rate_source(0.5);
    wrong.parameters.get_mut(RATE_PARAMETER).unwrap().dims = Time::DIMS;
    assert!(LangmuirKinetics::from_card(equilibrium(), wrong).is_err());
    let mut unsupported = rate_source(0.5);
    unsupported.validity = equilibrium_source().validity;
    assert!(LangmuirKinetics::from_card(equilibrium(), unsupported).is_err());
    let mut unsourced = rate_source(0.5);
    unsourced.sources.clear();
    assert!(LangmuirKinetics::from_card(equilibrium(), unsourced).is_err());
    for (mass, loading) in [(0.0, 0.0), (-1.0, 0.0), (1.0, 2.1), (1.0, -0.1)] {
        assert!(
            SorbentPopulation::new(
                "carbon",
                model(0.5),
                Mass::new(mass),
                Amount::new(loading),
                conditions()
            )
            .is_err()
        );
    }
    let owned = conditions()
        .with_quantity(
            ELAPSED_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.0,
        )
        .unwrap();
    assert!(
        SorbentPopulation::new(
            "carbon",
            model(0.5),
            Mass::new(1.0),
            Amount::new(0.0),
            owned
        )
        .is_err()
    );
    assert!(
        AdsorptionBatch::new(
            vec![],
            Amount::new(1.0),
            Volume::new(1.0),
            Temperature::new(298.15)
        )
        .is_err()
    );
    let duplicate = population("same", 1.0, 0.0, 0.5);
    assert!(
        AdsorptionBatch::new(
            vec![duplicate.clone(), duplicate],
            Amount::new(1.0),
            Volume::new(1.0),
            Temperature::new(298.15)
        )
        .is_err()
    );
    let copper = LangmuirKinetics::from_card(
        AqueousLangmuir::from_card(ElementId::new("Cu").unwrap(), equilibrium_source()).unwrap(),
        rate_source(0.5),
    )
    .unwrap();
    let copper = SorbentPopulation::new(
        "copper",
        copper,
        Mass::new(1.0),
        Amount::new(0.0),
        conditions(),
    )
    .unwrap();
    assert!(
        AdsorptionBatch::new(
            vec![population("gold", 1.0, 0.0, 0.5), copper],
            Amount::new(1.0),
            Volume::new(1.0),
            Temperature::new(298.15)
        )
        .is_err()
    );
    let mut bath = batch(vec![population("carbon", 3.0, 0.0, 0.5)], 8.0);
    let before = bath.clone();
    for duration in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            bath.advance(Time::new(duration), AdsorptionNumerics::default(), || false)
                .is_err()
        );
        assert_eq!(bath, before);
    }
    assert!(
        bath.advance(
            Time::new(1.0),
            AdsorptionNumerics {
                relative_amount_error: 0.0,
                ..AdsorptionNumerics::default()
            },
            || false
        )
        .is_err()
    );
    assert_eq!(bath, before);
}

#[test]
fn unsupported_trajectory_cancellation_and_work_failures_are_atomic() {
    let mut narrow = equilibrium_source();
    narrow.validity =
        narrow
            .validity
            .with_quantity(CONCENTRATION_AXIS, concentration_spec(), 3.0, 100.0);
    let model = LangmuirKinetics::from_card(
        AqueousLangmuir::from_card(ElementId::new("Au").unwrap(), narrow).unwrap(),
        rate_source(0.5),
    )
    .unwrap();
    let mut bath = batch(
        vec![
            SorbentPopulation::new(
                "carbon",
                model,
                Mass::new(3.0),
                Amount::new(0.0),
                conditions(),
            )
            .unwrap(),
        ],
        8.0,
    );
    let before = bath.clone();
    assert!(matches!(
        bath.advance(Time::new(20.0), AdsorptionNumerics::default(), || false),
        Err(AdsorptionKineticsError::OutsideConcentrationSupport)
    ));
    assert_eq!(bath, before);
    let mut bath = batch(vec![population("carbon", 3.0, 0.0, 0.5)], 8.0);
    let before = bath.clone();
    assert!(matches!(
        bath.advance(Time::new(1001.0), AdsorptionNumerics::default(), || false),
        Err(AdsorptionKineticsError::Source(
            AqueousAdsorptionError::OutsideValidity { .. }
        ))
    ));
    assert_eq!(bath, before);
    let mut polls = 0;
    assert_eq!(
        bath.advance(
            Time::new(1.0),
            AdsorptionNumerics {
                maximum_step: Time::new(0.01),
                ..AdsorptionNumerics::default()
            },
            || {
                polls += 1;
                polls > 2000
            }
        )
        .unwrap_err(),
        AdsorptionKineticsError::Cancelled
    );
    assert_eq!(bath, before);
    assert_eq!(
        bath.advance(
            Time::new(1.0),
            AdsorptionNumerics {
                max_work: 1,
                ..AdsorptionNumerics::default()
            },
            || false
        )
        .unwrap_err(),
        AdsorptionKineticsError::Budget
    );
    assert_eq!(bath, before);
    assert_eq!(
        bath.advance(
            Time::new(1.0),
            AdsorptionNumerics {
                max_root_iterations: 1,
                ..AdsorptionNumerics::default()
            },
            || false
        )
        .unwrap_err(),
        AdsorptionKineticsError::Budget
    );
    assert_eq!(bath, before);
}

#[test]
fn positive_range_loss_refuses_and_trace_inventory_remains_positive() {
    assert!(
        SorbentPopulation::new(
            "overflow",
            model(0.5),
            Mass::new(f64::MAX),
            Amount::new(0.0),
            conditions()
        )
        .is_err()
    );
    let mut bath = batch(
        vec![population("tiny-rate", 3.0, 0.0, f64::MIN_POSITIVE)],
        8.0,
    );
    let before = bath.clone();
    assert!(
        bath.advance(
            Time::new(f64::MIN_POSITIVE),
            AdsorptionNumerics::default(),
            || false
        )
        .is_err()
    );
    assert_eq!(bath, before);
    let mut trace = batch(vec![population("carbon", 3.0, 0.0, 0.5)], 1e-300);
    trace
        .advance(Time::new(0.1), AdsorptionNumerics::default(), || false)
        .unwrap();
    assert!(trace.dissolved().value() > 0.0 && trace.populations()[0].adsorbed.value() > 0.0);
    assert_balances(&trace);
}
