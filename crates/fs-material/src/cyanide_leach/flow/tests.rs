//! Manufactured local-rate fixtures, not measured ore or plant performance.
use super::*;
use crate::cyanide_leach::tests as fixture;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};

fn source(rate: f64) -> ConstitutiveModelCard {
    let artifact = hash_domain(
        "authored-local-cyanidation-cstr-v1",
        b"no measured kinetic calibration",
    );
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims) in local_axes() {
        let (lo, hi) = match axis {
            PH_AXIS => (-1.0, 15.0),
            TEMPERATURE_AXIS => (298.15, 298.15),
            VOLUME_AXIS | ACTIVITY_REFERENCE_AXIS => (0.001, 10.0),
            _ => (0.0, 10.0),
        };
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi);
    }
    ConstitutiveModelCard {
        law: LawId(LOCAL_CYANIDATION_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [(
            RATE_PARAMETER.into(),
            LawParameter {
                value: rate,
                dims: CyanidationRateConstant::DIMS,
            },
        )]
        .into(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: "authored local kinetic fixture; no ore calibration".into(),
            license: "project test fixture".into(),
            artifact: Some(artifact),
        },
    }
}
pub(super) fn model(rate: f64, buffer: bool) -> LocalGoldCyanidation {
    LocalGoldCyanidation::from_card(fixture::system(buffer), "uncomplexed-cyanide", source(rate))
        .unwrap()
}
pub(super) fn numerics() -> FlowNumerics {
    FlowNumerics {
        maximum_step: Time::new(0.125),
        absolute_amount_error: Amount::new(1e-9),
        relative_amount_error: Dimensionless::new(2e-6),
        maximum_work: 10_000_000,
        ..FlowNumerics::default()
    }
}
pub(super) fn feed(charge: &CyanidationCharge) -> CyanidationFeed {
    let concentration =
        |n: Amount| AmountConcentration::new(n.value() / charge.solution_volume.value());
    CyanidationFeed {
        accessible_gold: concentration(charge.accessible_gold),
        cyanide_ligand: concentration(charge.uncomplexed_cyanide),
        dissolved_oxygen: concentration(charge.dissolved_oxygen),
        gold_complex: concentration(charge.gold_complex),
        other_families: charge
            .other_families
            .iter()
            .map(|(id, n)| (id.clone(), concentration(*n)))
            .collect(),
        spectators: charge
            .spectators
            .iter()
            .map(|ion| SpectatorFeed {
                species: ion.species.clone(),
                charge: ion.charge,
                concentration: concentration(ion.amount),
            })
            .collect(),
        competing_metals: vec![],
    }
}
fn metals(total: f64) -> competition::CompetingMetals {
    let batch =
        fixture::competition::competing_batch(fixture::competition::copper_charge(total), total);
    let original = batch.layout.competition.as_ref().unwrap();
    competition::CompetingMetals {
        families: original.model.families().to_vec(),
        totals: original.totals.clone(),
        numerics: original.numerics,
    }
}
pub(super) fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance={tolerance:e}"
    );
}
fn vector(value: &CyanidationInventory) -> Vec<f64> {
    let mut values = vec![
        value.accessible_gold.value(),
        value.cyanide_ligand.value(),
        value.dissolved_oxygen.value(),
        value.gold_complex.value(),
    ];
    values.extend(value.other_families.iter().map(|(_, n)| n.value()));
    values.extend(value.spectators.iter().map(|s| s.amount.value()));
    values.extend(value.competing_metals.iter().map(|(_, n)| n.value()));
    values
}
fn balances(answer: &CyanidationFlowAdvance) {
    for residual in [
        &answer.interval_conservation,
        &answer.cumulative_conservation,
    ] {
        for amount in vector(&residual.components) {
            close(amount, 0.0, 2e-12);
        }
        close(residual.gold.value(), 0.0, 2e-12);
        close(residual.sodium.value(), 0.0, 2e-12);
        close(residual.dissolved_oxygen.value(), 0.0, 2e-12);
        close(
            residual.cyanide.value(),
            answer
                .state
                .competition
                .as_ref()
                .map_or(0.0, |s| s.ligand_balance_residual.value()),
            2e-12,
        );
    }
}
fn unchanged(current: &CyanidationCstr, before: &CyanidationCstr) {
    assert_eq!(current.raw, before.raw);
    assert_eq!(current.inlet, before.inlet);
    assert_eq!(current.outlet, before.outlet);
    assert_eq!(current.gold_dissolved, before.gold_dissolved);
    assert_eq!(current.accumulated_error, before.accumulated_error);
    assert_eq!(current.state.elapsed, before.state.elapsed);
    assert_eq!(current.state.conditions, before.state.conditions);
    assert_eq!(current.state.free_cyanide, before.state.free_cyanide);
    assert_eq!(current.initial_charge(), before.initial_charge());
}

#[test]
fn zero_reaction_exactly_transports_all_family_spectator_and_metal_pools() {
    let mut charge = fixture::competition::copper_charge(0.2);
    charge.gold_complex = Amount::new(0.1);
    charge.other_families = vec![("carbonate-buffer".into(), Amount::new(0.4))];
    let metal = metals(0.2);
    let metal_id = metal.totals[0].0.clone();
    let mut inlet = feed(&charge);
    inlet.accessible_gold = AmountConcentration::new(0.3);
    inlet.cyanide_ligand = AmountConcentration::new(1.0);
    inlet.dissolved_oxygen = AmountConcentration::new(0.1);
    inlet.gold_complex = AmountConcentration::new(0.05);
    inlet.other_families[0].1 = AmountConcentration::new(0.1);
    inlet.competing_metals = vec![(metal_id, AmountConcentration::new(0.1))];
    for ion in &mut inlet.spectators {
        ion.concentration = AmountConcentration::new(if ion.charge == 1 { 0.6 } else { 0.1 });
    }
    let mut tank = CyanidationCstr::new_with_competing_metals(
        model(0.0, true),
        charge.clone(),
        QueryPoint::new(),
        metal,
    )
    .unwrap();
    let original = tank.raw.clone();
    let answer = tank
        .advance(
            Time::new(1.2),
            VolumetricFlowRate::new(0.7),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    let exposure = 0.84_f64;
    let e = (-exposure).exp();
    for (((n, exported), supplied), initial) in vector(&answer.state.inventory)
        .into_iter()
        .zip(vector(&answer.outlet))
        .zip(vector(&answer.inlet))
        .zip(original)
    {
        let target = supplied / exposure;
        close(n, target + (initial - target) * e, 2e-13);
        // Independent integral of the continuously mixed outlet concentration.
        close(
            exported,
            target * exposure + (initial - target) * (1.0 - e),
            2e-13,
        );
    }
    assert_eq!(answer.interval_gold_dissolved.value(), 0.0);
    assert_eq!(answer.accepted_steps, 1);
    assert_eq!(answer.model().source(), &source(0.0));
    assert_eq!(answer.initial_charge(), &charge);
    assert!(
        answer
            .state
            .competition
            .as_ref()
            .unwrap()
            .bound_ligand
            .value()
            > 0.0
    );
    balances(&answer);
}

#[test]
fn stopped_flow_matches_the_original_batch_and_its_independent_known_answer() {
    let charge = fixture::charge(1.0);
    let inlet = feed(&charge);
    let mut tank = CyanidationCstr::new(
        model(4.0, false),
        charge.clone(),
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let mut batch = GoldCyanidationBatch::new(
        GoldCyanidation::from_card(
            fixture::system(false),
            "uncomplexed-cyanide",
            fixture::card(4.0),
        )
        .unwrap(),
        charge,
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let expected = batch
        .advance(Time::new(0.5), CyanidationNumerics::default(), || false)
        .unwrap();
    let answer = tank
        .advance(
            Time::new(0.5),
            VolumetricFlowRate::new(0.0),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    // Authored Ka=sqrt(Kw), CN=2Au and O2=Au/4 give dAu/dt=-Au^3.
    close(
        answer.state.inventory.accessible_gold.value(),
        1.0 / 2.0_f64.sqrt(),
        2e-6,
    );
    close(
        answer.cumulative_gold_dissolved.value(),
        expected.state.gold_dissolved.value(),
        2e-6,
    );
    close(
        answer.state.free_cyanide.value(),
        expected.state.free_cyanide.value(),
        2e-6,
    );
    close(
        answer.state.speciation.ph.value(),
        expected.state.speciation.ph.value(),
        1e-10,
    );
    assert!(vector(&answer.outlet).iter().all(|n| *n == 0.0));
    assert!(vector(&answer.inlet).iter().all(|n| *n == 0.0));
    balances(&answer);
}

#[test]
fn simultaneous_flow_and_reaction_keep_exports_distinct_from_reaction_progress() {
    let charge = fixture::charge(1.0);
    let inlet = feed(&charge);
    let mut tank = CyanidationCstr::new(
        model(4.0, false),
        charge.clone(),
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let answer = tank
        .advance(
            Time::new(0.5),
            VolumetricFlowRate::new(0.8),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    let x = answer.cumulative_gold_dissolved.value();
    assert!(x > 0.0 && answer.outlet.gold_complex.value() > 0.0);
    assert!(answer.outlet.accessible_gold.value() > 0.0);
    assert!(answer.state.inventory.gold_complex.value() < x);
    close(
        answer.state.inventory.gold_complex.value() + answer.outlet.gold_complex.value(),
        x,
        2e-13,
    );
    close(
        answer.state.inventory.accessible_gold.value() + answer.outlet.accessible_gold.value() + x,
        1.0 + answer.inlet.accessible_gold.value(),
        2e-13,
    );
    close(
        answer.state.inventory.cyanide_ligand.value()
            + answer.outlet.cyanide_ligand.value()
            + 2.0 * x,
        2.0 + answer.inlet.cyanide_ligand.value(),
        4e-13,
    );
    close(
        answer.state.inventory.dissolved_oxygen.value()
            + answer.outlet.dissolved_oxygen.value()
            + x / 4.0,
        0.25 + answer.inlet.dissolved_oxygen.value(),
        1e-13,
    );
    assert_eq!(
        answer.state.conditions.axes()[FREE_CYANIDE_AXIS],
        answer.state.free_cyanide.value()
    );
    assert!(
        !answer
            .state
            .conditions
            .axes()
            .contains_key(ELAPSED_TIME_AXIS)
    );
    balances(&answer);
}

#[test]
fn flowing_copper_repartitions_actual_free_cyanide_and_retains_its_source_equilibrium() {
    let charge = fixture::competition::copper_charge(0.2);
    let metal = metals(0.2);
    let id = metal.totals[0].0.clone();
    let mut inlet = feed(&charge);
    inlet.competing_metals = vec![(id, AmountConcentration::new(0.1))];
    let mut tank = CyanidationCstr::new_with_competing_metals(
        model(4.0, false),
        charge,
        QueryPoint::new(),
        metal,
    )
    .unwrap();
    let before = tank.state.free_cyanide;
    let answer = tank
        .advance(
            Time::new(0.125),
            VolumetricFlowRate::new(0.5),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    let competition = answer.state.competition.as_ref().unwrap();
    let expected_copper = 0.1 + 0.1 * (-0.0625_f64).exp();
    close(competition.metals[0].total.value(), expected_copper, 2e-13);
    close(
        competition.metals[0]
            .members
            .iter()
            .map(|m| m.amount.value())
            .sum(),
        expected_copper,
        2e-12,
    );
    close(
        competition.ligand_total.value(),
        answer.state.inventory.cyanide_ligand.value(),
        0.0,
    );
    close(
        answer.state.uncomplexed_cyanide.value() + competition.bound_ligand.value(),
        answer.state.inventory.cyanide_ligand.value(),
        3e-9,
    );
    assert_ne!(answer.state.free_cyanide, before);
    assert!(competition.bound_ligand.value() > 0.0);
    close(
        answer.state.rate.value(),
        4.0 * answer.state.inventory.accessible_gold.value()
            * answer.state.free_cyanide.value()
            * answer.state.inventory.dissolved_oxygen.value(),
        2e-13,
    );
    let activity = |name: &str| {
        competition.metals[0]
            .members
            .iter()
            .find(|m| m.species.as_str() == name)
            .unwrap()
            .activity
            .value()
    };
    close(
        activity("CuCN(aq)") / (activity("Cu[+](aq)") * answer.state.free_cyanide.value()),
        10.0,
        3e-8,
    );
    close(
        activity("Cu(CN)2[-](aq)") / (activity("CuCN(aq)") * answer.state.free_cyanide.value()),
        5.0,
        3e-8,
    );
    balances(&answer);
}

#[test]
fn feed_changes_stop_restart_and_partitioned_calls_preserve_complete_history() {
    let charge = fixture::charge(1.0);
    let mut inlet = feed(&charge);
    let mut single = CyanidationCstr::new(
        model(4.0, false),
        charge.clone(),
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let mut continued = single.clone();
    let one = single
        .advance(
            Time::new(0.5),
            VolumetricFlowRate::new(0.8),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    continued
        .advance(
            Time::new(0.25),
            VolumetricFlowRate::new(0.8),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    let two = continued
        .advance(
            Time::new(0.25),
            VolumetricFlowRate::new(0.8),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    for (a, b) in vector(&one.state.inventory)
        .into_iter()
        .zip(vector(&two.state.inventory))
    {
        close(a, b, 3e-6);
    }
    for (a, b) in vector(&one.cumulative_outlet)
        .into_iter()
        .zip(vector(&two.cumulative_outlet))
    {
        close(a, b, 3e-6);
    }
    let original_in = continued.cumulative_inlet();
    let original_out = continued.cumulative_outlet();
    let original_progress = continued.cumulative_gold_dissolved();
    let stopped = continued
        .advance(
            Time::new(0.125),
            VolumetricFlowRate::new(0.0),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    assert_eq!(continued.cumulative_inlet(), original_in);
    assert_eq!(continued.cumulative_outlet(), original_out);
    assert!(continued.cumulative_gold_dissolved().value() > original_progress.value());
    balances(&stopped);
    inlet.cyanide_ligand = AmountConcentration::new(3.0);
    inlet.spectators[0].concentration = AmountConcentration::new(1.5);
    let restarted = continued
        .advance(
            Time::new(0.125),
            VolumetricFlowRate::new(0.5),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    assert_eq!(restarted.previous_elapsed, Time::new(0.625));
    assert_eq!(restarted.elapsed, Time::new(0.75));
    assert_eq!(continued.initial_charge(), &charge);
    close(restarted.inlet.cyanide_ligand.value(), 0.1875, 0.0);
    balances(&restarted);
}

#[test]
fn batch_card_feed_source_budget_and_cancellation_refusals_are_atomic() {
    assert!(
        LocalGoldCyanidation::from_card(
            fixture::system(false),
            "uncomplexed-cyanide",
            fixture::card(4.0)
        )
        .is_err()
    );
    let mut aged = source(4.0);
    aged.validity = aged.validity.with_quantity(
        ELAPSED_TIME_AXIS,
        QuantitySpec::dimensional(Time::DIMS),
        0.0,
        10.0,
    );
    assert!(
        LocalGoldCyanidation::from_card(fixture::system(false), "uncomplexed-cyanide", aged)
            .is_err()
    );
    let charge = fixture::charge(1.0);
    let inlet = feed(&charge);
    let mut tank = CyanidationCstr::new(
        model(4.0, false),
        charge,
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let before = tank.clone();
    for kind in 0..4 {
        let mut bad = inlet.clone();
        match kind {
            0 => bad.spectators.clear(),
            1 => bad.spectators[0].charge = 2,
            2 => bad.spectators.push(bad.spectators[0].clone()),
            _ => bad.cyanide_ligand = AmountConcentration::new(f64::NAN),
        }
        assert!(
            tank.advance(
                Time::new(0.125),
                VolumetricFlowRate::new(0.5),
                &bad,
                numerics(),
                || false
            )
            .is_err()
        );
        unchanged(&tank, &before);
    }
    assert!(matches!(
        tank.advance(
            Time::new(0.125),
            VolumetricFlowRate::new(0.5),
            &inlet,
            FlowNumerics {
                maximum_work: 1,
                ..numerics()
            },
            || false
        ),
        Err(CyanidationError::Budget)
    ));
    unchanged(&tank, &before);
    let mut probe = tank.clone();
    let required = probe
        .advance(
            Time::new(0.125),
            VolumetricFlowRate::new(0.5),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap()
        .work;
    let mut polls = 0;
    assert!(matches!(
        tank.advance(
            Time::new(0.125),
            VolumetricFlowRate::new(0.5),
            &inlet,
            numerics(),
            || {
                polls += 1;
                polls >= required
            }
        ),
        Err(CyanidationError::Cancelled)
    ));
    unchanged(&tank, &before);
    let mut narrow = source(4.0);
    narrow.validity =
        narrow
            .validity
            .with_quantity(PH_AXIS, QuantitySpec::dimensional(Dims::NONE), 6.99, 7.01);
    let mut limited = CyanidationCstr::new(
        LocalGoldCyanidation::from_card(fixture::system(false), "uncomplexed-cyanide", narrow)
            .unwrap(),
        fixture::charge(1.0),
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let before = limited.clone();
    let mut unsupported = inlet.clone();
    unsupported.spectators[0].concentration = AmountConcentration::new(0.0);
    assert!(
        limited
            .advance(
                Time::new(0.25),
                VolumetricFlowRate::new(1.0),
                &unsupported,
                numerics(),
                || false
            )
            .is_err()
    );
    unchanged(&limited, &before);
}

#[test]
fn resolved_gross_transport_accepts_rounding_at_near_equilibrium() {
    let mut charge = fixture::charge(1.0);
    charge.accessible_gold = Amount::new(0.2);
    charge.solution_volume = Volume::new(1.0);
    for target in [0.2_f64.next_down(), 0.2_f64.next_up()] {
        let mut tank = CyanidationCstr::new(
            model(0.0, false),
            charge.clone(),
            QueryPoint::new(),
            fixture::options(),
        )
        .unwrap();
        let mut inlet = feed(&charge);
        inlet.accessible_gold = AmountConcentration::new(target);
        let report = tank
            .advance(
                Time::new(0.25),
                VolumetricFlowRate::new(0.5),
                &inlet,
                numerics(),
                || false,
            )
            .unwrap();
        // The exact net change is under half an ulp. Real streams still cross
        // the boundary and their independent analytical integral is retained.
        assert_eq!(
            report.state.inventory.accessible_gold,
            charge.accessible_gold
        );
        assert!(report.inlet.accessible_gold.value() > 0.0);
        assert!(report.outlet.accessible_gold.value() > 0.0);
        let x = 0.125_f64;
        let exported = target * x + (0.2 - target) * -(-x).exp_m1();
        close(report.outlet.accessible_gold.value(), exported, 2e-17);
        balances(&report);
        assert_eq!(report.interval_gold_dissolved, Amount::new(0.0));
    }
}
