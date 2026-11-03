//! G0/G1 manufactured solutions for the coupled flow, gas and reaction owner.
//! The gas fixture has independently authored K=1; no slurry calibration claim.
use super::*;
use crate::cyanide_leach::tests::gas_transfer::{henry_context, measured_source};
use crate::cyanide_leach::{
    aeration::{AeratedGoldCyanidationBatch, AerationNumerics, KLA_PARAMETER, OXYGEN_TRANSFER_LAW},
    tests as fixture,
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::Frequency;
use fs_thermochem::stoichiometry::thermodynamics::{
    FormationGibbsSource, gas_liquid::IdealGasLiquidPartition,
};

fn local(rate: f64) -> LocalGoldCyanidation {
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
    let artifact = hash_domain("local-aerated-cstr-known-answer-v1", b"authored local k");
    let source = ConstitutiveModelCard {
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
            source: "authored local kinetics; no ore calibration".into(),
            license: "project test fixture".into(),
            artifact: Some(artifact),
        },
    };
    LocalGoldCyanidation::from_card(fixture::system(false), "uncomplexed-cyanide", source).unwrap()
}
fn solubility(reference_pressure: f64) -> IdealGasLiquidPartition {
    let original = fixture::point(OXYGEN);
    let mut gas_at = QueryPoint::new();
    let mut liquid_at = QueryPoint::new();
    for (axis, value) in original.axes() {
        let value = if axis == "pressure" {
            reference_pressure
        } else {
            *value
        };
        liquid_at = liquid_at
            .with_quantity(axis, original.axis_quantities()[axis], value)
            .unwrap();
        if !matches!(
            axis.as_str(),
            "source-phase-aqueous" | "source-aqueous-unit-activity-infinite-dilution"
        ) {
            gas_at = gas_at
                .with_quantity(axis, original.axis_quantities()[axis], value)
                .unwrap();
        }
    }
    for axis in ["source-phase-gas", "source-pure-species"] {
        gas_at = gas_at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    let gas_pack = fixture::pack("O2(g)", 0.0, &gas_at);
    let liquid_pack = fixture::pack(OXYGEN, 0.0, &liquid_at);
    let gas = SpeciesId::new("O2(g)").unwrap();
    let liquid = SpeciesId::new(OXYGEN).unwrap();
    IdealGasLiquidPartition::from_sources(
        MineralReaction::OxygenDissolution.build().unwrap(),
        gas.clone(),
        liquid.clone(),
        &[
            FormationGibbsSource {
                species: &gas,
                pack: &gas_pack,
                at: &gas_at,
                reference_axis: fixture::REFERENCE,
            },
            FormationGibbsSource {
                species: &liquid,
                pack: &liquid_pack,
                at: &liquid_at,
                reference_axis: fixture::REFERENCE,
            },
        ],
    )
    .unwrap()
}
fn transfer_card(kla: f64) -> ConstitutiveModelCard {
    let artifact = hash_domain("aerated-cstr-transfer-known-answer-v1", b"authored kLa");
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, low, high) in [
        (OXYGEN_PARTIAL_PRESSURE_AXIS, Pressure::DIMS, 0.0, 200_000.0),
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 10.0),
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 100.0),
        (TEMPERATURE_AXIS, Temperature::DIMS, 298.15, 298.15),
        (VOLUME_AXIS, Volume::DIMS, 0.001, 10.0),
        (
            SATURATION_CONCENTRATION_AXIS,
            AmountConcentration::DIMS,
            0.0,
            10.0,
        ),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS, 0.001, 10.0),
        (PH_AXIS, Dims::NONE, -1.0, 15.0),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), low, high);
    }
    ConstitutiveModelCard {
        law: LawId(OXYGEN_TRANSFER_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [(
            KLA_PARAMETER.into(),
            LawParameter {
                value: kla,
                dims: Frequency::DIMS,
            },
        )]
        .into(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: "authored kLa; no aerator calibration".into(),
            license: "project test fixture".into(),
            artifact: Some(artifact),
        },
    }
}
fn transfer(kla: f64) -> OxygenTransfer {
    OxygenTransfer::from_card(solubility(100_000.0), transfer_card(kla)).unwrap()
}
fn bare(rate: f64, charge: CyanidationCharge) -> CyanidationCstr {
    CyanidationCstr::new(local(rate), charge, QueryPoint::new(), fixture::options()).unwrap()
}
fn feed(charge: &CyanidationCharge) -> CyanidationFeed {
    let c = |n: Amount| AmountConcentration::new(n.value() / charge.solution_volume.value());
    CyanidationFeed {
        accessible_gold: c(charge.accessible_gold),
        cyanide_ligand: c(charge.uncomplexed_cyanide),
        dissolved_oxygen: c(charge.dissolved_oxygen),
        gold_complex: c(charge.gold_complex),
        other_families: charge
            .other_families
            .iter()
            .map(|(id, n)| (id.clone(), c(*n)))
            .collect(),
        spectators: charge
            .spectators
            .iter()
            .map(|s| SpectatorFeed {
                species: s.species.clone(),
                charge: s.charge,
                concentration: c(s.amount),
            })
            .collect(),
        competing_metals: vec![],
    }
}
fn numerics() -> FlowNumerics {
    FlowNumerics {
        maximum_step: Time::new(0.25),
        absolute_amount_error: Amount::new(2e-10),
        relative_amount_error: Dimensionless::new(2e-8),
        maximum_work: 10_000_000,
        ..FlowNumerics::default()
    }
}
fn advance(
    tank: &mut AeratedCyanidationCstr,
    t: f64,
    q: f64,
    inlet: &CyanidationFeed,
    p: f64,
) -> AeratedFlowAdvance {
    tank.advance(
        Time::new(t),
        VolumetricFlowRate::new(q),
        inlet,
        Pressure::new(p),
        &QueryPoint::new(),
        numerics(),
        || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
    );
}
fn balanced(answer: &AeratedFlowAdvance) {
    for b in [
        &answer.flow.interval_conservation,
        &answer.flow.cumulative_conservation,
    ] {
        for value in [
            b.gold,
            b.dissolved_oxygen,
            b.sodium,
            b.components.accessible_gold,
            b.components.cyanide_ligand,
            b.components.dissolved_oxygen,
            b.components.gold_complex,
        ] {
            close(value.value(), 0.0, 2e-11);
        }
        let residual = answer
            .flow
            .state
            .competition
            .as_ref()
            .map_or(0.0, |c| c.ligand_balance_residual.value());
        close(b.cyanide.value(), residual, 2e-11);
        for (_, n) in b
            .components
            .other_families
            .iter()
            .chain(&b.components.competing_metals)
        {
            close(n.value(), 0.0, 2e-11);
        }
        for s in &b.components.spectators {
            close(s.amount.value(), 0.0, 2e-11);
        }
    }
    let f = &answer.flow;
    // Independently recompute the cumulative open-oxygen identity from the
    // public liquid streams, original charge, gas accounts and written equation.
    close(
        f.state.inventory.dissolved_oxygen.value()
            + f.cumulative_outlet.dissolved_oxygen.value()
            + answer.cumulative_oxygen_out.value()
            + f.cumulative_gold_dissolved.value() / 4.0
            - f.initial_charge().dissolved_oxygen.value()
            - f.cumulative_inlet.dissolved_oxygen.value()
            - answer.cumulative_oxygen_in.value(),
        0.0,
        2e-11,
    );
}

#[test]
fn g1_measured_henry_drives_continuous_flow_on_the_actual_volume_basis() {
    let source = measured_source(298.15);
    let context = henry_context(&source);
    let mut charge = fixture::charge(1.0);
    charge.solution_volume = Volume::new(0.4);
    charge.amount_at_unit_activity = Amount::new(0.2);
    charge.dissolved_oxygen = Amount::new(0.1);
    let mut inlet = feed(&charge);
    inlet.dissolved_oxygen = AmountConcentration::new(0.2);
    let mut tank = AeratedCyanidationCstr::from_tank(
        bare(0.0, charge),
        OxygenTransfer::from_henry_source(source.clone(), transfer_card(0.3)).unwrap(),
    )
    .unwrap();
    let answer = tank
        .advance(
            Time::new(0.7),
            VolumetricFlowRate::new(0.08),
            &inlet,
            Pressure::new(100_000.0),
            &context,
            numerics(),
            || false,
        )
        .unwrap();
    // Unsplit solution of dN/dt=Q*C_in+kLa*N_sat-(Q/V+kLa)*N.
    let equilibrium = (0.08 * 0.2 + 0.3 * 0.4) / (0.08 / 0.4 + 0.3);
    let expected = equilibrium + (0.1 - equilibrium) * (-0.35_f64).exp();
    close(answer.oxygen_saturation.value(), 0.4, 1e-15);
    close(
        answer.flow.state.inventory.dissolved_oxygen.value(),
        expected,
        2e-7,
    );
    assert!(answer.oxygen_in.value() > 0.0);
    assert!(answer.transfer().solubility().is_none());
    assert_eq!(
        answer.henry_answer.as_ref().unwrap().receipt.selected,
        source.answer().receipt.selected
    );
    balanced(&answer);
}
fn unchanged(tank: &AeratedCyanidationCstr, before: &AeratedCyanidationCstr) {
    assert_eq!(tank.tank.raw, before.tank.raw);
    assert_eq!(tank.tank.inlet, before.tank.inlet);
    assert_eq!(tank.tank.outlet, before.tank.outlet);
    assert_eq!(tank.tank.gold_dissolved, before.tank.gold_dissolved);
    assert_eq!(tank.tank.accumulated_error, before.tank.accumulated_error);
    assert_eq!(tank.tank.state.elapsed, before.tank.state.elapsed);
    assert_eq!(tank.tank.state.conditions, before.tank.state.conditions);
    assert_eq!(tank.tank.state.free_cyanide, before.tank.state.free_cyanide);
    assert_eq!(tank.oxygen_in, before.oxygen_in);
    assert_eq!(tank.oxygen_out, before.oxygen_out);
    assert_eq!(tank.tank.initial_charge(), before.tank.initial_charge());
    assert!(Arc::ptr_eq(&tank.transfer, &before.transfer));
}

#[test]
fn nonreacting_flow_and_gas_match_unsplit_solution_and_gross_flux_reversal() {
    // Independent linear ODE, including a crossing of n_sat during one call:
    // n'=q*n_feed+k*n_sat-(q+k)*n, n_inf=(q*n_feed+k*n_sat)/(q+k).
    for initial in [0.1, 0.72] {
        let mut charge = fixture::charge(1.0);
        charge.dissolved_oxygen = Amount::new(initial);
        let mut inlet = feed(&charge);
        inlet.dissolved_oxygen = AmountConcentration::new(1.2);
        let mut tank = AeratedCyanidationCstr::from_tank(bare(0.0, charge), transfer(0.3)).unwrap();
        let answer = advance(&mut tank, 6.0, 0.2, &inlet, 40_000.0);
        let (q, k, sat, time) = (0.2_f64, 0.3_f64, 0.4_f64, 6.0_f64);
        let lambda = q + k;
        let equilibrium = (q * 1.2 + k * sat) / lambda;
        let integral =
            |t: f64| equilibrium * t + (initial - equilibrium) * (-(-lambda * t).exp_m1()) / lambda;
        let net = |t: f64| k * (sat * t - integral(t));
        let expected = equilibrium + (initial - equilibrium) * (-lambda * time).exp();
        let (inward, outward) = if initial < sat {
            let crossing = -((sat - equilibrium) / (initial - equilibrium)).ln() / lambda;
            assert!(crossing > 0.0 && crossing < time);
            (net(crossing), net(crossing) - net(time))
        } else {
            (0.0, -net(time))
        };
        close(
            answer.flow.state.inventory.dissolved_oxygen.value(),
            expected,
            8e-7,
        );
        close(
            answer.flow.outlet.dissolved_oxygen.value(),
            q * integral(time),
            8e-7,
        );
        close(answer.oxygen_in.value(), inward, 8e-7);
        close(answer.oxygen_out.value(), outward, 8e-7);
        assert_eq!(answer.flow.interval_gold_dissolved.value(), 0.0);
        if initial < sat {
            assert!(answer.oxygen_in.value() > 0.0 && answer.oxygen_out.value() > 0.0);
        }
        balanced(&answer);
    }
}

#[test]
fn zero_transfer_matches_original_flow_and_retains_preattachment_history() {
    let original = fixture::charge(1.0);
    let inlet = feed(&original);
    let mut plain = bare(0.5, original.clone());
    plain
        .advance(
            Time::new(0.3),
            VolumetricFlowRate::new(0.2),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    let mut tank = AeratedCyanidationCstr::from_tank(plain.clone(), transfer(0.0)).unwrap();
    let expected = plain
        .advance(
            Time::new(1.0),
            VolumetricFlowRate::new(0.15),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    let answer = advance(&mut tank, 1.0, 0.15, &inlet, 40_000.0);
    for (a, b) in tank.tank.raw.iter().zip(&plain.raw) {
        close(*a, *b, 8e-8);
    }
    for (a, b) in tank.tank.outlet.iter().zip(&plain.outlet) {
        close(*a, *b, 8e-8);
    }
    close(
        answer.flow.cumulative_gold_dissolved.value(),
        expected.cumulative_gold_dissolved.value(),
        8e-8,
    );
    assert_eq!(answer.flow.initial_charge(), &original);
    assert_eq!(answer.flow.elapsed.value(), 1.3);
    assert_eq!(answer.oxygen_in.value(), 0.0);
    assert_eq!(answer.oxygen_out.value(), 0.0);
    assert_eq!(answer.transfer().source(), tank.transfer().source());
    assert!(
        !answer
            .flow
            .state
            .conditions
            .axes()
            .contains_key(ELAPSED_TIME_AXIS)
    );
    assert_eq!(answer.transfer_conditions.axes()[ELAPSED_TIME_AXIS], 1.3);
    balanced(&answer);
}

#[test]
fn stopped_flow_matches_batch_chemistry_and_keeps_all_liquid_exports_zero() {
    let mut charge = fixture::charge(1.0);
    charge.dissolved_oxygen = Amount::new(0.1);
    let inlet = feed(&charge);
    let mut tank =
        AeratedCyanidationCstr::from_tank(bare(0.5, charge.clone()), transfer(0.2)).unwrap();
    let batch_model = GoldCyanidation::from_card(
        fixture::system(false),
        "uncomplexed-cyanide",
        fixture::card(0.5),
    )
    .unwrap();
    let mut batch = AeratedGoldCyanidationBatch::new(
        batch_model,
        transfer(0.2),
        charge,
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let expected = batch
        .advance(
            Time::new(0.5),
            Pressure::new(40_000.0),
            &QueryPoint::new(),
            AerationNumerics {
                maximum_step: Time::new(0.02),
                relative_error: Dimensionless::new(1e-9),
                maximum_work: 10_000_000,
                ..AerationNumerics::default()
            },
            || false,
        )
        .unwrap();
    let answer = advance(&mut tank, 0.5, 0.0, &inlet, 40_000.0);
    close(
        answer.flow.state.inventory.accessible_gold.value(),
        expected.state.accessible_gold.value(),
        8e-8,
    );
    close(
        answer.flow.state.inventory.dissolved_oxygen.value(),
        expected.state.dissolved_oxygen.value(),
        8e-8,
    );
    close(answer.oxygen_in.value(), expected.oxygen_in.value(), 8e-8);
    assert!(answer.flow.interval_gold_dissolved.value() > 0.0);
    assert!(
        tank.tank
            .inlet
            .iter()
            .chain(&tank.tank.outlet)
            .all(|x| *x == 0.0)
    );
    balanced(&answer);
}

#[test]
fn gas_supply_enables_reaction_with_no_initial_or_feed_oxygen() {
    let mut charge = fixture::charge(1.0);
    charge.dissolved_oxygen = Amount::new(0.0);
    let inlet = feed(&charge);
    let mut plain = bare(1.0, charge.clone());
    let mut tank = AeratedCyanidationCstr::from_tank(bare(1.0, charge), transfer(0.3)).unwrap();
    let closed = plain
        .advance(
            Time::new(1.0),
            VolumetricFlowRate::new(0.1),
            &inlet,
            numerics(),
            || false,
        )
        .unwrap();
    let answer = advance(&mut tank, 1.0, 0.1, &inlet, 80_000.0);
    assert_eq!(closed.interval_gold_dissolved.value(), 0.0);
    assert!(answer.flow.interval_gold_dissolved.value() > 0.01);
    assert!(answer.flow.outlet.gold_complex.value() > 0.0);
    assert!(answer.flow.outlet.dissolved_oxygen.value() > 0.0);
    assert!(answer.oxygen_in.value() > answer.flow.interval_gold_dissolved.value() / 4.0);
    assert_eq!(answer.oxygen_out.value(), 0.0);
    balanced(&answer);
}

#[test]
fn copper_competition_survives_pressure_flow_and_feed_changes() {
    let charge = fixture::competition::copper_charge(0.2);
    let batch = fixture::competition::competing_batch(charge.clone(), 0.2);
    let bound = batch.layout.competition.as_ref().unwrap();
    let metals = competition::CompetingMetals {
        families: bound.model.families().to_vec(),
        totals: bound.totals.clone(),
        numerics: bound.numerics,
    };
    let metal_id = metals.totals[0].0.clone();
    let mut inlet = feed(&charge);
    inlet.competing_metals = vec![(metal_id, AmountConcentration::new(0.1))];
    let liquid = CyanidationCstr::new_with_competing_metals(
        local(0.5),
        charge.clone(),
        QueryPoint::new(),
        metals,
    )
    .unwrap();
    let sources = liquid.competing_chemistry().unwrap().families()[0].clone();
    let mut tank = AeratedCyanidationCstr::from_tank(liquid, transfer(0.3)).unwrap();
    let first = advance(&mut tank, 0.125, 0.5, &inlet, 60_000.0);
    assert!(first.oxygen_in.value() > 0.0);
    balanced(&first);
    let stopped = advance(&mut tank, 0.25, 0.0, &inlet, 0.0);
    assert!(stopped.oxygen_out.value() > 0.0);
    assert_eq!(stopped.flow.outlet.dissolved_oxygen.value(), 0.0);
    assert_eq!(stopped.flow.outlet.gold_complex.value(), 0.0);
    balanced(&stopped);
    inlet.cyanide_ligand = AmountConcentration::new(2.2);
    inlet
        .spectators
        .iter_mut()
        .find(|s| s.species.as_str() == dosing::SODIUM)
        .unwrap()
        .concentration = AmountConcentration::new(1.2);
    let last = advance(&mut tank, 0.125, 0.2, &inlet, 60_000.0);
    let copper = last.flow.state.competition.as_ref().unwrap();
    close(
        copper.metals[0].total.value(),
        0.1 + 0.1 * (-0.0875_f64).exp(),
        3e-12,
    );
    close(
        last.flow.state.uncomplexed_cyanide.value() + copper.bound_ligand.value(),
        last.flow.state.inventory.cyanide_ligand.value(),
        3e-9,
    );
    assert!(copper.bound_ligand.value() > 0.0);
    assert!(last.cumulative_oxygen_in.value() > 0.0 && last.cumulative_oxygen_out.value() > 0.0);
    assert_eq!(tank.tank.initial_charge(), &charge);
    assert_eq!(
        tank.tank.competing_chemistry().unwrap().families()[0].id(),
        sources.id()
    );
    balanced(&last);
}

#[test]
fn unsupported_boundaries_work_exhaustion_and_cancellation_publish_nothing() {
    let charge = fixture::charge(1.0);
    let inlet = feed(&charge);
    let mut tank =
        AeratedCyanidationCstr::from_tank(bare(0.5, charge.clone()), transfer(0.3)).unwrap();
    advance(&mut tank, 0.2, 0.1, &inlet, 60_000.0);
    let before = tank.clone();
    for (duration, q, pressure, extra, controls) in [
        (0.1, 0.1, 300_000.0, QueryPoint::new(), numerics()),
        (101.0, 0.1, 60_000.0, QueryPoint::new(), numerics()),
        (-0.1, 0.1, 60_000.0, QueryPoint::new(), numerics()),
        (0.1, -0.1, 60_000.0, QueryPoint::new(), numerics()),
        (
            0.1,
            0.1,
            60_000.0,
            QueryPoint::new()
                .with_quantity(
                    OXYGEN_AXIS,
                    QuantitySpec::dimensional(AmountConcentration::DIMS),
                    0.1,
                )
                .unwrap(),
            numerics(),
        ),
        (
            0.1,
            0.1,
            60_000.0,
            QueryPoint::new(),
            FlowNumerics {
                maximum_work: 1,
                ..numerics()
            },
        ),
        (
            0.1,
            0.1,
            60_000.0,
            QueryPoint::new(),
            FlowNumerics {
                maximum_attempts: 1,
                maximum_step: Time::new(0.01),
                ..numerics()
            },
        ),
    ] {
        assert!(
            tank.advance(
                Time::new(duration),
                VolumetricFlowRate::new(q),
                &inlet,
                Pressure::new(pressure),
                &extra,
                controls,
                || false
            )
            .is_err()
        );
        unchanged(&tank, &before);
    }
    let mut missing = inlet.clone();
    missing.spectators.clear();
    assert!(
        tank.advance(
            Time::new(0.1),
            VolumetricFlowRate::new(0.1),
            &missing,
            Pressure::new(60_000.0),
            &QueryPoint::new(),
            numerics(),
            || false
        )
        .is_err()
    );
    unchanged(&tank, &before);
    let mut calls = 0;
    assert!(matches!(
        tank.advance(
            Time::new(1.0),
            VolumetricFlowRate::new(0.1),
            &inlet,
            Pressure::new(60_000.0),
            &QueryPoint::new(),
            numerics(),
            || {
                calls += 1;
                calls >= 20
            }
        ),
        Err(CyanidationError::Cancelled)
    ));
    unchanged(&tank, &before);
    assert!(
        AeratedCyanidationCstr::from_tank(
            bare(0.5, charge),
            OxygenTransfer::from_card(solubility(101_000.0), transfer_card(0.3)).unwrap()
        )
        .is_err()
    );
}
