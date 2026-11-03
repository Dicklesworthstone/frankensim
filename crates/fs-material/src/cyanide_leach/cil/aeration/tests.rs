//! Coupled gas/CIL limits and open oxygen accounts.

use super::*;
use crate::cyanide_leach::{
    aeration::{KLA_PARAMETER, OXYGEN_TRANSFER_LAW},
    cil::tests::{cil, close, numerics},
    tests as fixture,
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::{AmountConcentration, Dims, Frequency, QuantitySpec, Temperature, Volume};
use fs_thermochem::stoichiometry::thermodynamics::{
    FormationGibbsSource, gas_liquid::IdealGasLiquidPartition,
};

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
    let artifact = hash_domain("aerated-cil-known-answer-v1", b"authored kLa");
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

#[test]
fn opening_oxygen_rechecks_the_expanded_counterion_demand() {
    let mut initial = fixture::charge(0.1);
    initial.dissolved_oxygen = Amount::new(0.0);
    let closed = cil(initial, 1.0);
    // Closed-bath progress is oxygen-limited to zero, so 0.1 mol Na+ is
    // sufficient for that owner. Aeration raises the Au/CN ceiling to 1 mol
    // and must refuse rather than invent the missing conserved counterion.
    assert!(matches!(
        AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.3)),
        Err(CilError::Invalid(
            "insufficient counterion reserve for aerated CIL gold loading"
        ))
    ));
}

#[test]
fn zero_transfer_matches_original_simultaneous_cil() {
    let initial = fixture::charge(1.0);
    let closed = cil(initial, 1.0);
    let mut expected = closed.clone();
    let expected = expected
        .advance(Time::new(0.4), numerics(1.0), || false)
        .unwrap();
    let mut actual = AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.0)).unwrap();
    let answer = actual
        .advance(
            Time::new(0.4),
            Pressure::new(40_000.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .unwrap();
    close(
        answer.cil.state.solution.gold_dissolved.value(),
        expected.state.solution.gold_dissolved.value(),
        3e-8,
    );
    close(
        answer.cil.state.solution.gold_complex.value(),
        expected.state.solution.gold_complex.value(),
        3e-8,
    );
    for (a, b) in answer
        .cil
        .state
        .populations
        .iter()
        .zip(&expected.state.populations)
    {
        close(a.adsorbed.value(), b.adsorbed.value(), 3e-8);
    }
    assert_eq!(answer.interval_oxygen_in.value(), 0.0);
    assert_eq!(answer.interval_oxygen_out.value(), 0.0);
    close(answer.cil.conservation.dissolved_oxygen.value(), 0.0, 2e-11);
}

#[test]
fn gas_supply_enables_gold_dissolution_and_carbon_loading_without_liquid_oxygen() {
    let mut initial = fixture::charge(1.0);
    initial.dissolved_oxygen = Amount::new(0.0);
    let closed = cil(initial, 1.0);
    let mut baseline = closed.clone();
    let no_gas = baseline
        .advance(Time::new(0.5), numerics(1.0), || false)
        .unwrap();
    assert_eq!(no_gas.interval_gold_dissolved.value(), 0.0);

    let mut aerated = AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.3)).unwrap();
    let answer = aerated
        .advance(
            Time::new(0.5),
            Pressure::new(80_000.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .unwrap();
    assert!(answer.cil.interval_gold_dissolved.value() > 0.0);
    assert!(answer.cil.interval_product_capture.value() > 0.0);
    assert!(answer.interval_oxygen_in.value() > 0.0);
    assert_eq!(answer.interval_oxygen_out.value(), 0.0);
    close(answer.cil.conservation.gold.value(), 0.0, 2e-11);
    close(answer.cil.conservation.counterion.value(), 0.0, 2e-11);
    close(answer.cil.conservation.dissolved_oxygen.value(), 0.0, 2e-11);
}

#[test]
fn zero_pressure_strips_oxygen_while_retaining_carbon_and_atomic_history() {
    let mut initial = fixture::charge(1.0);
    initial.accessible_gold = Amount::new(0.0);
    let closed = cil(initial, 1.0);
    let mut aerated = AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.4)).unwrap();
    let first = aerated
        .advance(
            Time::new(0.2),
            Pressure::new(0.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .unwrap();
    assert_eq!(first.cil.interval_gold_dissolved.value(), 0.0);
    assert_eq!(first.interval_oxygen_in.value(), 0.0);
    assert!(first.interval_oxygen_out.value() > 0.0);
    let before = format!("{:?}", aerated.state());
    let cumulative_out = aerated.oxygen_out();
    assert!(
        aerated
            .advance(
                Time::new(101.0),
                Pressure::new(0.0),
                &QueryPoint::new(),
                numerics(1.0),
                || false,
            )
            .is_err()
    );
    assert_eq!(format!("{:?}", aerated.state()), before);
    assert_eq!(aerated.oxygen_out(), cumulative_out);
}
