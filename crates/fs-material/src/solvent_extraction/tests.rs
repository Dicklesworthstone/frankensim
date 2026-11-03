use super::*;
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::Dims;
use std::collections::BTreeMap;

fn copper() -> ElementId {
    ElementId::new("Cu").unwrap()
}
fn source(knots: &[(f64, f64)]) -> ConstitutiveModelCard {
    let mut parameters = BTreeMap::new();
    for (i, &(x, y)) in knots.iter().enumerate() {
        for (prefix, value) in [(AQUEOUS_KNOT_PREFIX, x), (ORGANIC_KNOT_PREFIX, y)] {
            parameters.insert(
                format!("{prefix}{i}"),
                LawParameter {
                    value,
                    dims: AmountConcentration::DIMS,
                },
            );
        }
    }
    ConstitutiveModelCard {
        law: LawId(ISOTHERM_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                AQUEOUS_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                knots[0].0,
                knots.last().unwrap().0,
            )
            .with_quantity(
                ORGANIC_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                knots[0].1,
                knots.last().unwrap().1,
            )
            .with_quantity(
                VOLUME_RATIO_AXIS,
                QuantitySpec::dimensional(Dimensionless::DIMS),
                0.01,
                100.0,
            )
            .with_quantity(
                TEMPERATURE_AXIS,
                QuantitySpec::dimensional(Temperature::DIMS),
                298.15,
                298.15,
            )
            .with_quantity(
                "authored-chemistry",
                QuantitySpec::dimensional(Dims::NONE),
                7.0,
                7.0,
            ),
        sources: vec![ContentHash([83; 32])],
        provenance: Provenance {
            source: "authored equilibrium table; not ore or extractant calibration".into(),
            license: "project test fixture".into(),
            artifact: Some(ContentHash([83; 32])),
        },
    }
}
fn model(knots: &[(f64, f64)]) -> SolventExtractionIsotherm {
    SolventExtractionIsotherm::from_card(copper(), source(knots)).unwrap()
}
fn phase(volume: f64, amount: f64) -> PhaseInventory {
    PhaseInventory::new(copper(), Volume::new(volume), Amount::new(amount)).unwrap()
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "authored-chemistry",
            QuantitySpec::dimensional(Dims::NONE),
            7.0,
        )
        .unwrap()
}
fn run(
    model: &SolventExtractionIsotherm,
    aqueous: &PhaseInventory,
    organic: &PhaseInventory,
) -> Result<EquilibriumStage> {
    model.equilibrate(aqueous, organic, Temperature::new(298.15), &conditions())
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-14 * expected.abs().max(1e-300),
        "{actual} != {expected}"
    );
}

#[test]
fn loaded_unequal_phases_match_linear_balance_and_keep_original_source() {
    let original = source(&[(0.0, 0.0), (1.0, 2.0), (4.0, 8.0)]);
    let model = SolventExtractionIsotherm::from_card(copper(), original.clone()).unwrap();
    let aqueous = phase(2.0, 10.0); // Initial C=5 exceeds equilibrium-table support; valid feed.
    let organic = phase(3.0, 4.0);
    let a_before = aqueous.clone();
    let o_before = organic.clone();
    let result = run(&model, &aqueous, &organic).unwrap();
    // Independent balance: x = (10+4)/(2+3*2) = 1.75, y=2x=3.5.
    close(result.aqueous.concentration().value(), 1.75);
    close(result.organic.concentration().value(), 3.5);
    close(result.aqueous.amount().value(), 3.5);
    close(result.organic.amount().value(), 10.5);
    close(result.transfer_to_organic.value(), 6.5);
    assert_eq!(result.conservation_residual.value(), 0.0);
    assert_eq!(result.model().source(), &original);
    assert_eq!(result.model().solute(), &copper());
    assert_eq!(result.conditions.axes()[VOLUME_RATIO_AXIS], 1.5);
    assert_eq!(result.segment_index, 1);
    assert_eq!((aqueous, organic), (a_before, o_before));
}

#[test]
fn nonlinear_segment_flat_loading_exact_knot_and_empty_contact_have_known_answers() {
    let nonlinear = model(&[(0.0, 0.0), (1.0, 4.0), (4.0, 7.0)]);
    let r = run(&nonlinear, &phase(2.0, 12.0), &phase(1.0, 0.0)).unwrap();
    // Second segment y=x+3: 2x+(x+3)=12 -> x=3, y=6.
    assert_eq!(r.segment_index, 1);
    close(r.aqueous.amount().value(), 6.0);
    close(r.organic.amount().value(), 6.0);
    let knot = run(&nonlinear, &phase(2.0, 6.0), &phase(1.0, 0.0)).unwrap();
    assert_eq!(knot.segment_index, 0); // Deterministic left segment at an internal knot.
    close(knot.aqueous.concentration().value(), 1.0);
    close(knot.organic.concentration().value(), 4.0);
    let flat = model(&[(0.0, 0.0), (1.0, 4.0), (3.0, 4.0)]);
    let r = run(&flat, &phase(2.0, 8.0), &phase(1.0, 0.0)).unwrap();
    close(r.aqueous.concentration().value(), 2.0);
    close(r.organic.concentration().value(), 4.0);
    let empty = run(&flat, &phase(2.0, 0.0), &phase(1.0, 0.0)).unwrap();
    assert_eq!(empty.aqueous.amount().value(), 0.0);
    assert_eq!(empty.organic.amount().value(), 0.0);
}

#[test]
fn stripping_reverses_transfer_and_common_inventory_volume_scaling_preserves_equilibrium() {
    let extract = model(&[(0.0, 0.0), (4.0, 16.0)]);
    let loaded = run(&extract, &phase(2.0, 12.0), &phase(1.0, 0.0)).unwrap();
    close(loaded.organic.amount().value(), 8.0); // x=12/(2+4)=2; y=8.
    let strip = model(&[(0.0, 0.0), (16.0, 4.0)]); // Same y-versus-x axis convention.
    let stripped = run(&strip, &phase(0.25, 0.0), &loaded.organic).unwrap();
    close(stripped.aqueous.concentration().value(), 16.0);
    close(stripped.aqueous.amount().value(), 4.0);
    close(stripped.organic.amount().value(), 4.0);
    close(stripped.transfer_to_organic.value(), -4.0);
    assert_eq!(stripped.conservation_residual.value(), 0.0);
    let repeat = run(&strip, &stripped.aqueous, &stripped.organic).unwrap();
    assert_eq!(repeat.transfer_to_organic.value(), 0.0);
    for scale in [0.001, 1_000.0] {
        let r = run(
            &extract,
            &phase(2.0 * scale, 12.0 * scale),
            &phase(scale, 0.0),
        )
        .unwrap();
        close(
            r.aqueous.concentration().value(),
            loaded.aqueous.concentration().value(),
        );
        close(r.organic.amount().value(), 8.0 * scale);
    }
}

#[test]
fn no_extrapolation_inferred_origin_or_spoofed_source_conditions() {
    let model = model(&[(1.0, 2.0), (4.0, 8.0)]);
    for amount in [0.0, 2.0, 13.0] {
        // Supported total is [3,12] at unit volumes.
        assert_eq!(
            run(&model, &phase(1.0, amount), &phase(1.0, 0.0)),
            Err(SolventExtractionError::OutsideIsotherm)
        );
    }
    let a = phase(1.0, 6.0);
    let o = phase(1.0, 0.0);
    assert_eq!(
        model.equilibrate(&a, &o, Temperature::new(300.0), &conditions()),
        Err(SolventExtractionError::OutsideValidity)
    );
    assert_eq!(
        model.equilibrate(&a, &o, Temperature::new(298.15), &QueryPoint::new()),
        Err(SolventExtractionError::OutsideValidity)
    );
    let untyped = QueryPoint::new().with("authored-chemistry", 7.0).unwrap();
    assert_eq!(
        model.equilibrate(&a, &o, Temperature::new(298.15), &untyped),
        Err(SolventExtractionError::OutsideValidity)
    );
    let spoofed = conditions().with(AQUEOUS_CONCENTRATION_AXIS, 1.0).unwrap();
    assert!(matches!(
        model.equilibrate(&a, &o, Temperature::new(298.15), &spoofed),
        Err(SolventExtractionError::Invalid(_))
    ));
    assert_eq!(
        run(&model, &phase(1.0, 2.0), &phase(0.001, 0.0)),
        Err(SolventExtractionError::OutsideValidity)
    );
}

#[test]
fn malformed_tables_dimensions_and_missing_evidence_are_refused_without_repair() {
    for mutation in 0..7 {
        let mut s = source(&[(0.0, 0.0), (1.0, 2.0), (4.0, 8.0)]);
        match mutation {
            0 => s.parameters.get_mut("aqueous-knot-1").unwrap().value = 0.0,
            1 => s.parameters.get_mut("organic-knot-2").unwrap().value = 1.0,
            2 => s.parameters.get_mut("organic-knot-1").unwrap().dims = Volume::DIMS,
            3 => s.sources.clear(),
            4 => s.law_version = 2,
            5 => {
                s.parameters.remove("organic-knot-1");
            }
            _ => s.validity = ValidityDomain::unconstrained(),
        }
        assert!(SolventExtractionIsotherm::from_card(copper(), s).is_err());
    }
    let m = model(&[(0.0, 0.0), (4.0, 8.0)]);
    let iron = PhaseInventory::new(
        ElementId::new("Fe").unwrap(),
        Volume::new(1.0),
        Amount::new(2.0),
    )
    .unwrap();
    assert!(matches!(
        run(&m, &phase(1.0, 1.0), &iron),
        Err(SolventExtractionError::Invalid(_))
    ));
}

#[test]
fn finite_inventory_checks_and_small_positive_phase_loading_are_explicit() {
    for (volume, amount) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (1.0, -1.0),
        (1.0, f64::NAN),
        (1e-300, 1e300),
    ] {
        assert!(PhaseInventory::new(copper(), Volume::new(volume), Amount::new(amount)).is_err());
    }
    let m = model(&[(0.0, 0.0), (1.0, 1e-20)]);
    let r = run(&m, &phase(1.0, 1.0), &phase(1.0, 0.0)).unwrap();
    // A small extracted population is not computed as 1 minus a rounded recovery.
    close(r.organic.amount().value(), 1e-20);
    assert!(r.conservation_residual.value().abs() <= 32.0 * f64::EPSILON);
    assert!(matches!(
        run(&m, &phase(1.0, 1e308), &phase(1.0, 1e308)),
        Err(SolventExtractionError::Arithmetic(_))
    ));
}
