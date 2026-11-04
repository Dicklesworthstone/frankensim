use super::*;
use crate::phase::{
    EnthalpyPhaseKnot, HeatCapacityKnot, ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_thermal::crystallization::CrystallizationSolve;
use fs_thermal::enthalpy::crystallization::{
    CrystallizationControls, CrystallizationMaterial, CrystallizationSimulation,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, TimedEvent,
};

fn branch(phase: SolidLiquidPhase, cp: f64, h0: f64, rho: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([41; 32]),
        phase,
        h0,
        &[
            HeatCapacityKnot {
                temperature_k: 280.0,
                specific_heat_capacity_j_kg_k: cp,
                bulk_density_kg_m3: rho,
            },
            HeatCapacityKnot {
                temperature_k: 380.0,
                specific_heat_capacity_j_kg_k: cp,
                bulk_density_kg_m3: rho,
            },
        ],
        1e-6,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_distinct_original_cp_and_phase_volumes_survive_crystallization() {
    let parent = branch(SolidLiquidPhase::Liquid, 2000.0, 100_000.0, 800.0);
    let crystal = branch(SolidLiquidPhase::Solid, 1000.0, 0.0, 1000.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    assert_eq!(pair.parent().identity(), parent.identity());
    assert_eq!(pair.crystal().identity(), crystal.identity());
    assert_eq!(pair.maximum_temperature_slope(), 0.001);
    assert_eq!(pair.latent_enthalpy_j_kg(300.0).unwrap(), 120_000.0);
    assert_eq!(pair.latent_enthalpy_j_kg(350.0).unwrap(), 170_000.0);

    let rate = |_: f64| Ok(0.2);
    let law = pair
        .law(&rate, ContentHash([42; 32]), 2.0, 0.6, 400.0)
        .unwrap();
    let initial = law.state_at_temperature(300.0, 0.0).unwrap();
    let result = law
        .advance_specific_enthalpy(
            initial,
            initial.specific_enthalpy_j_kg(),
            1.0,
            CrystallizationSolve::default(),
            || false,
        )
        .unwrap();
    let state = result.state;
    let x = 0.6 * (1.0 - (-0.04_f64).exp());
    let expected_t = 280.0 + (140_000.0 - (1.0 - x) * 100_000.0) / (2000.0 - 1000.0 * x);
    close(state.crystalline_mass_fraction(), x, 2e-15);
    close(state.temperature_k(), expected_t, 2e-10);
    assert_eq!(
        state.specific_enthalpy_j_kg(),
        initial.specific_enthalpy_j_kg()
    );
    let volumes = pair.volume_state(state.temperature_k(), x).unwrap();
    let v = (1.0 - x) / 800.0 + x / 1000.0;
    close(volumes.specific_volume_m3_kg, v, 1e-18);
    close(volumes.crystalline_volume_fraction, (x / 1000.0) / v, 2e-15);
    assert!(volumes.crystalline_volume_fraction < x);
    assert!(volumes.bulk_density_kg_m3 > 800.0);
}

#[test]
fn g0_never_infers_metastability_from_equilibrium_plateau_or_uncovered_sources() {
    let parent = branch(SolidLiquidPhase::Liquid, 2000.0, 100_000.0, 800.0);
    let crystal = branch(SolidLiquidPhase::Solid, 1000.0, 0.0, 1000.0);
    assert!(CrystallizationCaloricPair::new(&parent, &crystal, [279.0, 380.0]).is_err());
    assert!(CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 381.0]).is_err());
    assert!(CrystallizationCaloricPair::new(&parent, &parent, [280.0, 380.0]).is_err());
    let wrong_datum = branch(SolidLiquidPhase::Solid, 2000.0, 120_000.0, 1000.0);
    assert!(CrystallizationCaloricPair::new(&parent, &wrong_datum, [280.0, 380.0]).is_err());
    let plateau = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([41; 32]),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 280.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 100.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 1000.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 800.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 1100.0,
                temperature_k: 380.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 800.0,
            },
        ],
    )
    .unwrap();
    assert!(CrystallizationCaloricBranch::new(&plateau).is_err());
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    for x in [-0.1, 1.1, f64::NAN] {
        assert!(pair.volume_state(300.0, x).is_err());
    }
    assert!(pair.volume_state(500.0, 0.5).is_err());
    let rate = |_: f64| Ok(0.2);
    assert!(
        pair.law(&rate, ContentHash([42; 32]), 2.0, 0.6, 370.0)
            .is_err()
    );
}

#[test]
fn g3_energy_datum_and_kinetic_source_identity_are_separate_from_physics() {
    let rate = |_: f64| Ok(0.2);
    let run = |datum| {
        let parent = branch(SolidLiquidPhase::Solid, 2000.0, datum + 100_000.0, 800.0);
        let crystal = branch(SolidLiquidPhase::Solid, 1000.0, datum, 1000.0);
        let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
        let law = pair
            .law(&rate, ContentHash([42; 32]), 2.0, 0.6, 400.0)
            .unwrap();
        let other = pair
            .law(&rate, ContentHash([43; 32]), 2.0, 0.6, 400.0)
            .unwrap();
        assert_ne!(law.identity(), other.identity());
        let initial = law.state_at_temperature(300.0, 0.0).unwrap();
        assert!(
            other
                .advance_specific_enthalpy(
                    initial,
                    initial.specific_enthalpy_j_kg(),
                    1.0,
                    CrystallizationSolve::default(),
                    || false
                )
                .is_err()
        );
        law.advance_specific_enthalpy(
            initial,
            initial.specific_enthalpy_j_kg() - 500.0,
            1.0,
            CrystallizationSolve {
                specific_enthalpy_tolerance_j_kg: 1e-6,
                ..CrystallizationSolve::default()
            },
            || false,
        )
        .unwrap()
        .state
    };
    let a = run(0.0);
    let b = run(1e7);
    close(a.temperature_k(), b.temperature_k(), 3e-9);
    assert_eq!(a.crystalline_mass_fraction(), b.crystalline_mass_fraction());
}

#[test]
fn g1_original_materials_drive_cooling_and_recalescence_on_shared_timeline() {
    let parent = branch(SolidLiquidPhase::Liquid, 2000.0, 100_000.0, 800.0);
    let crystal = branch(SolidLiquidPhase::Solid, 1000.0, 0.0, 1000.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let rate = |_: f64| Ok(0.2);
    let law = pair
        .law(&rate, ContentHash([42; 32]), 2.0, 0.6, 400.0)
        .unwrap();
    let original_h = pair.parent().enthalpy_at_temperature(300.0).unwrap();
    let cell = DepositedCell {
        id: CellId(19),
        segment_id: 3,
        body_id: 7,
        material_id: MaterialId(11),
        mass_kg: 0.01,
        specific_h_j_kg: original_h,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [0.01, 0.0, 0.0],
            volume_m3: 0.01 / 800.0,
            surface_area_m2: 0.01,
        },
        boundary: Boundary {
            ambient_k: 280.0,
            htc_w_m2_k: 5.0,
            reservoir: None,
            source_w: 0.0,
        },
    };
    let thermal = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(11), pair.parent())],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let simulation = CrystallizationSimulation::new(
        thermal,
        &[CrystallizationMaterial {
            material_id: MaterialId(11),
            amorphous: pair.parent(),
            law: &law,
            maximum_temperature_slope: pair.maximum_temperature_slope(),
        }],
        || Ok(()),
    )
    .unwrap();
    let (cooled, audit) = simulation
        .advance_to(2.0, CrystallizationControls::default(), || Ok(()))
        .unwrap();
    let state = cooled.states()[&CellId(19)];
    assert!(audit.thermal.external_heat_j < 0.0);
    assert!(state.specific_enthalpy_j_kg() < original_h);
    assert!(
        state.temperature_k() > 300.0,
        "latent heat must feed back into actual T"
    );
    close(
        audit.crystallized_mass_kg,
        0.01 * state.crystalline_mass_fraction(),
        1e-15,
    );
    close(
        audit.thermal.stored_energy_change_j,
        audit.thermal.external_heat_j,
        1e-8,
    );
    assert_eq!(cooled.mass_kg(), simulation.mass_kg());
    assert_eq!(cooled.cells()[&CellId(19)].body_id, 7);
    let original_volume = pair.volume_state(300.0, 0.0).unwrap();
    let current_volume = pair
        .volume_state(state.temperature_k(), state.crystalline_mass_fraction())
        .unwrap();
    assert!(current_volume.specific_volume_m3_kg < original_volume.specific_volume_m3_kg);
    // Reference transport geometry is intentionally distinct from free volume.
    assert_eq!(
        cooled.cells()[&CellId(19)].geometry,
        simulation.cells()[&CellId(19)].geometry
    );
}

#[test]
fn g0_sourced_heat_capacity_and_density_receipts_remain_original() {
    use crate::phase::HeatCapacitySourceQuery;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    };
    use fs_qty::{Dims, QuantitySpec, Temperature};
    let cp_dims = Dims([2, 0, -2, -1, 0, 0]);
    let rho_dims = Dims([-3, 1, 0, 0, 0, 0]);
    let rate_dims = Dims([0, 0, -1, 0, 0, 0]);
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("parent-cp", 2000.0, cp_dims),
        ("crystal-cp", 1000.0, cp_dims),
        ("parent-rho", 800.0, rho_dims),
        ("crystal-rho", 1000.0, rho_dims),
        ("K", 0.1, rate_dims),
        ("invalid-K", -0.1, rate_dims),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: if name == "K" {
                        vec![(280.0, 0.0), (320.0, 0.3), (350.0, 0.15), (380.0, 0.0)]
                    } else {
                        vec![(280.0, value), (380.0, value)]
                    },
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    280.0,
                    380.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: fs_matdb::UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured crystallization fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "declared parent and crystal".into(),
            process: "caloric fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let low = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 280.0)
        .unwrap();
    let high = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 380.0)
        .unwrap();
    let parent = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("parent-cp", cp_dims),
            density: &PropertyKey::new("parent-rho", rho_dims),
            temperature_axis: "T",
            lower: &low,
            upper: &high,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        100_000.0,
        1e-6,
    )
    .unwrap();
    let crystal = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("crystal-cp", cp_dims),
            density: &PropertyKey::new("crystal-rho", rho_dims),
            temperature_axis: "T",
            lower: &low,
            upper: &high,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    let pair =
        CrystallizationCaloricPair::new(parent.curve(), crystal.curve(), [280.0, 380.0]).unwrap();
    assert!(!pair.parent().source_receipts().is_empty());
    assert_eq!(pair.parent().source_receipts(), parent.receipts());
    assert_eq!(pair.crystal().source_receipts(), crystal.receipts());
    for t in [280.0, 290.0, 333.0, 380.0] {
        let h = pair.parent().enthalpy_at_temperature(t).unwrap();
        close(
            parent.curve().temperature_at_specific_enthalpy(h).unwrap(),
            t,
            1e-13,
        );
        assert_eq!(pair.parent().density_at_temperature(t).unwrap(), 800.0);
    }
    use super::kinetics::PreparedCrystallizationRate;
    use fs_matdb::ClaimSelection;
    let selected = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let rate = PreparedCrystallizationRate::new(
        &card,
        &PropertyKey::new("K", rate_dims),
        "T",
        &low,
        &high,
        selected,
    )
    .unwrap();
    assert_eq!(rate.source().knots.len(), 4);
    close(rate.rate_s_inv(335.0).unwrap(), 0.225, 1e-15);
    assert_eq!(rate.rate_s_inv(380.0).unwrap(), 0.0);
    assert!(rate.rate_s_inv(400.0).is_err());
    assert!(
        PreparedCrystallizationRate::new(
            &card,
            &PropertyKey::new("parent-cp", cp_dims),
            "T",
            &low,
            &high,
            selected
        )
        .is_err()
    );
    assert!(
        PreparedCrystallizationRate::new(
            &card,
            &PropertyKey::new("invalid-K", rate_dims),
            "T",
            &low,
            &high,
            selected
        )
        .is_err()
    );
    let law = pair.law_with_sourced_rate(&rate, 2.0, 0.6, 380.0).unwrap();
    assert_eq!(law.rate_s_inv(380.0).unwrap(), 0.0);
    let initial = law.state_at_temperature(335.0, 0.0).unwrap();
    let final_state = law.advance_isothermal(initial, 1.0).unwrap().state;
    close(
        final_state.crystalline_mass_fraction(),
        0.6 * (1.0 - (-0.225_f64.powi(2)).exp()),
        1e-15,
    );
}
