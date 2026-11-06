//! Original ternary calorics and conductivity on the native thermal timeline.

use super::*;
use fs_material::liquid::solution::multicomponent::caloric::storage::{
    MulticomponentLiquidThermalStorage, conduction::MulticomponentLiquidConductionLaw,
};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionLaw,
    ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

#[path = "storage/stream.rs"]
mod stream;

fn thermal_controls() -> MulticomponentCaloricSolveV1 {
    MulticomponentCaloricSolveV1 {
        temperature_bracket: [Temperature::new(300.0), Temperature::new(600.0)],
        ..controls()
    }
}
fn variable_sources(datum: f64) -> Vec<(SpeciesId, IsobaricHeatCapacityLaw)> {
    (0..3)
        .map(|i| {
            (
                species(i),
                cp_source(
                    i,
                    [[80.0, 100.0, 150.0], [120.0; 3], [200.0; 3]][i],
                    100_000.0,
                    datum,
                    0,
                ),
            )
        })
        .collect()
}
// Original k is a composition/temperature surface, independent of the pure
// source Cp knots. Along each fixed composition k=base(T)+2*x_A.
fn thermal_card(pressure_supported: bool) -> MaterialCard {
    let original = card(CompositionBasis::MoleFraction, 1.0);
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut source = source.clone();
        if source.key.name() == THERMAL_CONDUCTIVITY_PROPERTY {
            source.value = PropertyValue::Surface2D {
                axes: ["T".into(), "x_a".into()],
                axis_dims: [Temperature::DIMS, Dimensionless::DIMS],
                samples: vec![
                    ([300.0, 0.0], 1.0),
                    ([450.0, 0.0], 4.0),
                    ([600.0, 0.0], 2.0),
                    ([300.0, 1.0], 3.0),
                    ([450.0, 1.0], 6.0),
                    ([600.0, 1.0], 4.0),
                ],
                triangles: vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
                dims: source.key.quantity().dims(),
            };
            if !pressure_supported {
                source.validity = ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        300.0,
                        600.0,
                    )
                    .with_quantity("x_a", composition(CompositionBasis::MoleFraction), 0.0, 1.0)
                    .with_quantity("x_b", composition(CompositionBasis::MoleFraction), 0.0, 1.0);
            }
        }
        claims.insert_claim(source).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored ternary thermal surface".into(),
            phase: "liquid".into(),
            process: "synthetic conductivity, not measured data".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (tolerance {tolerance:.3e})"
    );
}
fn h(storage: &MulticomponentLiquidThermalStorage, t: f64) -> f64 {
    storage
        .specific_enthalpy_j_kg(Temperature::new(t), || false)
        .unwrap()
}
fn k(t: f64, x_a: f64) -> (f64, f64) {
    if t < 450.0 {
        (1.0 + 0.02 * (t - 300.0) + 2.0 * x_a, 0.02)
    } else {
        (4.0 - (t - 450.0) / 75.0 + 2.0 * x_a, -1.0 / 75.0)
    }
}

#[test]
fn homogeneous_storage_preserves_original_cp_entropy_mass_and_reference_changes() {
    let card = thermal_card(true);
    let liquid = liquid(&card);
    let law =
        MulticomponentLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0), "p").unwrap();
    let shifted =
        MulticomponentLiquidCaloricLaw::try_new(&liquid, variable_sources(350.0), "p").unwrap();
    let storage = law
        .homogeneous_thermal_storage(&[1.0, 2.0, 1.0], thermal_controls(), || false)
        .unwrap();
    let doubled = shifted
        .homogeneous_thermal_storage(&[2.0, 4.0, 2.0], thermal_controls(), || false)
        .unwrap();
    assert_eq!(storage.sources(), law.sources());
    assert_eq!(storage.component_moles(), &[1.0, 2.0, 1.0]);
    assert_eq!(storage.caloric_identity(), None);
    assert!(storage.advected_composition().is_none());
    let mass = 0.210;
    for t in [300.0, 350.0, 400.0, 470.0, 600.0] {
        let d = t - 300.0;
        let expected_h = (520.0 * d + 0.1 * d * d) / mass;
        let actual_h = h(&storage, t);
        close(actual_h, expected_h, 2e-9);
        close(h(&doubled, t), expected_h, 2e-9);
        close(storage.temperature_k(actual_h).unwrap(), t, 2e-9);
        let cp = storage
            .specific_heat_cp_j_kg_k(Temperature::new(t), || false)
            .unwrap();
        close(cp, (460.0 + 0.2 * t) / mass, 1e-10);
        assert!(1.0 / cp <= storage.maximum_temperature_slope());
    }
    let entropy = storage
        .specific_entropy_change_j_kg_k(h(&storage, 350.0), h(&storage, 570.0))
        .unwrap()
        .unwrap();
    let expected = (460.0 * (570.0_f64 / 350.0).ln() + 0.2 * 220.0) / mass;
    close(entropy, expected, 2e-8);
    close(
        storage
            .specific_entropy_change_j_kg_k(h(&storage, 570.0), h(&storage, 350.0))
            .unwrap()
            .unwrap(),
        -expected,
        2e-8,
    );
    // The original surface evaluator admits coordinates down to 2^-200;
    // keep a trace well below unit-scale roundoff but inside that source limit.
    let trace = law
        .homogeneous_thermal_storage(&[1.0, 1e-50, 1.0], thermal_controls(), || false)
        .unwrap();
    assert_eq!(
        trace
            .caloric_state_at_specific_enthalpy(h(&trace, 480.0), || false)
            .unwrap()
            .component_moles(),
        &[1.0, 1e-50, 1.0]
    );
    assert!(matches!(
        law.homogeneous_thermal_storage(&[1.0, 1e-200, 1.0], thermal_controls(), || false),
        Err(MulticomponentLiquidCaloricError::Material(
            MulticomponentLiquidError::Source(_)
        ))
    ));
}

#[test]
fn conduction_retains_actual_composition_sources_and_exact_nonlinear_tangents() {
    let card = thermal_card(true);
    let liquid = liquid(&card);
    let law =
        MulticomponentLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0), "p").unwrap();
    for amounts in [[1.0, 2.0, 1.0], [2.0, 1.0, 1.0]] {
        let storage = law
            .homogeneous_thermal_storage(&amounts, thermal_controls(), || false)
            .unwrap();
        let conduction = MulticomponentLiquidConductionLaw::try_new(&card, &storage).unwrap();
        assert!(core::ptr::eq(conduction.caloric(), &storage));
        assert!(
            conduction
                .source_line()
                .knots
                .iter()
                .any(|(t, _)| *t == 450.0)
        );
        let selected = liquid
            .reference_state()
            .material()
            .property(THERMAL_CONDUCTIVITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected;
        assert!(
            conduction
                .source_line()
                .knots
                .iter()
                .all(|(_, answer)| answer.receipt.selected == selected)
        );
        for t in [300.0, 370.0, 420.0, 510.0, 600.0] {
            let energy = h(&storage, t);
            let state = conduction.state_at_specific_enthalpy(energy).unwrap();
            let cp = storage
                .specific_heat_cp_j_kg_k(Temperature::new(t), || false)
                .unwrap();
            let (expected_k, dk_dt) = k(t, amounts[0] / 4.0);
            close(state.temperature_k, t, 2e-9);
            close(state.temperature_slope_kg_k_j, 1.0 / cp, 2e-15);
            close(state.conductivity_w_m_k, expected_k, 2e-10);
            close(state.conductivity_slope_kg_w_m_k_j, dk_dt / cp, 2e-16);
            close(
                conduction.conductivity_w_m_k(energy).unwrap(),
                expected_k,
                2e-10,
            );
            if t > 300.0 && t < 600.0 {
                let a = conduction.state_at_specific_enthalpy(energy - 0.1).unwrap();
                let b = conduction.state_at_specific_enthalpy(energy + 0.1).unwrap();
                close(
                    (b.temperature_k - a.temperature_k) / 0.2,
                    state.temperature_slope_kg_k_j,
                    2e-10,
                );
                close(
                    (b.conductivity_w_m_k - a.conductivity_w_m_k) / 0.2,
                    state.conductivity_slope_kg_w_m_k_j,
                    2e-12,
                );
            }
        }
    }
}

fn birth(id: u64, mass_kg: f64, h: f64) -> TimedEvent {
    TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: 1,
            material_id: MaterialId(1),
            mass_kg,
            specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [id as f64, 0.0, 0.0],
                end_m: [id as f64, 1.0, 0.0],
                volume_m3: mass_kg / 1000.0,
                surface_area_m2: 1.0,
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        }),
    }
}

#[test]
fn native_implicit_heat_exchange_matches_independent_nonlinear_energy_and_entropy() {
    let card = thermal_card(true);
    let liquid = liquid(&card);
    let law =
        MulticomponentLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0), "p").unwrap();
    let storage = law
        .homogeneous_thermal_storage(&[1.0, 2.0, 1.0], thermal_controls(), || false)
        .unwrap();
    let adapter = MulticomponentLiquidConductionLaw::try_new(&card, &storage).unwrap();
    let timeline = DepositionSimulation::new(
        vec![
            birth(1, 0.210, h(&storage, 350.0)),
            birth(2, 0.420, h(&storage, 580.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.1,
                    conductance_w_k: 0.0,
                }),
            },
        ],
        &[(MaterialId(1), &storage)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 8,
        },
    )
    .unwrap();
    let (timeline, _) = timeline
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-6,
            },
            || Ok(()),
        )
        .unwrap();
    let conduction = ConductionObserver::new(
        &timeline,
        &[(MaterialId(1), &adapter as &dyn ConductivityLaw)],
        &[ConductiveLink {
            link_id: 1,
            half_distances_m: [0.1; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    let implicit = ImplicitConductionChain::new(
        &timeline,
        &conduction,
        &[(MaterialId(1), &adapter as &dyn ImplicitConductionLaw)],
        || Ok(()),
    )
    .unwrap();
    let dt = 1000.0;
    let energy = |t: f64| 520.0 * (t - 300.0) + 0.1 * (t - 300.0).powi(2);
    let temperature = |e: f64| 300.0 + 2.0 * e / (520.0 + (270400.0 + 0.4 * e).sqrt());
    let (cold, hot) = (energy(350.0), energy(580.0));
    let mut bracket = [0.0, (2.0 / 3.0) * (hot - cold)];
    for _ in 0..100 {
        let q = 0.5 * (bracket[0] + bracket[1]);
        let t = [temperature(cold + q), temperature(hot - q / 2.0)];
        let g = 1.0 / (1.0 / k(t[0], 0.25).0 + 1.0 / k(t[1], 0.25).0);
        bracket[usize::from(q - dt * g * (t[1] - t[0]) > 0.0)] = q;
    }
    let q = 0.5 * (bracket[0] + bracket[1]);
    let expected_t = [temperature(cold + q), temperature(hot - q / 2.0)];
    let config = ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 2e-6,
        total_tolerance_j: 1e-5,
        ..ImplicitConductionControls::default()
    };
    let original = timeline.cells().clone();
    let (next, report) = implicit
        .advance_to(
            &timeline,
            dt,
            config,
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    for (i, expected) in expected_t.into_iter().enumerate() {
        close(
            next.temperature_k(CellId(i as u64 + 1)).unwrap(),
            expected,
            2e-7,
        );
    }
    close(next.energy_j(), timeline.energy_j(), 1e-5);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(report.energy.external_heat_j, 0.0);
    let entropy = |a: f64, b: f64| 460.0 * (b / a).ln() + 0.2 * (b - a);
    close(
        report.entropy.thermal_storage_change_j_k,
        entropy(350.0, expected_t[0]) + 2.0 * entropy(580.0, expected_t[1]),
        2e-7,
    );
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    assert_eq!(next.mass_kg(), timeline.mass_kg());
    assert!(
        implicit
            .advance_to(
                &timeline,
                dt,
                config,
                &mut ImplicitConductionWorkspace::new(2).unwrap(),
                || Err(DepositionError::Cancelled)
            )
            .is_err()
    );
    assert_eq!(timeline.cells(), &original);
    assert_eq!(timeline.time_s(), 0.0);
}

#[test]
fn unsupported_sources_bounds_and_cancellation_refuse_before_publication() {
    let card = thermal_card(true);
    let liquid = liquid(&card);
    let law = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    assert!(
        law.homogeneous_thermal_storage(&[1.0; 3], controls(), || false)
            .is_err()
    );
    assert!(matches!(
        law.homogeneous_thermal_storage(&[1.0; 3], thermal_controls(), || true),
        Err(MulticomponentLiquidCaloricError::Cancelled)
    ));
    let storage = law
        .homogeneous_thermal_storage(&[1.0; 3], thermal_controls(), || false)
        .unwrap();
    for invalid in [
        -1.0,
        storage.specific_enthalpy_domain_j_kg()[1].next_up(),
        f64::NAN,
    ] {
        assert!(storage.temperature_k(invalid).is_err());
        assert!(
            storage
                .specific_entropy_change_j_kg_k(0.0, invalid)
                .is_err()
        );
    }
    let other = thermal_card(false);
    assert!(MulticomponentLiquidConductionLaw::try_new(&other, &storage).is_err());
    let other_liquid = super::liquid(&other);
    let other_law = MulticomponentLiquidCaloricLaw::try_new(&other_liquid, sources(), "p").unwrap();
    let unsupported = other_law
        .homogeneous_thermal_storage(&[1.0; 3], thermal_controls(), || false)
        .unwrap();
    assert!(MulticomponentLiquidConductionLaw::try_new(&other, &unsupported).is_err());
}
