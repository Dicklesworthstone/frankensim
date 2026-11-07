//! Native heat exchange with independently reconstructed ternary phase response.

use super::super::{
    AMOUNTS, GAS_CONSTANT, MASSES, MOLAR_CP_AT_400, candidate, card, controls, enthalpic_model,
    enthalpy_oracle, entropy_oracle, law, phase_composition, sources,
};
use fs_material::liquid::solution::multicomponent::caloric::{
    MulticomponentLiquidCaloricLaw,
    equilibrium::storage::{
        MulticomponentLiquidEquilibriumThermalStorage,
        conduction::MulticomponentLiquidEquilibriumConductionLaw,
    },
};
use fs_material::thermal::conduction::ConductivityCurve;
use fs_qty::Temperature;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionLaw,
    ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLink, TimedEvent,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

fn prepare(
    caloric: &MulticomponentLiquidCaloricLaw<'_, '_>,
) -> MulticomponentLiquidEquilibriumThermalStorage {
    let resolved = caloric
        .liquid_law()
        .resolve_equilibrium(
            candidate(caloric.model().thermodynamics(), 400.0, AMOUNTS, true),
            || false,
        )
        .unwrap();
    let state = caloric
        .equilibrium_caloric_state(resolved, Default::default(), || false)
        .unwrap();
    caloric
        .equilibrium_thermal_storage(&state, controls(), || false)
        .unwrap()
}

fn bulk_curve() -> ConductivityCurve {
    // An explicitly authored bulk response. No phase morphology is inferred.
    ConductivityCurve::new(vec![[380.0, 1.4], [450.0, 2.8], [490.0, 2.0]]).unwrap()
}

fn bulk_k(t: f64) -> (f64, f64) {
    if t < 450.0 {
        (1.4 + 0.02 * (t - 380.0), 0.02)
    } else {
        (2.8 - 0.02 * (t - 450.0), -0.02)
    }
}

fn cp_oracle(t: f64) -> f64 {
    // Differentiate the independent common-tangent equation
    // ln((1-2q)/q) = (1600/T)(1-3q), then the actual phase mixing enthalpy.
    let (r, q) = phase_composition(t);
    let df_dq = -2.0 / r - 1.0 / q + 3.0 * 1600.0 / t;
    let dq_dt = -(1600.0 / (t * t)) * (r - q) / df_dq;
    let pure: f64 = AMOUNTS
        .iter()
        .zip(MOLAR_CP_AT_400)
        .map(|(n, cp)| n * cp)
        .sum();
    let scale = if t < 450.0 {
        0.004 * t - 0.6
    } else {
        0.006 * t - 1.5
    };
    pure * scale + 10.0 * GAS_CONSTANT * 1600.0 * (2.0 - 6.0 * q) * dq_dt
}

#[test]
fn g1_equilibrium_conduction_uses_redistribution_cp_and_one_inverse_tangent() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let storage = prepare(&caloric);
    let adapter =
        MulticomponentLiquidEquilibriumConductionLaw::try_new(&storage, bulk_curve()).unwrap();
    assert!(core::ptr::eq(adapter.caloric(), &storage));
    assert_eq!(adapter.conductivity(), &bulk_curve());
    let mass: f64 = AMOUNTS.iter().zip(MASSES).map(|(n, m)| n * m).sum();
    for t in [385.0, 410.0, 440.0, 460.0, 480.0] {
        let h = (enthalpy_oracle(t) - enthalpy_oracle(380.0)) / mass;
        let actual = adapter.state_at_specific_enthalpy(h).unwrap();
        let expected_slope = mass / cp_oracle(t);
        let (expected_k, dk_dt) = bulk_k(t);
        close(actual.temperature_k, t, 2e-7);
        close(actual.temperature_slope_kg_k_j, expected_slope, 2e-11);
        close(actual.conductivity_w_m_k, expected_k, 5e-9);
        close(
            actual.conductivity_slope_kg_w_m_k_j,
            dk_dt * expected_slope,
            5e-13,
        );
        close(adapter.conductivity_w_m_k(h).unwrap(), expected_k, 5e-9);
        let minus = adapter.state_at_specific_enthalpy(h - 1.0).unwrap();
        let plus = adapter.state_at_specific_enthalpy(h + 1.0).unwrap();
        close(
            (plus.temperature_k - minus.temperature_k) / 2.0,
            expected_slope,
            2e-8,
        );
        let state = storage
            .caloric_state_at_specific_enthalpy(h, || false)
            .unwrap();
        assert_eq!(state.candidate().phases().len(), 3);
        assert!(state.redistribution_heat_capacity_j_per_k() > 0.0);
        assert!(expected_slope < mass / state.frozen_heat_capacity_j_per_k());
        assert_eq!(
            adapter.time_observables(h).unwrap().phase_mass_fraction,
            None
        );
    }
    let incomplete = ConductivityCurve::new(vec![[390.0, 1.0], [490.0, 1.0]]).unwrap();
    assert!(MulticomponentLiquidEquilibriumConductionLaw::try_new(&storage, incomplete).is_err());
    for h in [
        -1.0,
        storage.specific_enthalpy_domain_j_kg()[1].next_up(),
        f64::NAN,
    ] {
        assert!(adapter.state_at_specific_enthalpy(h).is_err());
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
                ambient_k: 380.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        }),
    }
}

fn temperature_oracle(enthalpy_j: f64) -> f64 {
    let mut bracket = [380.0, 490.0];
    for _ in 0..70 {
        let t = 0.5 * (bracket[0] + bracket[1]);
        bracket[usize::from(enthalpy_oracle(t) > enthalpy_j)] = t;
    }
    0.5 * (bracket[0] + bracket[1])
}

#[test]
fn g1_native_heat_transport_conserves_energy_and_entropy_while_phases_redistribute() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let storage = prepare(&caloric);
    let adapter =
        MulticomponentLiquidEquilibriumConductionLaw::try_new(&storage, bulk_curve()).unwrap();
    let h = |t| {
        storage
            .specific_enthalpy_j_kg(Temperature::new(t), || false)
            .unwrap()
    };
    let mass = storage.mass_kg();
    let timeline = DepositionSimulation::new(
        vec![
            birth(1, mass, h(440.0)),
            birth(2, 2.0 * mass, h(480.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.2,
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
            half_distances_m: [0.05; 2],
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

    // Independent backward-Euler heat-exchange equation in the physical
    // transferred energy Q. Neither storage inversion nor its phase optimizer
    // is used by this scalar oracle.
    let dt = 120.0;
    let (cold, hot) = (enthalpy_oracle(440.0), enthalpy_oracle(480.0));
    let mut bracket = [0.0, (2.0 / 3.0) * (hot - cold)];
    for _ in 0..70 {
        let q = 0.5 * (bracket[0] + bracket[1]);
        let t = [
            temperature_oracle(cold + q),
            temperature_oracle(hot - q / 2.0),
        ];
        let conductance = 0.2 / (0.05 / bulk_k(t[0]).0 + 0.05 / bulk_k(t[1]).0);
        bracket[usize::from(q > dt * conductance * (t[1] - t[0]))] = q;
    }
    let q = 0.5 * (bracket[0] + bracket[1]);
    let expected = [
        temperature_oracle(cold + q),
        temperature_oracle(hot - q / 2.0),
    ];
    assert!(
        expected[0] > 450.0,
        "the cold cell must cross the original Cp and k knots"
    );
    let config = ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 5e-6,
        total_tolerance_j: 2e-5,
        ..Default::default()
    };
    let (next, report) = implicit
        .advance_to(
            &timeline,
            dt,
            config,
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    for (i, target) in expected.into_iter().enumerate() {
        close(
            next.temperature_k(CellId(i as u64 + 1)).unwrap(),
            target,
            4e-7,
        );
    }
    close(next.energy_j(), timeline.energy_j(), 2e-5);
    assert_eq!(next.mass_kg(), timeline.mass_kg());
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(report.energy.external_heat_j, 0.0);
    let expected_entropy = entropy_oracle(expected[0]) - entropy_oracle(440.0)
        + 2.0 * (entropy_oracle(expected[1]) - entropy_oracle(480.0));
    close(
        report.entropy.thermal_storage_change_j_k,
        expected_entropy,
        4e-7,
    );
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    let original_phase = storage
        .caloric_state_at_specific_enthalpy(h(440.0), || false)
        .unwrap();
    let evolved_phase = storage
        .caloric_state_at_specific_enthalpy(h(expected[0]), || false)
        .unwrap();
    assert_ne!(
        original_phase.candidate().phases()[0].total_moles(),
        evolved_phase.candidate().phases()[0].total_moles()
    );
    assert_eq!(evolved_phase.candidate().phases().len(), 3);

    let original = timeline.cells().clone();
    assert!(
        implicit
            .advance_to(
                &timeline,
                dt,
                config,
                &mut ImplicitConductionWorkspace::new(2).unwrap(),
                || Err(DepositionError::Cancelled),
            )
            .is_err()
    );
    assert_eq!(timeline.cells(), &original);
    assert_eq!(timeline.time_s(), 0.0);
}
