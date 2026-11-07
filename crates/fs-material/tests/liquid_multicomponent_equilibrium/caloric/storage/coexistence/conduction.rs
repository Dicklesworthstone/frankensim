//! Native spatial heat transfer into a genuine isothermal coexistence interval.

use super::{
    GAS_CONSTANT, MulticomponentLiquidCaloricLaw, card, enthalpic_model, joined, law,
    source_enthalpy, source_entropy, sources, transition_temperature,
};
use fs_material::liquid::solution::multicomponent::caloric::equilibrium::storage::coexistence::conduction::MulticomponentLiquidCoexistenceConductionLaw;
use fs_material::thermal::conduction::ConductivityCurve;
use fs_qty::Temperature;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionLaw,
    ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

fn bulk_k(t: f64) -> f64 {
    1.0 + 0.01 * (t - 380.0)
}

fn bulk_curve() -> ConductivityCurve {
    ConductivityCurve::new(vec![[380.0, bulk_k(380.0)], [595.0, bulk_k(595.0)]]).unwrap()
}

#[test]
fn g1_coexistence_conduction_keeps_the_exact_latent_tangent_and_transition_fraction() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let storage = joined(&caloric);
    let adapter =
        MulticomponentLiquidCoexistenceConductionLaw::try_new(&storage, bulk_curve()).unwrap();
    assert!(core::ptr::eq(adapter.caloric(), &storage));
    assert_eq!(adapter.conductivity(), &bulk_curve());
    let transition = transition_temperature();
    let [_, latent] = storage.coexistence_specific_enthalpy_domain_j_kg();
    for fraction in [0.0, 0.2, 0.8, 1.0] {
        let h = fraction * latent;
        let state = adapter.state_at_specific_enthalpy(h).unwrap();
        assert_eq!(state.temperature_k, transition);
        assert_eq!(state.temperature_slope_kg_k_j, 0.0);
        assert_eq!(state.conductivity_slope_kg_w_m_k_j, 0.0);
        close(state.conductivity_w_m_k, bulk_k(transition), 1e-14);
        let observed = adapter.time_observables(h).unwrap();
        assert_eq!(observed.temperature_k, transition);
        close(
            observed.phase_mass_fraction.unwrap(),
            fraction,
            4.0 * f64::EPSILON,
        );
    }
    let h = (source_enthalpy(590.0) - source_enthalpy(transition)) / storage.mass_kg();
    let state = adapter.state_at_specific_enthalpy(h).unwrap();
    let slope = storage.mass_kg() / (220.0 * (0.006 * 590.0 - 1.5));
    close(state.temperature_k, 590.0, 2e-7);
    close(state.temperature_slope_kg_k_j, slope, 2e-12);
    close(state.conductivity_slope_kg_w_m_k_j, 0.01 * slope, 2e-14);
    assert_eq!(
        adapter.time_observables(h).unwrap().phase_mass_fraction,
        Some(1.0)
    );
    let cold = storage
        .specific_enthalpy_j_kg(Temperature::new(400.0), || false)
        .unwrap();
    assert_eq!(
        adapter.time_observables(cold).unwrap().phase_mass_fraction,
        Some(0.0)
    );
    let incomplete = ConductivityCurve::new(vec![[381.0, 1.0], [595.0, 1.0]]).unwrap();
    assert!(MulticomponentLiquidCoexistenceConductionLaw::try_new(&storage, incomplete).is_err());
    let [lo, hi] = storage.specific_enthalpy_domain_j_kg();
    for h in [lo.next_down(), hi.next_up(), f64::NAN] {
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

fn hot_temperature(enthalpy_j: f64) -> f64 {
    let mut bracket = [transition_temperature(), 595.0];
    for _ in 0..70 {
        let t = 0.5 * (bracket[0] + bracket[1]);
        bracket[usize::from(source_enthalpy(t) > enthalpy_j)] = t;
    }
    0.5 * (bracket[0] + bracket[1])
}

#[test]
fn g1_native_heat_transfer_changes_phase_inventory_at_constant_coexistence_temperature() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let storage = joined(&caloric);
    let adapter =
        MulticomponentLiquidCoexistenceConductionLaw::try_new(&storage, bulk_curve()).unwrap();
    let mass = storage.mass_kg();
    let transition = transition_temperature();
    let latent_j = 400.0 * GAS_CONSTANT;
    let [_, latent] = storage.coexistence_specific_enthalpy_domain_j_kg();
    let initial_cold = 0.2 * latent;
    let initial_hot = storage
        .specific_enthalpy_j_kg(Temperature::new(590.0), || false)
        .unwrap();
    let timeline = DepositionSimulation::new(
        vec![
            birth(1, mass, initial_cold),
            birth(2, mass, initial_hot),
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

    // Independent backward Euler balance in transferred heat Q: the cold
    // temperature is fixed by coexistence, while the hot original Cp law cools.
    let dt = 20.0;
    let hot_energy = source_enthalpy(590.0);
    let mut bracket = [0.0, 0.8 * latent_j];
    for _ in 0..70 {
        let q = 0.5 * (bracket[0] + bracket[1]);
        let hot = hot_temperature(hot_energy - q);
        let conductance = 0.2 / (0.05 / bulk_k(transition) + 0.05 / bulk_k(hot));
        bracket[usize::from(q > dt * conductance * (hot - transition))] = q;
    }
    let q = 0.5 * (bracket[0] + bracket[1]);
    assert!(q > 0.0 && q < 0.8 * latent_j);
    let expected_hot = hot_temperature(hot_energy - q);
    let expected_fraction = 0.2 + q / latent_j;
    let (next, report) = implicit
        .advance_to(
            &timeline,
            dt,
            ImplicitConductionControls {
                maximum_dt_s: dt,
                cell_tolerance_j: 5e-6,
                total_tolerance_j: 2e-5,
                ..Default::default()
            },
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), transition);
    close(next.temperature_k(CellId(2)).unwrap(), expected_hot, 3e-7);
    close(next.energy_j(), timeline.energy_j(), 2e-5);
    assert_eq!(next.mass_kg(), timeline.mass_kg());
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(report.energy.external_heat_j, 0.0);
    let final_cold = next.cells()[&CellId(1)].specific_h_j_kg;
    let final_hot = next.cells()[&CellId(2)].specific_h_j_kg;
    close(mass * (final_cold - initial_cold), q, 8e-6);
    close(mass * (initial_hot - final_hot), q, 8e-6);
    let state = storage
        .caloric_state_at_specific_enthalpy(final_cold, || false)
        .unwrap();
    close(state.upper_branch_mass_fraction(), expected_fraction, 3e-9);
    let phases = state.coexistence_state().unwrap();
    assert_eq!(phases.phases().count(), 4);
    let mut conserved = [0.0; 3];
    for phase in phases.phases() {
        for (sum, amount) in conserved.iter_mut().zip(phase.component_moles()) {
            *sum += amount;
        }
    }
    for amount in conserved {
        close(amount, 1.0, 2e-12);
    }
    let expected_entropy = q / transition + source_entropy(expected_hot) - source_entropy(590.0);
    close(
        report.entropy.thermal_storage_change_j_k,
        expected_entropy,
        3e-8,
    );
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    let original_entropy = mass
        * (storage
            .specific_source_entropy_change_j_kg_k(initial_cold, final_cold, || false)
            .unwrap()
            + storage
                .specific_source_entropy_change_j_kg_k(initial_hot, final_hot, || false)
                .unwrap());
    let heat_integral = mass
        * (storage
            .specific_entropy_change_j_kg_k(initial_cold, final_cold)
            .unwrap()
            .unwrap()
            + storage
                .specific_entropy_change_j_kg_k(initial_hot, final_hot)
                .unwrap()
                .unwrap());
    let admitted_residual = (state.upper_branch_mass_fraction() - 0.2)
        * storage.coexistence().latent_clausius_residual_j_per_k();
    close(original_entropy - heat_integral, admitted_residual, 2e-12);
}
