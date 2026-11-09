//! Original isolated equilibrium -> original tangent/LU -> independent analytic
//! response and every native power row. No substitute transport or time solver.
use super::super::equilibrium::isolated::IsolatedThermalControls;
use super::super::harmonic::ThermalImpedanceBudget;
use super::super::harmonic::isolated::{
    IsolatedThermalImpedanceControls, IsolatedThermalImpedanceError,
    ThermalEnergyConstraintError,
};
use super::*;
use fs_math::c64::C64;

fn budget() -> ThermalImpedanceBudget {
    ThermalImpedanceBudget {
        maximum_cells: 16,
        maximum_dense_bytes: 100_000,
        maximum_dense_work: 100_000,
        equilibrium_tolerance_w: 1e-9,
        relative_residual_tolerance: 1e-12,
    }
}
fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:?} != {expected:?}");
}
fn unequal_pair<'a>(law: &'a Law) -> (DepositionSimulation<'a>, ConductionObserver<'a>) {
    let (mut thermal, _) = setup(law, &[300.0; 2], None, Vec::new());
    thermal.cells.get_mut(&CellId(0)).unwrap().mass_kg = 2.0;
    thermal.cells.get_mut(&CellId(1)).unwrap().mass_kg = 3.0;
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), law)],
        &[ConductiveLink {
            link_id: 0,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    ).unwrap();
    (thermal, conduction)
}

#[test]
fn closed_dc_heat_transfer_preserves_unequal_original_capacities_and_inventory() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (thermal, conduction) = unequal_pair(&law);
    let before = thermal.cells().clone();
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap();
    let equilibrium = chain.solve_isolated_equilibrium(
        &thermal, IsolatedThermalControls::default(), || Ok(()),
    ).unwrap();
    let model = equilibrium.prepare_thermal_impedance(
        budget(), IsolatedThermalImpedanceControls::default(), || Ok(()),
    ).unwrap();
    assert_eq!(model.components(), &[vec![CellId(0), CellId(1)]]);
    assert!(model.original().factor(0.0, || Ok(())).is_err());
    let factor = model.factor(0.0, || Ok(())).unwrap();
    let transfer = factor.heat_transfer_impedance(CellId(0), CellId(1), || Ok(())).unwrap();
    // G=1, C0=2, C1=3: deltaT0-deltaT1=1 and 2deltaT0+3deltaT1=0.
    complex_near(transfer[0], C64::from_re(0.6), 2e-15);
    complex_near(transfer[1], C64::from_re(-0.4), 2e-15);
    let response = factor.response(&[C64::ONE, -C64::ONE], || Ok(())).unwrap();
    complex_near(response.field.energy_j[0], C64::from_re(1.2), 3e-15);
    complex_near(response.field.energy_j[1], C64::from_re(-1.2), 3e-15);
    assert!(response.field.maximum_residual_w < 1e-14);
    assert!(response.components[0].energy_residual_j.abs() < 1e-14);
    let changed = factor.response_with_component_energy(
        &[C64::ZERO; 2], &[C64::from_re(10.0)], || Ok(()),
    ).unwrap();
    for value in changed.field.temperature_k {
        complex_near(value, C64::from_re(2.0), 2e-15);
    }
    assert!(matches!(factor.impedance_column(CellId(0), || Ok(())),
        Err(IsolatedThermalImpedanceError::Conservation(
            ThermalEnergyConstraintError::IncompatibleHeat { .. }
        ))));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(equilibrium.field().operating_time_s(), 0.0);
}

#[test]
fn closed_harmonic_response_matches_native_lu_and_has_a_regular_balanced_dc_limit() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (thermal, conduction) = unequal_pair(&law);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap();
    let equilibrium = chain.solve_isolated_equilibrium(
        &thermal, IsolatedThermalControls::default(), || Ok(()),
    ).unwrap();
    let model = equilibrium.prepare_thermal_impedance(
        budget(), IsolatedThermalImpedanceControls::default(), || Ok(()),
    ).unwrap();
    for omega in [1.0, 1e-12, 1e-100] {
        let factor = model.factor(omega, || Ok(())).unwrap();
        let response = factor.response(&[C64::ONE, -C64::ONE], || Ok(())).unwrap();
        let energy = C64::new(5.0 / 6.0, -omega).recip();
        complex_near(response.field.energy_j[0], energy, 2e-15);
        complex_near(response.field.energy_j[1], -energy, 2e-15);
        assert!(response.components[0].energy_j.abs() < 1e-14);
        assert!(response.field.maximum_residual_w < 1e-14);
    }
    let omega = 0.7;
    let factor = model.factor(omega, || Ok(())).unwrap();
    let heat = [C64::new(0.3, -0.2), C64::new(0.7, 0.1)];
    let closed = factor.response(&heat, || Ok(())).unwrap();
    let ordinary = model.original().factor(omega, || Ok(())).unwrap()
        .response(&heat, || Ok(())).unwrap();
    for (actual, expected) in closed.field.energy_j.iter().zip(ordinary.energy_j) {
        complex_near(*actual, expected, 3e-15);
    }
    complex_near(
        closed.components[0].energy_j,
        C64::new(0.1 / omega, 1.0 / omega),
        3e-15,
    );
    assert!(factor.response_with_component_energy(&heat, &[C64::ZERO], || Ok(())).is_err());
}

#[test]
fn disconnected_native_graph_keeps_separate_dc_energy_rows() {
    let law = Law { nonlinear_k: false, plateau: false };
    let events = [300.0, 400.0].into_iter().enumerate().map(|(i, h)| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(i as u64), material_id: MaterialId(0), body_id: 0, segment_id: 0,
            mass_kg: 1.0, specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [i as f64, 0.0, 0.0], end_m: [i as f64 + 1.0, 0.0, 0.0],
                volume_m3: 1.0, surface_area_m2: 2.0,
            },
            boundary: boundary(None),
        }),
    }).collect();
    let thermal = DepositionSimulation::new(
        events, &[(MaterialId(0), &law)],
        DepositionBudget { max_cells: 2, max_links: 2, max_events: 4, max_substeps: 10 },
    ).unwrap().advance_to(
        0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 0.0 }, || Ok(()),
    ).unwrap().0;
    let conduction = ConductionObserver::new(
        &thermal, &[(MaterialId(0), &law)], &[], &[], || Ok(()),
    ).unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal, &conduction, &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(), || Ok(()),
    ).unwrap();
    let equilibrium = graph.solve_isolated_equilibrium(
        &thermal, IsolatedThermalControls::default(), || Ok(()),
    ).unwrap();
    let model = equilibrium.prepare_thermal_impedance(
        budget(), IsolatedThermalImpedanceControls::default(), || Ok(()),
    ).unwrap();
    assert_eq!(model.components(), &[vec![CellId(0)], vec![CellId(1)]]);
    let factor = model.factor(0.0, || Ok(())).unwrap();
    assert!(factor.heat_transfer_impedance(CellId(0), CellId(1), || Ok(())).is_err());
    let changes = [C64::from_re(2.0), C64::from_re(-3.0)];
    let response = factor.response_with_component_energy(
        &[C64::ZERO; 2], &changes, || Ok(()),
    ).unwrap();
    assert_eq!(response.field.energy_j, changes);
    assert_eq!(response.field.temperature_k, changes);
    let harmonic = model.factor(2.0, || Ok(())).unwrap().response(
        &[C64::ONE, -C64::ONE], || Ok(()),
    ).unwrap();
    assert_eq!(harmonic.field.energy_j, [C64::new(0.0, 0.5), C64::new(0.0, -0.5)]);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cursor, 2);
}

#[test]
fn sharp_latent_storage_survives_frequency_response_without_choosing_dc_phases() {
    let law = Law { nonlinear_k: false, plateau: true };
    let (thermal, conduction) = setup(&law, &[350.0, 450.0], None, Vec::new());
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap();
    let equilibrium = chain.solve_isolated_equilibrium(
        &thermal, IsolatedThermalControls::default(), || Ok(()),
    ).unwrap();
    let model = equilibrium.prepare_thermal_impedance(
        budget(), IsolatedThermalImpedanceControls::default(), || Ok(()),
    ).unwrap();
    assert!(matches!(model.factor(0.0, || Ok(())),
        Err(IsolatedThermalImpedanceError::Conservation(
            ThermalEnergyConstraintError::Underdetermined { latent_cells: 2, .. }
        ))));
    let response = model.factor(2.0, || Ok(())).unwrap()
        .response(&[C64::ONE; 2], || Ok(())).unwrap();
    assert_eq!(response.field.temperature_k, [C64::ZERO; 2]);
    assert_eq!(response.field.energy_j, [C64::new(0.0, 0.5); 2]);
    assert_eq!(equilibrium.field().specific_enthalpies_j_kg(), &[350.0, 450.0]);
}

#[test]
fn closed_frequency_refusals_leave_the_prepared_factor_and_history_reusable() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (thermal, conduction) = unequal_pair(&law);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap();
    let equilibrium = chain.solve_isolated_equilibrium(
        &thermal, IsolatedThermalControls::default(), || Ok(()),
    ).unwrap();
    let mut too_small = budget();
    too_small.maximum_dense_bytes = 900; // Ordinary payload fits; closed metadata does not.
    assert!(equilibrium.prepare_thermal_impedance(
        too_small, IsolatedThermalImpedanceControls::default(), || Ok(()),
    ).is_err());
    let model = equilibrium.prepare_thermal_impedance(
        budget(), IsolatedThermalImpedanceControls::default(), || Ok(()),
    ).unwrap();
    let factor = model.factor(1.0, || Ok(())).unwrap();
    let expected = factor.response(&[C64::ONE; 2], || Ok(())).unwrap();
    assert!(factor.response(&[C64::ONE; 2], || Err(DepositionError::Cancelled)).is_err());
    assert!(factor.response(&[C64::from_re(f64::NAN); 2], || Ok(())).is_err());
    assert!(factor.impedance_column(CellId(99), || Ok(())).is_err());
    assert_eq!(factor.response(&[C64::ONE; 2], || Ok(())).unwrap(), expected);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 300.0);
}
