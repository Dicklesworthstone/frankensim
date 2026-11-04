//! Analytic phase data exercise the actual thermal timeline and mass transfers.

use super::*;
use crate::enthalpy::advection::{MassBoundary, MassFlow};
use crate::enthalpy::advection::reconstruction::{AxialCellStencil, LimitedEnthalpyReconstruction};
use crate::enthalpy::backward_euler::BackwardEulerStep;
use crate::enthalpy::{
    Boundary, CellGeometry, CellId, DepositionBudget, Event, StepConfig, TimedEvent,
};

const MAT: MaterialId = MaterialId(1);

struct PhaseLaw { datum: f64, identity: u8, melting_range_k: f64 }
impl ThermalLaw for PhaseLaw {
    fn caloric_identity(&self) -> Option<[u8; 32]> { Some([self.identity; 32]) }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let h = h - self.datum;
        if !h.is_finite() || !(-100.0..=300.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 0.0 { 500.0 + 0.5 * h }
            else if h <= 100.0 { 500.0 + self.melting_range_k * (h / 100.0) }
            else { 500.0 + self.melting_range_k + 0.5 * (h - 100.0) })
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.5 }
}
fn source_knots(law: &PhaseLaw) -> Vec<LiquidCarrierKnot> {
    [(-100.0, 0.0), (0.0, 0.0), (50.0, 0.5), (100.0, 1.0), (300.0, 1.0)]
        .into_iter().map(|(h, liquid)| LiquidCarrierKnot {
            specific_enthalpy_j_kg: h + law.datum,
            temperature_k: law.temperature_k(h + law.datum).unwrap(),
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: 2.0,
        }).collect()
}
fn carrier(law: &PhaseLaw) -> LiquidEnthalpyCarrier {
    LiquidEnthalpyCarrier::from_source_chart(MAT, law, &source_knots(law), || Ok(())).unwrap()
}
fn prepare<'a>(law: &'a PhaseLaw, relative_h: &[f64], density: f64) -> DepositionSimulation<'a> {
    let events = relative_h.iter().enumerate().map(|(i, &h)| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(i as u64), segment_id: i as u64, body_id: 0, material_id: MAT,
            mass_kg: 1.0, specific_h_j_kg: h + law.datum,
            geometry: CellGeometry {
                start_m: [i as f64, 0.0, 0.0], end_m: [i as f64 + 1.0, 0.0, 0.0],
                volume_m3: 1.0 / density, surface_area_m2: 1.0,
            },
            boundary: Boundary { ambient_k: 500.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }),
    }).collect();
    DepositionSimulation::new(events, &[(MAT, law as &dyn ThermalLaw)], DepositionBudget {
        max_cells: 3, max_links: 6, max_events: 3, max_substeps: 1000,
    }).unwrap().advance_to(0.0, config(1.0), || Ok(())).unwrap().0
}
fn config(dt: f64) -> StepConfig {
    StepConfig { max_dt_s: dt, energy_tolerance_j: 1e-10 }
}
fn flows(count: usize, rate: f64) -> Vec<MassFlow> {
    (0..count).map(|i| MassFlow {
        id: i as u64, a: CellId(i as u64), b: CellId(((i + 1) % count) as u64), kg_s: rate,
    }).collect()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual} != {expected}");
}

#[test]
fn isothermal_liquid_circulation_preserves_stationary_solid_fractions_and_datum() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = PhaseLaw { datum, identity: 7, melting_range_k: 0.0 };
        let initial = prepare(&law, &[20.0, 50.0, 80.0], 2.0);
        for rate in [-2.0, 2.0] {
            let flow = Advected::new((), flows(3, rate), 3).unwrap()
                .with_liquid_carrier(carrier(&law)).unwrap();
            let endpoint = flow.endpoint_transport(&initial, || Ok(())).unwrap();
            assert!(endpoint.transfers.iter().all(|transfer| transfer.power_w == 0.0));
            let (next, observed, audit) = initial.advance_observed(
                &flow, 2.0, config(10.0), || Ok(()),
            ).unwrap();
            assert_eq!(next.cells(), initial.cells());
            assert_eq!(audit.external_heat_j, 0.0);
            assert_eq!(audit.stored_energy_change_j, 0.0);
            assert_eq!(observed.audit().external_enthalpy_j, 0.0);
            assert!(observed.audit().maximum_courant <= 0.8);
            close(observed.audit().transported_mass_kg, 12.0, 1e-12);
            for &id in next.cells().keys() { assert_eq!(next.temperature_k(id).unwrap(), 500.0); }
        }
        if datum == 0.0 {
            let homogeneous = Advected::new((), flows(3, 2.0), 3).unwrap();
            let (mixed, _, _) = initial.advance_observed(&homogeneous, 0.1, config(0.1), || Ok(())).unwrap();
            assert_ne!(mixed.cells(), initial.cells());
        }
    }
}

#[test]
fn superheated_liquid_melts_stationary_solid_in_the_original_implicit_energy_balance() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let initial = prepare(&law, &[150.0, 25.0], 2.0);
    let flow = Advected::new((), flows(2, 0.5), 2).unwrap()
        .with_liquid_carrier(carrier(&law)).unwrap();
    let step = BackwardEulerStep::new(&initial, 1.0, || Ok(())).unwrap();
    // m=1, dt*q=0.5, h_l=100. The hot node loses its superheat through
    // q*(h_hot-h_l); the receiving node absorbs exactly that latent energy.
    let heat = 50.0 * 0.5 / 1.5;
    let trial = step.evaluate_with_transport(
        &[-heat, heat], |state| flow.endpoint_transport(state, || Ok(())), || Ok(()),
    ).unwrap();
    let (next, (), audit) = trial.accept(1e-10, 1e-10, || Ok(())).unwrap();
    close(next.cells()[&CellId(0)].specific_h_j_kg, 150.0 - heat, 1e-12);
    close(next.cells()[&CellId(1)].specific_h_j_kg, 25.0 + heat, 1e-12);
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 500.0);
    assert!(next.temperature_k(CellId(0)).unwrap() < initial.temperature_k(CellId(0)).unwrap());
    assert_eq!(audit.external_heat_j, 0.0);
    close(audit.stored_energy_change_j, 0.0, 1e-12);
    close(next.energy_j(), initial.energy_j(), 1e-12);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn fully_liquid_transport_preserves_the_original_trajectory_and_accounts_exactly() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let initial = prepare(&law, &[130.0, 170.0, 250.0], 2.0);
    for rate in [-0.2, 0.2] {
        let original = Advected::new((), flows(3, rate), 3).unwrap();
        let selected = original.clone().with_liquid_carrier(carrier(&law)).unwrap();
        let (ordinary, original_flow, original_energy) = initial.advance_observed(
            &original, 1.0, config(0.1), || Ok(()),
        ).unwrap();
        let (liquid, selected_flow, selected_energy) = initial.advance_observed(
            &selected, 1.0, config(0.1), || Ok(()),
        ).unwrap();
        assert_eq!(liquid.cells(), ordinary.cells());
        assert_eq!(selected_flow.audit(), original_flow.audit());
        assert_eq!(selected_energy, original_energy);
    }
}

#[test]
fn liquid_carrier_refuses_solid_flow_changed_source_and_inconsistent_inventory() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let bound = carrier(&law);
    for rate in [-0.5, 0.5] {
        let solid = prepare(&law, &[-20.0, 25.0, 75.0], 2.0);
        let observer = Advected::new((), flows(3, rate), 3).unwrap().with_liquid_carrier(bound).unwrap();
        assert!(matches!(solid.advance_observed(&observer, 0.1, config(0.1), || Ok(())),
            Err(DepositionError::InvalidInput("nonzero liquid mass flow crosses a fully solid cell"))));
        assert_eq!(observer.audit().substeps, 0);
        assert_eq!(solid.time_s(), 0.0);
    }
    let solid = prepare(&law, &[-20.0; 3], 2.0);
    let stationary = Advected::new((), flows(3, 0.0), 3).unwrap().with_liquid_carrier(bound).unwrap();
    let (same, _, _) = solid.advance_observed(&stationary, 0.1, config(0.1), || Ok(())).unwrap();
    assert_eq!(same.cells(), solid.cells());
    let changed = PhaseLaw { identity: 8, ..law };
    let wrong_source = prepare(&changed, &[25.0; 3], 2.0);
    let flowing = Advected::new((), flows(3, 0.5), 3).unwrap().with_liquid_carrier(bound).unwrap();
    assert!(matches!(flowing.endpoint_transport(&wrong_source, || Ok(())),
        Err(DepositionError::InvalidInput("liquid-carrier caloric identity or storage basis changed"))));
    let wrong_mass = prepare(&law, &[25.0; 3], 4.0);
    assert!(matches!(flowing.endpoint_transport(&wrong_mass, || Ok(())),
        Err(DepositionError::InvalidInput("liquid carriage changed original material or fixed-density inventory"))));
}

#[test]
fn phase_source_admission_and_incompatible_transport_modes_refuse() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let alloy = PhaseLaw { melting_range_k: 10.0, ..law };
    assert!(LiquidEnthalpyCarrier::from_source_chart(MAT, &alloy, &source_knots(&alloy), || Ok(())).is_err());
    let mut density = source_knots(&law);
    density[3].bulk_density_kg_m3 = 1.0;
    assert!(LiquidEnthalpyCarrier::from_source_chart(MAT, &law, &density, || Ok(())).is_err());
    let mut nonlinear_fraction = source_knots(&law);
    nonlinear_fraction[2].liquid_mass_fraction = 0.7;
    assert!(LiquidEnthalpyCarrier::from_source_chart(MAT, &law, &nonlinear_fraction, || Ok(())).is_err());
    assert!(matches!(LiquidEnthalpyCarrier::from_source_chart(MAT, &law, &source_knots(&law), || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)));

    let initial = prepare(&law, &[25.0; 3], 2.0);
    let original = Advected::new((), flows(3, 0.5), 3).unwrap();
    let selected = original.clone().with_liquid_carrier(carrier(&law)).unwrap();
    assert!(original.clone().with_entropy_audit().with_liquid_carrier(carrier(&law)).is_err());
    assert!(selected.clone().with_entropy_audit().endpoint_transport(&initial, || Ok(())).is_err());
    assert!(selected.clone().with_stationary_species(vec!["liquid".into()], vec![], 3, 1).is_err());
    let boundaries = vec![
        MassBoundary::Inlet { id: 10, cell: CellId(0), kg_s: 0.5, material_id: MAT, specific_h_j_kg: 100.0 },
        MassBoundary::Outlet { id: 11, cell: CellId(0), kg_s: 0.5 },
    ];
    assert!(selected.clone().with_boundaries(boundaries.clone(), 5).is_err());
    assert!(original.with_boundaries(boundaries, 5).unwrap().with_liquid_carrier(carrier(&law)).is_err());
    let reconstruction = LimitedEnthalpyReconstruction::try_new(&initial,
        &[AxialCellStencil { previous: CellId(0), cell: CellId(1), next: CellId(2) }], 1, || Ok(())).unwrap();
    assert!(matches!(selected.endpoint_transport_with_reconstruction(&initial, &reconstruction, || Ok(())),
        Err(DepositionError::InvalidInput("bulk-enthalpy reconstruction has no liquid-carrier phase closure"))));
}
