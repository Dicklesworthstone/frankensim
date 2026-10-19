//! Analytic caloric fixtures exercise the real shared enthalpy timeline.
use super::*;
use crate::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId,
    StepConfig, ThermalLaw, TimedEvent,
};

struct MixtureLaw {
    names: Vec<String>,
    fractions: Vec<f64>,
    reference: &'static str,
}
impl MixtureLaw {
    fn new(fractions: Vec<f64>) -> Self {
        Self { names: vec!["A".into(), "B".into()], fractions, reference: "analytic-species-reference" }
    }
    fn cp(&self) -> f64 { 2.0 * self.fractions[0] + 6.0 * self.fractions[1] }
    fn offset(&self) -> f64 { 1000.0 * self.fractions[0] - 400.0 * self.fractions[1] }
    fn h(&self, t: f64) -> f64 { self.cp() * t + self.offset() }
}
impl ThermalLaw for MixtureLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = (h - self.offset()) / self.cp();
        if t.is_finite() && t > 0.0 { Ok(t) } else { Err(DepositionError::MaterialDomain) }
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 / self.cp() }
    fn advected_composition(&self) -> Option<CaloricComposition<'_>> {
        Some(CaloricComposition { species: &self.names, mass_fractions: &self.fractions,
            enthalpy_reference: self.reference })
    }
}
fn cfg() -> StepConfig { StepConfig { max_dt_s: 0.2, energy_tolerance_j: 1e-9 } }
fn prepared<'a>(laws: &'a [MixtureLaw], initial_t: f64) -> DepositionSimulation<'a> {
    let cell = DepositedCell {
        id: CellId(0), segment_id: 0, body_id: 0, material_id: MaterialId(0),
        mass_kg: 2.0, specific_h_j_kg: laws[0].h(initial_t),
        geometry: CellGeometry { start_m: [0.0, 0.0, 0.0], end_m: [1.0, 0.0, 0.0],
            volume_m3: 1.0, surface_area_m2: 1.0 },
        boundary: Boundary { ambient_k: initial_t, htc_w_m2_k: 0.0,
            reservoir: None, source_w: 0.0 },
    };
    let materials: Vec<_> = laws.iter().enumerate()
        .map(|(i, law)| (MaterialId(i as u64), law as &dyn ThermalLaw)).collect();
    DepositionSimulation::new(
        vec![TimedEvent { time_s: 0.0, event: Event::Birth(cell) }], &materials,
        DepositionBudget { max_cells: 3, max_links: 1, max_events: 1, max_substeps: 1000 },
    ).unwrap().advance_to(0.0, cfg(), || Ok(())).unwrap().0
}
fn mixed(laws: &[MixtureLaw]) -> Advected<()> {
    Advected::new((), vec![], 0).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 10, cell: CellId(0), kg_s: 0.25,
            material_id: MaterialId(1), specific_h_j_kg: laws[1].h(400.0) },
        MassBoundary::Inlet { id: 11, cell: CellId(0), kg_s: 0.75,
            material_id: MaterialId(2), specific_h_j_kg: laws[2].h(300.0) },
        MassBoundary::Outlet { id: 12, cell: CellId(0), kg_s: 1.0 },
    ], 3).unwrap().with_stationary_species(vec!["A".into(), "B".into()], vec![
        InletSpecies { boundary_id: 10, mass_fractions: vec![1.0, 0.0] },
        InletSpecies { boundary_id: 11, mass_fractions: vec![0.0, 1.0] },
    ], 1, 2).unwrap()
}

#[test]
fn species_specific_enthalpies_share_the_original_rate_and_energy_owner() {
    let laws = [MixtureLaw::new(vec![0.25, 0.75]), MixtureLaw::new(vec![1.0, 0.0]),
        MixtureLaw::new(vec![0.0, 1.0])];
    let observer = mixed(&laws);
    assert_eq!(observer.stationary_species().unwrap().cells()[&CellId(0)], vec![0.25, 0.75]);
    let initial = prepared(&laws, 280.0);
    let balance = initial.power_balance_observed(&observer, || Ok(())).unwrap();
    // 0.25*h_A(400)+0.75*h_B(300)-h_mix(280) = 150 W.
    assert_eq!(balance.net_power_w[&CellId(0)], 150.0);
    let (next, accepted, audit) = initial.advance_observed(&observer, 0.1, cfg(), || Ok(())).unwrap();
    assert_eq!(audit.substeps, 1);
    assert!((next.temperature_k(CellId(0)).unwrap() - 281.5).abs() < 1e-12);
    assert!((next.energy_j() - initial.energy_j() - 15.0).abs() < 1e-10);
    assert!((audit.external_heat_j - 15.0).abs() < 1e-10);
    assert!((accepted.audit().external_enthalpy_j - 15.0).abs() < 1e-10);
    // Equilibrium is 310 K. Averaging supplied temperatures would give 325 K.
    let equilibrium = prepared(&laws, 310.0);
    let balance = equilibrium.power_balance_observed(&observer, || Ok(())).unwrap();
    assert_eq!(balance.net_power_w[&CellId(0)], 0.0);
    assert_eq!(balance.external_power_w, 0.0);
    assert_eq!(observer.audit().substeps, 0);
}

#[test]
fn unequal_reference_or_wrong_cell_composition_is_refused_before_any_step() {
    let mut laws = [MixtureLaw::new(vec![0.25, 0.75]), MixtureLaw::new(vec![1.0, 0.0]),
        MixtureLaw::new(vec![0.0, 1.0])];
    laws[1].reference = "different-energy-datum";
    let observer = mixed(&laws);
    let initial = prepared(&laws, 280.0);
    assert!(initial.power_balance_observed(&observer, || Ok(())).is_err());
    assert!(initial.advance_observed(&observer, 0.1, cfg(), || Ok(())).is_err());
    drop(initial);
    laws[1].reference = "analytic-species-reference";
    laws[0].fractions = vec![0.5, 0.5];
    let initial = prepared(&laws, 280.0);
    assert!(initial.power_balance_observed(&observer, || Ok(())).is_err());
}

#[test]
fn stationary_species_refuses_cycles_unbalanced_flow_and_authored_renormalization() {
    let flows = vec![
        MassFlow { id: 0, a: CellId(0), b: CellId(1), kg_s: 1.0 },
        MassFlow { id: 1, a: CellId(1), b: CellId(0), kg_s: 1.0 },
    ];
    assert!(Advected::new((), flows, 2).unwrap()
        .with_stationary_species(vec!["A".into()], vec![], 2, 1).is_err());
    let make = |outlet| Advected::new((), vec![], 0).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 0, cell: CellId(0), kg_s: 1.0,
            material_id: MaterialId(0), specific_h_j_kg: 1.0 },
        MassBoundary::Outlet { id: 1, cell: CellId(0), kg_s: outlet },
    ], 2).unwrap();
    assert!(make(2.0).with_stationary_species(vec!["A".into()],
        vec![InletSpecies { boundary_id: 0, mass_fractions: vec![1.0] }], 1, 1).is_err());
    assert!(make(1.0).with_stationary_species(vec!["A".into()],
        vec![InletSpecies { boundary_id: 0, mass_fractions: vec![0.99] }], 1, 1).is_err());
    assert!(make(1.0).with_stationary_species(vec!["A".into()], vec![], 1, 1).is_err());
}

#[test]
fn uniform_streams_retain_authored_fractions_across_signed_connections() {
    let fractions = vec![0.13, 0.87];
    let observer = Advected::new((), vec![
        MassFlow { id: 0, a: CellId(1), b: CellId(0), kg_s: -0.07 },
    ], 1).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 1, cell: CellId(0), kg_s: 0.07,
            material_id: MaterialId(0), specific_h_j_kg: 1.0 },
        MassBoundary::Outlet { id: 2, cell: CellId(1), kg_s: 0.07 },
    ], 3).unwrap().with_stationary_species(vec!["A".into(), "B".into()],
        vec![InletSpecies { boundary_id: 1, mass_fractions: fractions.clone() }], 2, 2).unwrap();
    let stationary = observer.stationary_species().unwrap();
    assert_eq!(stationary.cells()[&CellId(0)], fractions);
    assert_eq!(stationary.cells()[&CellId(1)], fractions);
    assert!(stationary.maximum_residual_kg_s() <= 1e-16);
}
