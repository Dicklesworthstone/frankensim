//! Analytic caloric fixtures exercise the real shared enthalpy timeline.
use super::*;
use crate::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId,
    StepConfig, ThermalLaw, TimedEvent,
};

#[test]
fn generated_simplex_closes_roundoff_oscillation_without_distorting_trace_species() {
    // Actual generated reactor states: correcting only the largest entry
    // alternates between adjacent sums on either side of one.
    for original in [
        [0.6237127367872349, 0.126287263212765, 0.24999999999999992],
        [0.8623144793038832, 0.12469282590474756, 0.012992694791369117],
    ] {
        for trace in [0.0, 1e-30] {
            let mut fractions = original.to_vec();
            fractions.extend([trace, 0.0]);
            let before = fractions.clone();
            close_generated_simplex(&mut fractions).unwrap();
            assert_eq!(fractions.iter().sum::<f64>(), 1.0);
            assert_eq!(fractions[3].to_bits(), trace.to_bits());
            assert_eq!(fractions[4], 0.0);
            let tolerance = 32.0 * f64::EPSILON * fractions.len() as f64;
            for (&old, &new) in before.iter().zip(&fractions) {
                assert!((new - old).abs() <= tolerance * old);
                assert_eq!(new > 0.0, old > 0.0);
            }
            let closed = fractions.clone();
            close_generated_simplex(&mut fractions).unwrap();
            assert_eq!(fractions, closed);
        }
    }
    assert!(close_generated_simplex(&mut [0.2, 0.2, 0.2]).is_err());
}

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

fn selective_stream(selected_rate: f64) -> Advected<()> {
    Advected::new((), vec![], 0).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 10, cell: CellId(0), kg_s: 1.0,
            material_id: MaterialId(1), specific_h_j_kg: 1900.0 },
        MassBoundary::Outlet { id: 11, cell: CellId(0), kg_s: 1.0 - selected_rate },
        MassBoundary::Outlet { id: 12, cell: CellId(0), kg_s: selected_rate },
    ], 3).unwrap().with_selective_outlets(vec![selected_outlet(12)], 2).unwrap()
}

fn selected_outlet(boundary_id: u64) -> SelectiveOutlet {
    SelectiveOutlet { boundary_id, mass_fractions: vec![0.0, 1.0],
        specific_h_j_kg: -2000.0, enthalpy_reference: "analytic-species-reference".into() }
}

fn solve_selective(observer: Advected<()>) -> Result<Advected<()>, DepositionError> {
    observer.with_stationary_species(vec!["A".into(), "B".into()],
        vec![InletSpecies { boundary_id: 10, mass_fractions: vec![0.5, 0.5] }], 1, 2)
}

#[test]
fn selective_species_and_phase_enthalpy_use_the_original_thermal_account() {
    let observer = solve_selective(selective_stream(0.25)).unwrap();
    let fractions = observer.stationary_species().unwrap().cells()[&CellId(0)].clone();
    assert_eq!(fractions, vec![2.0 / 3.0, 1.0 / 3.0]);
    for (index, incoming) in [0.5, 0.5].into_iter().enumerate() {
        let selected = if index == 1 { 0.25 } else { 0.0 };
        assert!((incoming - 0.75 * fractions[index] - selected).abs() < 1e-15);
    }
    let laws = [MixtureLaw::new(fractions), MixtureLaw::new(vec![0.5, 0.5])];
    let initial = prepared(&laws, 280.0);
    // The carried liquid/selected-phase h is outside the gas inversion's
    // domain. Its explicit common reference still makes the flux meaningful.
    assert!(laws[1].temperature_k(-2000.0).is_err());
    let power = initial.power_balance_observed(&observer, || Ok(())).unwrap();
    // 1*h_in - 0.75*h_remaining(280 K) - 0.25*h_selected = 1300 W.
    assert!((power.net_power_w[&CellId(0)] - 1300.0).abs() < 1e-10);
    let (next, accepted, audit) = initial.advance_observed(&observer, 0.1, cfg(), || Ok(())).unwrap();
    assert_eq!(audit.substeps, 1);
    assert!((next.energy_j() - initial.energy_j() - 130.0).abs() < 1e-9);
    assert!((audit.external_heat_j - 130.0).abs() < 1e-9);
    assert!((accepted.audit().external_enthalpy_j - 130.0).abs() < 1e-9);
    assert_eq!(accepted.audit().inlet_mass_kg, 0.1);
    assert_eq!(accepted.audit().outlet_mass_kg, 0.1);
    assert!(accepted.audit().maximum_courant <= 0.8);
    assert!(initial.power_balance_observed(&observer.clone().with_entropy_audit(), || Ok(())).is_err());

    let mut wrong = selected_outlet(12);
    wrong.enthalpy_reference = "unrelated-liquid-zero".into();
    let wrong = solve_selective(selective_stream(0.25).with_selective_outlets(vec![wrong], 2).unwrap()).unwrap();
    assert!(initial.advance_observed(&wrong, 0.1, cfg(), || Ok(())).is_err());
    assert_eq!(observer.audit().substeps, 0);
}

#[test]
fn selective_streams_propagate_the_remaining_composition_downstream() {
    let observer = Advected::new((), vec![
        MassFlow { id: 0, a: CellId(0), b: CellId(1), kg_s: 0.75 },
    ], 1).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 10, cell: CellId(0), kg_s: 1.0,
            material_id: MaterialId(0), specific_h_j_kg: 1.0 },
        MassBoundary::Outlet { id: 11, cell: CellId(1), kg_s: 0.625 },
        MassBoundary::Outlet { id: 12, cell: CellId(0), kg_s: 0.25 },
        MassBoundary::Outlet { id: 13, cell: CellId(1), kg_s: 0.125 },
    ], 5).unwrap().with_selective_outlets(vec![selected_outlet(13), selected_outlet(12)], 2).unwrap()
        .with_stationary_species(vec!["A".into(), "B".into()],
            vec![InletSpecies { boundary_id: 10, mass_fractions: vec![0.5, 0.5] }], 2, 2).unwrap();
    let solved = observer.stationary_species().unwrap();
    assert_eq!(solved.cells()[&CellId(0)], vec![2.0 / 3.0, 1.0 / 3.0]);
    assert_eq!(solved.cells()[&CellId(1)], vec![0.8, 0.2]);
    assert!(solved.maximum_residual_kg_s() < 1e-14);
}

#[test]
fn selective_outlets_refuse_overdraw_unknown_bindings_and_stale_species() {
    assert!(solve_selective(selective_stream(0.75)).is_err());
    let valid = solve_selective(selective_stream(0.25)).unwrap();
    for id in [10, 99] {
        assert!(valid.clone().with_selective_outlets(vec![selected_outlet(id)], 2).is_err());
    }
    assert!(valid.clone().with_selective_outlets(vec![selected_outlet(12), selected_outlet(12)], 2).is_err());
    let mut bad = selected_outlet(12);
    bad.mass_fractions = vec![0.0, 0.99];
    assert!(valid.clone().with_selective_outlets(vec![bad], 2).is_err());
    let mut bad = selected_outlet(12);
    bad.specific_h_j_kg = f64::NAN;
    assert!(valid.clone().with_selective_outlets(vec![bad], 2).is_err());
    let stale = valid.with_selective_outlets(vec![selected_outlet(12)], 2).unwrap();
    let laws = [MixtureLaw::new(vec![2.0 / 3.0, 1.0 / 3.0]), MixtureLaw::new(vec![0.5, 0.5])];
    assert!(prepared(&laws, 280.0).power_balance_observed(&stale, || Ok(())).is_err());
}

#[test]
fn selective_outlets_refuse_closed_bulk_flow_and_recirculation() {
    let no_bulk = Advected::new((), vec![], 0).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 10, cell: CellId(0), kg_s: 1.0,
            material_id: MaterialId(0), specific_h_j_kg: 1.0 },
        MassBoundary::Outlet { id: 12, cell: CellId(0), kg_s: 1.0 },
    ], 2).unwrap().with_selective_outlets(vec![SelectiveOutlet {
        mass_fractions: vec![0.5, 0.5], ..selected_outlet(12)
    }], 2).unwrap();
    assert!(solve_selective(no_bulk).is_err());
    let cycle = recirculating().with_selective_outlets(vec![selected_outlet(4)], 2).unwrap()
        .with_stationary_species(vec!["A".into(), "B".into()], vec![
            InletSpecies { boundary_id: 2, mass_fractions: vec![1.0, 0.0] },
            InletSpecies { boundary_id: 3, mass_fractions: vec![0.0, 1.0] },
        ], 2, 2);
    assert!(matches!(cycle, Err(DepositionError::InvalidInput(
        "selective outlets require an acyclic stationary flow graph"))));
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
fn stationary_species_refuses_unfed_cycles_unbalanced_flow_and_authored_renormalization() {
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

fn recirculating() -> Advected<()> {
    Advected::new((), vec![
        MassFlow { id: 0, a: CellId(0), b: CellId(1), kg_s: 2.0 },
        MassFlow { id: 1, a: CellId(1), b: CellId(0), kg_s: 1.0 },
    ], 2).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 2, cell: CellId(0), kg_s: 1.0,
            material_id: MaterialId(2), specific_h_j_kg: 1800.0 },
        MassBoundary::Inlet { id: 3, cell: CellId(1), kg_s: 1.0,
            material_id: MaterialId(3), specific_h_j_kg: 1400.0 },
        MassBoundary::Outlet { id: 4, cell: CellId(1), kg_s: 2.0 },
    ], 5).unwrap()
}

#[test]
fn g1_inlet_fed_recirculation_solves_actual_multispecies_balances() {
    // Fresh inlet 0 carries half A, half B; inlet 1 carries C. Recycling one
    // of three outgoing kg/s makes cell 0 75% first inlet and cell 1 50%.
    let observer = recirculating().with_stationary_species(
        vec!["A".into(), "B".into(), "C".into()], vec![
            InletSpecies { boundary_id: 2, mass_fractions: vec![0.5, 0.5, 0.0] },
            InletSpecies { boundary_id: 3, mass_fractions: vec![0.0, 0.0, 1.0] },
        ], 2, 3).unwrap();
    let solved = observer.stationary_species().unwrap();
    for (id, expected) in [(CellId(0), [0.375, 0.375, 0.25]), (CellId(1), [0.25, 0.25, 0.5])] {
        for (&actual, reference) in solved.cells()[&id].iter().zip(expected) {
            assert!((actual - reference).abs() < 1e-14);
        }
    }
    assert!(solved.maximum_residual_kg_s() < 1e-14);
}

#[test]
fn g3_disconnected_uniform_cycles_preserve_each_inlet_and_ignore_zero_flow() {
    let first = vec![0.13, 0.87];
    let second = vec![0.3, 0.7];
    let mut flows = Vec::new();
    let mut boundaries = Vec::new();
    let mut inlets = Vec::new();
    for (base, fractions) in [(0, &first), (10, &second)] {
        flows.push(MassFlow { id: base, a: CellId(base + 1), b: CellId(base), kg_s: -2.0 });
        flows.push(MassFlow { id: base + 1, a: CellId(base + 1), b: CellId(base), kg_s: 1.0 });
        boundaries.push(MassBoundary::Inlet { id: base + 2, cell: CellId(base), kg_s: 1.0,
            material_id: MaterialId(base), specific_h_j_kg: 1.0 });
        boundaries.push(MassBoundary::Outlet { id: base + 3, cell: CellId(base + 1), kg_s: 1.0 });
        inlets.push(InletSpecies { boundary_id: base + 2, mass_fractions: fractions.clone() });
    }
    // As on the original DAG path, zero rates do not introduce a flowing cell
    // or supply a composition boundary condition to another component.
    flows.push(MassFlow { id: 20, a: CellId(0), b: CellId(90), kg_s: 0.0 });
    flows.reverse();
    boundaries.reverse();
    let observer = Advected::new((), flows, 5).unwrap().with_boundaries(boundaries, 9).unwrap()
        .with_stationary_species(vec!["A".into(), "B".into()], inlets, 4, 2).unwrap();
    let solved = observer.stationary_species().unwrap();
    assert_eq!(solved.cells().len(), 4);
    for id in [CellId(0), CellId(1)] { assert_eq!(solved.cells()[&id], first); }
    for id in [CellId(10), CellId(11)] { assert_eq!(solved.cells()[&id], second); }
}

#[test]
fn g1_cycle_between_acyclic_prefix_and_tail_uses_the_same_species_owner() {
    let observer = Advected::new((), vec![
        MassFlow { id: 0, a: CellId(0), b: CellId(1), kg_s: 1.0 },
        MassFlow { id: 1, a: CellId(1), b: CellId(2), kg_s: 2.0 },
        MassFlow { id: 2, a: CellId(2), b: CellId(1), kg_s: 1.0 },
        MassFlow { id: 3, a: CellId(2), b: CellId(3), kg_s: 2.0 },
    ], 4).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 4, cell: CellId(0), kg_s: 1.0,
            material_id: MaterialId(0), specific_h_j_kg: 1.0 },
        MassBoundary::Inlet { id: 5, cell: CellId(2), kg_s: 1.0,
            material_id: MaterialId(1), specific_h_j_kg: 1.0 },
        MassBoundary::Outlet { id: 6, cell: CellId(3), kg_s: 2.0 },
    ], 7).unwrap().with_stationary_species(vec!["A".into(), "B".into()], vec![
        InletSpecies { boundary_id: 4, mass_fractions: vec![1.0, 0.0] },
        InletSpecies { boundary_id: 5, mass_fractions: vec![0.0, 1.0] },
    ], 4, 2).unwrap();
    let solved = observer.stationary_species().unwrap();
    assert_eq!(solved.cells()[&CellId(0)], vec![1.0, 0.0]);
    for (id, expected) in [(CellId(1), [0.75, 0.25]),
        (CellId(2), [0.5, 0.5]), (CellId(3), [0.5, 0.5])] {
        for (&actual, reference) in solved.cells()[&id].iter().zip(expected) {
            assert!((actual - reference).abs() < 1e-14);
        }
    }
}

#[test]
fn disconnected_closed_cycle_cannot_inherit_an_unrelated_inlet() {
    let flows = vec![
        MassFlow { id: 10, a: CellId(10), b: CellId(11), kg_s: 1.0 },
        MassFlow { id: 11, a: CellId(11), b: CellId(10), kg_s: 1.0 },
        MassFlow { id: 12, a: CellId(0), b: CellId(10), kg_s: 0.0 },
    ];
    let result = Advected::new((), flows, 3).unwrap().with_boundaries(vec![
        MassBoundary::Inlet { id: 0, cell: CellId(0), kg_s: 1.0,
            material_id: MaterialId(0), specific_h_j_kg: 1.0 },
        MassBoundary::Outlet { id: 1, cell: CellId(0), kg_s: 1.0 },
    ], 5).unwrap().with_stationary_species(vec!["A".into()],
        vec![InletSpecies { boundary_id: 0, mass_fractions: vec![1.0] }], 3, 1);
    assert!(matches!(result, Err(DepositionError::InvalidInput(
        "stationary species component has no inlet"))));
}

#[test]
fn dense_cycle_budget_and_nearly_sealed_conditioning_refuse() {
    let ring = |n: usize, feed: f64| {
        let flows = (0..n).map(|i| MassFlow { id: i as u64,
            a: CellId(i as u64), b: CellId(((i + 1) % n) as u64), kg_s: 1.0 }).collect();
        Advected::new((), flows, n).unwrap().with_boundaries(vec![
            MassBoundary::Inlet { id: n as u64, cell: CellId(0), kg_s: feed,
                material_id: MaterialId(0), specific_h_j_kg: 1.0 },
            MassBoundary::Outlet { id: n as u64 + 1, cell: CellId(0), kg_s: feed },
        ], n + 2).unwrap().with_stationary_species(vec!["A".into(), "B".into()],
            vec![InletSpecies { boundary_id: n as u64, mass_fractions: vec![0.5, 0.5] }], n, 2)
    };
    assert!(matches!(ring(MAX_RECIRCULATION_CELLS + 1, 1.0),
        Err(DepositionError::Budget("stationary recirculation cells"))));
    assert!(matches!(ring(2, 1e-12), Err(DepositionError::InvalidInput(_))));
    assert!(ring(2, 0.01).is_ok());
}

#[test]
fn g1_recirculating_mixtures_advance_original_enthalpy_and_external_work_once() {
    let observer = recirculating().with_stationary_species(vec!["A".into(), "B".into()], vec![
        InletSpecies { boundary_id: 2, mass_fractions: vec![1.0, 0.0] },
        InletSpecies { boundary_id: 3, mass_fractions: vec![0.0, 1.0] },
    ], 2, 2).unwrap();
    // Consumers bind their caloric laws to the ACTUAL generated simplex, not
    // independently rounded nominal fractions with a different bit identity.
    let compositions = observer.stationary_species().unwrap().cells();
    let laws = [MixtureLaw::new(compositions[&CellId(0)].clone()),
        MixtureLaw::new(compositions[&CellId(1)].clone()),
        MixtureLaw::new(vec![1.0, 0.0]), MixtureLaw::new(vec![0.0, 1.0])];
    let prepare = |temperatures: [f64; 2]| {
        let events = temperatures.iter().enumerate().map(|(i, &t)| TimedEvent {
            time_s: 0.0, event: Event::Birth(DepositedCell {
                id: CellId(i as u64), segment_id: i as u64, body_id: 0,
                material_id: MaterialId(i as u64), mass_kg: 2.0, specific_h_j_kg: laws[i].h(t),
                geometry: CellGeometry { start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0], volume_m3: 1.0, surface_area_m2: 1.0 },
                boundary: Boundary { ambient_k: t, htc_w_m2_k: 0.0,
                    reservoir: None, source_w: 0.0 },
            }) }).collect();
        let materials: Vec<_> = laws.iter().enumerate()
            .map(|(i, law)| (MaterialId(i as u64), law as &dyn ThermalLaw)).collect();
        DepositionSimulation::new(events, &materials, DepositionBudget {
            max_cells: 4, max_links: 2, max_events: 2, max_substeps: 100,
        }).unwrap().advance_to(0.0, cfg(), || Ok(())).unwrap().0
    };
    let initial = prepare([300.0, 300.0]);
    let balance = initial.power_balance_observed(&observer, || Ok(())).unwrap();
    assert!((balance.net_power_w[&CellId(0)] - 200.0).abs() < 1e-10);
    assert!(balance.net_power_w[&CellId(1)].abs() < 1e-10);
    assert!((balance.external_power_w - 200.0).abs() < 1e-10);
    let (next, accepted, audit) = initial.advance_observed(&observer, 0.1, cfg(), || Ok(())).unwrap();
    assert_eq!(audit.substeps, 1);
    assert!((next.temperature_k(CellId(0)).unwrap() - (300.0 + 10.0 / 3.0)).abs() < 1e-12);
    assert!((next.temperature_k(CellId(1)).unwrap() - 300.0).abs() < 1e-12);
    assert!((next.energy_j() - initial.energy_j() - 20.0).abs() < 1e-10);
    assert!((audit.external_heat_j - 20.0).abs() < 1e-10);
    assert!((accepted.audit().external_enthalpy_j - 20.0).abs() < 1e-10);
    assert_eq!(observer.audit().substeps, 0);
    // Two different mixture caloric laws reach different temperatures. Using
    // mixture-blind temperature averaging or a second heat integrator cannot
    // satisfy both local balances and this common enthalpy datum.
    let equilibrium = prepare([350.0, 325.0]);
    let balance = equilibrium.power_balance_observed(&observer, || Ok(())).unwrap();
    assert!(balance.net_power_w[&CellId(0)].abs() < 1e-10);
    assert!(balance.net_power_w[&CellId(1)].abs() < 1e-10);
    assert!(balance.external_power_w.abs() < 1e-10);
}
