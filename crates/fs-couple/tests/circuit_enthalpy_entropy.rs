//! Physical entropy on the original finite electrical/thermal history.
#[path = "circuit_enthalpy/common.rs"]
mod common;
use common::*;
use fs_couple::circuit_enthalpy::{self, *};
use fs_couple::circuit_enthalpy::entropy::{self, CircuitEntropyAudit};
use fs_thermal::enthalpy::*;
use std::cell::Cell;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:.16e} != {b:.16e}; tolerance {tolerance:e}");
}
fn balance(a: &CircuitEntropyAudit) {
    close(a.thermal_storage_change_j_k, a.reservoir_entropy_input_j_k
        + a.conduction_generation_j_k + a.electrical_loss_generation_j_k
        + a.unclassified_source_at_cell_j_k + a.temporal_defect_j_k
        + a.algebraic_residual_j_k + a.accounting_residual_j_k, 2e-14);
}
fn timeline<'a>(laws: [&'a dyn ThermalLaw; 2], enthalpies: [f64; 2], boundaries: [Boundary; 2], g: f64)
    -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (i, id) in [CellId(30), CellId(10)].into_iter().enumerate() {
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id, segment_id: i as u64, body_id: i as u64, material_id: MaterialId(i as u64),
            mass_kg: 1.0, specific_h_j_kg: enthalpies[i],
            geometry: CellGeometry { start_m: [i as f64, 0.0, 0.0],
                end_m: [i as f64+1.0, 0.0, 0.0], volume_m3: 1.0, surface_area_m2: 6.0 },
            boundary: boundaries[i],
        }) });
    }
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: 7, a: CellId(30), b: CellId(10), area_m2: 1.0, conductance_w_k: g,
    }) });
    DepositionSimulation::new(events, &[(MaterialId(0), laws[0]), (MaterialId(1), laws[1])],
        DepositionBudget { max_cells: 8, max_links: 8, max_events: 16, max_substeps: 1000 })
        .unwrap().advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-7 },
            || Ok(())).unwrap().0
}

#[test]
fn finite_capacitor_uses_original_loss_and_caloric_entropy_without_extra_circuit_trials() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb];
    let old = initial(&curves, 2.0, None);
    let electrical = circuit(&old, &cards, &curves);
    let plain = circuit_enthalpy::step_candidate(&old, 1.0, input(&electrical), config(), &[0.0; 2],
        |s| coefficients(s, &cards, &curves), || false).unwrap();
    let calls = Cell::new(0);
    let (next, audit) = entropy::step_candidate(&old, 1.0, input(&electrical), config(), &[0.0; 2],
        |s| { calls.set(calls.get()+1); coefficients(s, &cards, &curves) }, || false).unwrap();
    assert_eq!(next.circuit(), plain.circuit());
    assert_eq!(next.thermal().cells(), plain.thermal().cells());
    assert_eq!(next.report().thermal, plain.report().thermal);
    assert_eq!(calls.get(), plain.report().circuit_trials);
    assert_eq!(calls.get(), next.report().circuit_trials);
    let ta = next.thermal().temperature_k(CellId(30)).unwrap();
    let tb = next.thermal().temperature_k(CellId(10)).unwrap();
    let qa = next.report().electrical.resistors[0].heat_j;
    let qb = next.report().electrical.resistors[1].heat_j;
    let storage = (ta/300.0).ln() + 2.0*(tb/300.0).ln();
    let temporal = storage - (ta-300.0)/ta - 2.0*(tb-300.0)/tb;
    close(audit.thermal_storage_change_j_k, storage, 2e-14);
    close(audit.electrical_loss_generation_j_k, qa/ta+qb/tb, 2e-14);
    close(audit.conduction_generation_j_k, 2.0*(ta-tb).powi(2)/(ta*tb), 2e-14);
    close(audit.temporal_defect_j_k, temporal, 2e-14);
    let residual = ((ta-300.0)-qa-2.0*(tb-ta))/ta
        + (2.0*(tb-300.0)-qb-2.0*(ta-tb))/tb;
    close(audit.algebraic_residual_j_k, residual, 2e-14);
    assert!(audit.temporal_defect_j_k > 0.0);
    assert_eq!(audit.reservoir_entropy_input_j_k, 0.0);
    assert_eq!(audit.unclassified_source_at_cell_j_k, 0.0);
    assert!(!audit.has_unclassified_sources);
    close(audit.accounting_residual_j_k, 0.0, 2e-14);
    balance(&audit);
    assert_eq!(electrical.charges_c(), &[10.0]);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn real_reservoir_temperatures_and_other_sources_remain_separate_from_joule_loss() {
    let (a, ca) = material(1.0, 1.0, 0.005);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb];
    let baths = [Boundary { ambient_k: 400.0, htc_w_m2_k: 0.2, ..boundary() },
        Boundary { ambient_k: 280.0, htc_w_m2_k: 0.1, ..boundary() }];
    let old = timeline([curves[0].curve(), curves[1].curve()], [100.0, 200.0], baths, 2.0);
    let electrical = circuit(&old, &cards, &curves);
    let dt = 0.1;
    let (next, audit) = entropy::step_candidate(&old, dt, input(&electrical), config(), &[0.0; 2], |s| {
        let mut c = coefficients(s, &cards, &curves)?;
        c.sources_w.insert(CellId(30), 2.0);
        c.sources_w.insert(CellId(10), -1.0);
        Ok::<_, fs_couple::cavity_phs::piston::thermal::power::joule::JouleHeatingError>(c)
    }, || false).unwrap();
    let ta = next.thermal().temperature_k(CellId(30)).unwrap();
    let tb = next.thermal().temperature_k(CellId(10)).unwrap();
    let qa = dt*(400.0-ta);
    let qb = dt*0.5*(280.0-tb);
    close(audit.reservoir_entropy_input_j_k, qa/400.0+qb/280.0, 2e-14);
    close(audit.conduction_generation_j_k,
        dt*2.0*(ta-tb).powi(2)/(ta*tb)+qa*(1.0/ta-1.0/400.0)+qb*(1.0/tb-1.0/280.0), 2e-14);
    close(audit.unclassified_source_at_cell_j_k, dt*(2.0/ta-1.0/tb), 2e-14);
    assert!(audit.has_unclassified_sources);
    assert!(audit.electrical_loss_generation_j_k > 0.0);
    close(next.report().electrical.resistors[0].material.conductor.resistance_ohm(),
        1.0+0.005*(ta-300.0), 2e-12);
    balance(&audit);
}

#[test]
fn cancelling_independent_sources_do_not_claim_a_closed_entropy_system() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb];
    let old = timeline([curves[0].curve(), curves[1].curve()], [100.0, 200.0],
        [Boundary { source_w: 2.0, ..boundary() }, boundary()], 0.0);
    let electrical = circuit(&old, &cards, &curves);
    let (_, audit) = entropy::step_candidate(&old, 0.25, input(&electrical), config(), &[0.0; 2], |s| {
        let mut c = coefficients(s, &cards, &curves)?;
        c.sources_w.insert(CellId(30), -2.0);
        Ok::<_, fs_couple::cavity_phs::piston::thermal::power::joule::JouleHeatingError>(c)
    }, || false).unwrap();
    assert!(audit.has_unclassified_sources);
    assert_eq!(audit.unclassified_source_at_cell_j_k, 0.0);
    balance(&audit);
}

#[test]
fn physical_electrical_loss_enters_latent_storage_without_artificial_temperature_rise() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if (0.0..=1000.0).contains(&h) { Ok(300.0) }
            else { Err(DepositionError::MaterialDomain) }
        }
        fn maximum_temperature_slope(&self) -> f64 { 0.0 }
        fn specific_entropy_change_j_kg_k(&self, a: f64, b: f64) -> Result<Option<f64>, DepositionError> {
            self.temperature_k(a)?; self.temperature_k(b)?; Ok(Some((b-a)/300.0))
        }
    }
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb];
    let initial_sources = initial(&curves, 0.0, None);
    let electrical = circuit(&initial_sources, &cards, &curves);
    let old = timeline([&Plateau, &Plateau], [100.0; 2], [boundary(); 2], 2.0);
    let (next, audit) = entropy::step_candidate(&old, 1.0, input(&electrical), config(), &[0.0; 2], |_| {
        Ok::<_, &'static str>(CircuitThermalCoefficients {
            resistors: electrical.resistors().to_vec(), sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(), source: (),
        })
    }, || false).unwrap();
    let expected = next.report().deposited_joule_j/300.0;
    close(audit.thermal_storage_change_j_k, expected, 2e-10);
    close(audit.electrical_loss_generation_j_k, expected, 2e-14);
    close(audit.temporal_defect_j_k, 0.0, 2e-14);
    assert_eq!(audit.conduction_generation_j_k, 0.0);
    for id in [CellId(30), CellId(10)] { assert_eq!(next.thermal().temperature_k(id).unwrap(), 300.0); }
    balance(&audit);
}

#[test]
fn missing_entropy_and_late_entropy_cancellation_publish_neither_owner() {
    struct Caloric<'a> { cp: f64, calls: &'a AtomicUsize, supported: bool }
    impl ThermalLaw for Caloric<'_> {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            let t = 200.0+h/self.cp;
            if t.is_finite() && t > 0.0 { Ok(t) } else { Err(DepositionError::MaterialDomain) }
        }
        fn maximum_temperature_slope(&self) -> f64 { 1.0/self.cp }
        fn specific_entropy_change_j_kg_k(&self, a: f64, b: f64) -> Result<Option<f64>, DepositionError> {
            self.calls.fetch_add(1, Ordering::SeqCst);
            if !self.supported { return Ok(None); }
            let ta = self.temperature_k(a)?;
            self.temperature_k(b)?;
            Ok(Some(self.cp*((b-a)/(self.cp*ta)).ln_1p()))
        }
    }
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb];
    for supported in [false, true] {
        let calls = AtomicUsize::new(0);
        let laws = [Caloric { cp: 1.0, calls: &calls, supported }, Caloric { cp: 2.0, calls: &calls, supported }];
        let old = timeline([&laws[0], &laws[1]], [100.0, 200.0], [boundary(); 2], 2.0);
        let electrical = circuit(&old, &cards, &curves);
        let saved = old.cells().clone();
        circuit_enthalpy::step_candidate(&old, 1.0, input(&electrical), config(), &[0.0; 2],
            |s| coefficients(s, &cards, &curves), || false).unwrap();
        assert_eq!(calls.load(Ordering::SeqCst), 0, "energy-only path must not request entropy");
        let result = entropy::step_candidate(&old, 1.0, input(&electrical), config(), &[0.0; 2],
            |s| coefficients(s, &cards, &curves), || supported && calls.load(Ordering::SeqCst) != 0);
        if supported { assert!(matches!(result, Err(CircuitEnthalpyError::Cancelled))); }
        else { assert!(matches!(result, Err(CircuitEnthalpyError::Thermal(DepositionError::InvalidInput(_))))); }
        assert!(calls.load(Ordering::SeqCst) > 0);
        assert_eq!(old.cells(), &saved);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(electrical.charges_c(), &[10.0]);
    }
}

#[test]
fn time_refinement_reduces_numerical_entropy_and_preserves_total_stored_energy() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb];
    let mut defects = Vec::new();
    for count in [16, 32, 64] {
        let mut thermal = initial(&curves, 0.0, None);
        let mut electrical = circuit(&thermal, &cards, &curves);
        let energy0 = thermal.energy_j()+50.0;
        let (mut storage, mut defect) = (0.0, 0.0);
        for i in 1..=count {
            let (next, audit) = entropy::step_candidate(&thermal, i as f64/count as f64,
                input(&electrical), config(), &[0.0; 2], |s| coefficients(s, &cards, &curves), || false).unwrap();
            balance(&audit);
            storage += audit.thermal_storage_change_j_k;
            defect += audit.temporal_defect_j_k;
            close(next.thermal().energy_j()+next.report().electrical.stored_after_j, energy0, 1e-6);
            let (t, c, _) = next.into_parts(); thermal = t; electrical = c;
        }
        let expected = (thermal.temperature_k(CellId(30)).unwrap()/300.0).ln()
            + 2.0*(thermal.temperature_k(CellId(10)).unwrap()/300.0).ln();
        close(storage, expected, 1e-12);
        assert!(defect > 0.0);
        defects.push(defect);
    }
    for pair in defects.windows(2) { assert!((0.45..0.55).contains(&(pair[1]/pair[0])), "{defects:?}"); }
}
