//! G0/G1/G3/G4 deposition regressions; synthetic SI fixtures, not filament data.
mod support;

use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::deposition::implicit::*;
use fs_conduction::{ResistanceUncertainty, SeriesThermalResistance, ThermalResistanceTerm};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use support::{with_cancelled_cx, with_cx};

fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([0x41; 32]), SolidLiquidPhase::Solid,
        vec![knot(0.0, 200.0, 0.0), knot(600_000.0, 800.0, 0.0)]).unwrap()
}
fn knot(h: f64, t: f64, f: f64) -> EnthalpyPhaseKnot {
    EnthalpyPhaseKnot { specific_enthalpy_j_kg: h, temperature_k: t,
        liquid_mass_fraction: f, bulk_density_kg_m3: 1000.0 }
}
fn budget() -> DepositionBudget {
    DepositionBudget { max_cells: 1024, max_contacts: 4096, max_materials: 4,
        max_events: 1024, max_steps: 10000, max_operator_solves: 1_000_000 }
}
fn config(dt: f64) -> DepositionStepConfig {
    DepositionStepConfig { maximum_step_s: dt, solver_energy_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-7 }
}
fn simulation(curve: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_, ()> {
    DepositionSimulation::new(vec![(MaterialId(7), DepositionMaterial {
        curve, conductivity_w_m_k: 10.0, emissivity: 0.0 })], budget()).unwrap()
}
fn birth(id: u64, temperature: f64, convection: f64) -> MaterialBirth<()> {
    MaterialBirth { id: CellId(id), body: BodyId(1), segment: SegmentId(id), material: MaterialId(7),
        mass_kg: 0.001, specific_enthalpy_j_kg: (temperature - 200.0) * 1000.0,
        geometry: CellGeometry { position_m: [0.0, 0.0, id as f64 * 1e-3],
            material_axes: [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]],
            volume_m3: 1e-6, surface_area_m2: 1e-3, detail: () },
        environment: Environment { ambient_k: 300.0, radiation_k: 300.0,
            convection_w_m2_k: convection, source_w: 0.0 } }
}
fn event(id: u64, time: f64, births: Vec<MaterialBirth<()>>) -> DepositionEvent<()> {
    DepositionEvent { id: EventId(id), time_s: time, births, contacts: vec![], environments: vec![] }
}
fn contact(id: u64, a: u64, b: u64, r: f64) -> ContactChange {
    ContactChange::Set(ThermalContact { id: InterfaceId(id), a: CellId(a), b: CellId(b), area_m2: 1e-4,
        resistance: SeriesThermalResistance::new(vec![ThermalResistanceTerm::declared("synthetic",
            r, ResistanceUncertainty::Unstated, "analytic regression; not measured adhesion data").unwrap()]).unwrap() })
}
fn temperature(sim: &DepositionSimulation<'_, ()>, id: u64) -> f64 {
    sim.cells()[&CellId(id)].phase.temperature_k()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual-expected).abs() <= tolerance, "actual={actual:e} expected={expected:e}");
}

#[test]
fn g0_birth_preserves_old_history_and_counts_incoming_energy_once() {
    let law = curve(); let mut sim = simulation(&law);
    with_cx(|cx| {
        let first = event(1, 0.0, vec![birth(1, 400.0, 0.0)]);
        let r = sim.apply_event(cx, &first, 1e-9).unwrap();
        close(r.incoming_mass_kg, 0.001, 1e-15); close(r.incoming_enthalpy_j, 200.0, 1e-12);
        let old = sim.cells()[&CellId(1)].phase;
        sim.apply_event(cx, &event(2, 0.0, vec![birth(2, 600.0, 0.0)]), 1e-9).unwrap();
        assert_eq!(old, sim.cells()[&CellId(1)].phase);
        close(sim.total_mass_kg(), 0.002, 1e-15); close(sim.stored_enthalpy_j(), 600.0, 1e-9);
        assert!(sim.apply_event(cx, &first, 1e-9).is_err());
        close(sim.total_mass_kg(), 0.002, 1e-15);
    });
}

#[test]
fn g1_single_bead_matches_backward_euler_and_energy_balance() {
    let law = curve(); let mut sim = simulation(&law);
    with_cx(|cx| {
        sim.apply_event(cx, &event(1, 0.0, vec![birth(1, 500.0, 10.0)]), 1e-9).unwrap();
        let report = sim.advance_to(cx, 1.0, config(1.0)).unwrap();
        // C = m Cp = 1 J/K, h A = 0.01 W/K, T_new=(500+0.01*300)/1.01.
        close(temperature(&sim, 1), 300.0 + 200.0 / 1.01, 1e-8);
        close(report.stored_energy_change_j, report.boundary_heat_j, 1e-8);
        assert!(report.boundary_heat_j < 0.0); assert_eq!(report.steps, 1);
        assert_eq!(sim.cells()[&CellId(1)].reference_temperature_k, 500.0);
    });
}

#[test]
fn g1_fixed_domain_time_refinement_approaches_exponential() {
    let law = curve();
    let errors: Vec<_> = [4.0, 2.0, 1.0].into_iter().map(|dt| {
        let mut sim = simulation(&law);
        with_cx(|cx| {
            sim.apply_event(cx, &event(1, 0.0, vec![birth(1, 500.0, 10.0)]), 1e-9).unwrap();
            sim.advance_to(cx, 20.0, config(dt)).unwrap();
        });
        (temperature(&sim, 1) - (300.0 + 200.0 * (-0.2_f64).exp())).abs()
    }).collect();
    assert!(errors[0] / errors[1] > 1.9); assert!(errors[1] / errors[2] > 1.9);
}

#[test]
fn g1_successive_layer_reheats_first_and_replaces_exposed_area() {
    let law = curve(); let mut sim = simulation(&law);
    with_cx(|cx| {
        sim.apply_event(cx, &event(1, 0.0, vec![birth(1, 400.0, 10.0)]), 1e-9).unwrap();
        sim.advance_to(cx, 10.0, config(1.0)).unwrap();
        let cooled = temperature(&sim, 1); assert!(cooled < 400.0);
        let old = sim.cells()[&CellId(1)].phase;
        let mut layer = event(2, 10.0, vec![birth(2, 600.0, 0.0)]);
        layer.contacts.push(contact(1, 1, 2, 100.0));
        // Disable ambient here to isolate inter-layer reheating analytically.
        layer.environments.push((CellId(1), birth(1, 400.0, 0.0).environment));
        sim.apply_event(cx, &layer, 1e-9).unwrap();
        assert_eq!(old, sim.cells()[&CellId(1)].phase);
        close(sim.exposed_areas_m2().unwrap()[&CellId(1)], 9e-4, 1e-15);
        close(sim.exposed_areas_m2().unwrap()[&CellId(2)], 9e-4, 1e-15);
        let report = sim.advance_to(cx, 11.0, config(1.0)).unwrap();
        let q = (600.0 - cooled) / 102.0; // Q=dt*dT/(R+dt/Ca+dt/Cb).
        close(temperature(&sim, 1), cooled + q, 1e-8);
        close(temperature(&sim, 2), 600.0 - q, 1e-8);
        close(report.stored_energy_change_j, 0.0, 1e-8);
        close(sim.contacts()[&InterfaceId(1)].heat_a_to_b_j, -q, 1e-8);
    });
}

#[test]
fn g0_opening_contact_restores_surfaces_without_erasing_history() {
    let law = curve(); let mut sim = simulation(&law);
    with_cx(|cx| {
        let mut e = event(1, 0.0, vec![birth(1, 600.0, 0.0), birth(2, 300.0, 0.0)]);
        e.contacts.push(contact(1, 1, 2, 100.0)); sim.apply_event(cx, &e, 1e-9).unwrap();
        sim.advance_to(cx, 1.0, config(1.0)).unwrap();
        let heat = sim.contacts()[&InterfaceId(1)].heat_a_to_b_j;
        let mut open = event(2, 1.0, vec![]); open.contacts.push(ContactChange::Open(InterfaceId(1)));
        sim.apply_event(cx, &open, 1e-9).unwrap();
        close(sim.exposed_areas_m2().unwrap()[&CellId(1)], 1e-3, 1e-15);
        let t = temperature(&sim, 1); sim.advance_to(cx, 2.0, config(1.0)).unwrap();
        close(temperature(&sim, 1), t, 1e-8);
        assert_eq!(sim.contacts()[&InterfaceId(1)].heat_a_to_b_j, heat);
    });
}

#[test]
fn g1_source_and_boundary_heat_are_separate_from_material_birth() {
    let law = curve(); let mut sim = simulation(&law);
    with_cx(|cx| {
        let mut b = birth(1, 400.0, 10.0); b.environment.source_w = 2.0;
        sim.apply_event(cx, &event(1, 0.0, vec![b]), 1e-9).unwrap();
        let r = sim.advance_to(cx, 5.0, config(0.5)).unwrap();
        close(r.source_heat_j, 10.0, 1e-12); assert_eq!(r.incoming_enthalpy_j, 0.0);
        close(r.stored_energy_change_j, r.source_heat_j + r.boundary_heat_j, 1e-8);
    });
}

#[test]
fn g0_failed_event_does_not_partially_activate_material_or_advance_replay_cursor() {
    let law = curve(); let mut sim = simulation(&law);
    with_cx(|cx| {
        let mut invalid = event(1, 0.0, vec![birth(1, 400.0, 0.0)]);
        invalid.contacts.push(contact(1, 1, 999, 100.0));
        assert!(sim.apply_event(cx, &invalid, 1e-9).is_err()); assert!(sim.cells().is_empty());
        sim.apply_event(cx, &event(1, 0.0, vec![birth(1, 400.0, 0.0)]), 1e-9).unwrap();
    });
}

#[test]
fn g4_cancelled_or_work_limited_advance_is_atomic() {
    let law = curve();
    let mut sim = DepositionSimulation::new(vec![(MaterialId(7), DepositionMaterial {
        curve: &law, conductivity_w_m_k: 10.0, emissivity: 0.0 })],
        DepositionBudget { max_steps: 1, ..budget() }).unwrap();
    with_cx(|cx| { sim.apply_event(cx, &event(1, 0.0, vec![birth(1, 500.0, 10.0)]), 1e-9).unwrap(); });
    let phase = sim.cells()[&CellId(1)].phase;
    with_cancelled_cx(|cx| { assert!(sim.advance_to(cx, 1.0, config(1.0)).is_err()); });
    assert_eq!(sim.time_s(), 0.0); assert_eq!(sim.cells()[&CellId(1)].phase, phase);
    // Exhaust the next step after one real solve, without a long test loop.
    with_cx(|cx| { assert!(sim.advance_to(cx, 1.0, config(1e-8)).is_err()); });
    assert_eq!(sim.time_s(), 0.0); assert_eq!(sim.cells()[&CellId(1)].phase, phase);
}

#[test]
fn g0_rejects_overcoverage_invalid_frame_and_aggregate_biot() {
    let law = curve(); let mut sim = simulation(&law);
    with_cx(|cx| {
        let mut reflected = birth(1, 400.0, 0.0); reflected.geometry.material_axes[0][0] = -1.0;
        assert!(sim.apply_event(cx, &event(1, 0.0, vec![reflected]), 1e-9).is_err());
        let mut e = event(1, 0.0, vec![birth(1, 600.0, 0.0), birth(2, 300.0, 0.0), birth(3, 300.0, 0.0)]);
        // Each contact alone has Bi=0.06, but their sum at cell 1 exceeds 0.1.
        e.contacts = vec![contact(1, 1, 2, 1.0/0.06), contact(2, 1, 3, 1.0/0.06)];
        sim.apply_event(cx, &e, 1e-9).unwrap();
        assert!(sim.advance_to(cx, 1.0, config(1.0)).is_err()); assert_eq!(sim.time_s(), 0.0);
        let mut over = event(2, 0.0, vec![]);
        let ContactChange::Set(mut c) = contact(3, 2, 3, 100.0) else { unreachable!() };
        c.area_m2 = 1.0; over.contacts.push(ContactChange::Set(c));
        assert!(sim.apply_event(cx, &over, 1e-9).is_err()); assert_eq!(sim.contacts().len(), 2);
    });
}

#[test]
fn g1_equilibrium_latent_enthalpy_is_not_replaced_by_temperature_averaging() {
    let law = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x42; 32]), vec![
        knot(0.0, 200.0, 0.0), knot(200_000.0, 400.0, 0.0),
        knot(300_000.0, 400.0, 1.0), knot(700_000.0, 800.0, 1.0)]).unwrap();
    let mut sim = simulation(&law);
    with_cx(|cx| {
        let mut a = birth(1, 400.0, 0.0); a.specific_enthalpy_j_kg = 250_000.0;
        let mut b = birth(2, 600.0, 0.0); b.specific_enthalpy_j_kg = 500_000.0;
        let mut e = event(1, 0.0, vec![a,b]); e.contacts.push(contact(1, 1, 2, 100.0));
        sim.apply_event(cx, &e, 1e-9).unwrap(); let before = sim.stored_enthalpy_j();
        sim.advance_to(cx, 1.0, config(1.0)).unwrap();
        close(temperature(&sim, 1), 400.0, 1e-8);
        assert!(sim.cells()[&CellId(1)].phase.liquid_mass_fraction() > 0.5);
        close(sim.stored_enthalpy_j(), before, 1e-8);
    });
}
