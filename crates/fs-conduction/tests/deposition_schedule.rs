//! Shared resolved schedules, exact birth times, rollback, and native work bounds.
mod support;

use fs_blake3::ContentHash;
use fs_conduction::transient::deposition as resolved;
use fs_conduction::transient::enthalpy::deposition::implicit::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use support::with_cx;

type Simulation<'a> = DepositionSimulation<'a, resolved::CellGeometry>;

fn curve(offset: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([0x53; 32]), SolidLiquidPhase::Solid,
        [200.0, 800.0].into_iter().map(|t| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: offset + (t - 200.0) * 1000.0, temperature_k: t,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0,
        }).collect()).unwrap()
}
fn budget() -> DepositionBudget {
    DepositionBudget { max_cells: 1024, max_contacts: 2048, max_materials: 4,
        max_events: 4096, max_steps: 10000, max_operator_solves: 1_000_000 }
}
fn config(dt: f64) -> DepositionStepConfig {
    DepositionStepConfig { maximum_step_s: dt, solver_energy_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-7 }
}
fn simulation(curve: &EquilibriumEnthalpyPhaseCurve, emissivity: f64) -> Simulation<'_> {
    DepositionSimulation::new(vec![(MaterialId(7), DepositionMaterial {
        curve, conductivity_w_m_k: 10.0, emissivity })], budget()).unwrap()
}
fn boundary(htc: f64) -> resolved::Boundary {
    resolved::Boundary { ambient_k: 300.0, htc_w_m2_k: htc, reservoir: None, source_w: 0.0 }
}
fn cell(id: u64, t: f64, htc: f64) -> resolved::DepositedCell {
    resolved::DepositedCell { id: resolved::CellId(id), body_id: 5, segment_id: 100 + id,
        material_id: resolved::MaterialId(7), mass_kg: 0.001, specific_h_j_kg: (t-200.0)*1000.0,
        geometry: resolved::CellGeometry { start_m: [id as f64 * 0.001, 0.0, 0.0],
            end_m: [(id + 1) as f64 * 0.001, 0.0, 0.0], volume_m3: 1e-6, surface_area_m2: 1e-3 },
        boundary: boundary(htc) }
}
fn at(time_s: f64, event: resolved::Event) -> resolved::TimedEvent {
    resolved::TimedEvent { time_s, event }
}
fn born(time_s: f64, id: u64, t: f64, htc: f64) -> resolved::TimedEvent {
    at(time_s, resolved::Event::Birth(cell(id, t, htc)))
}
fn linked(time_s: f64, id: u64, a: u64, b: u64, area: f64) -> resolved::TimedEvent {
    at(time_s, resolved::Event::Contact(resolved::ThermalLink {
        id, a: resolved::CellId(a), b: resolved::CellId(b), area_m2: area, conductance_w_k: 0.01 }))
}
fn t(sim: &Simulation<'_>, id: u64) -> f64 { sim.cells()[&CellId(id)].phase.temperature_k() }
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a:e} != {b:e}"); }

#[test]
fn exact_horizon_birth_has_no_prebirth_cooling_and_future_cells_are_absent() {
    let law = curve(0.0); let mut sim = simulation(&law, 0.0);
    let events = vec![born(0.37, 1, 500.0, 10.0), born(0.91, 2, 600.0, 10.0)];
    with_cx(|cx| {
        let first = sim.run_resolved_until(cx, &events, EventId(1), 0.37, 300.0, config(1.0)).unwrap();
        assert_eq!(first.events, 1); assert_eq!(first.steps, 0); assert_eq!(first.operator_solves, 0);
        assert_eq!(sim.cells().len(), 1); assert_eq!(t(&sim, 1), 500.0);
        let next = sim.run_resolved_until(cx, &events[1..], EventId(2), 1.17, 300.0, config(1.0)).unwrap();
        assert_eq!(next.events, 1); assert_eq!(next.steps, 2); assert_eq!(next.operator_solves, 3);
        close(t(&sim, 1), 300.0 + 200.0 / 1.0054 / 1.0026, 1e-8);
        close(t(&sim, 2), 300.0 + 300.0 / 1.0026, 1e-8);
        assert_eq!(sim.cells()[&CellId(2)].birth_time_s, 0.91);
        assert_eq!(sim.cells()[&CellId(2)].birth.segment, SegmentId(102));
        assert_eq!(sim.cells()[&CellId(2)].birth.geometry.detail, cell(2, 600.0, 10.0).geometry);
        sim.run_resolved_until(cx, &[at(1.17, resolved::Event::Boundary(resolved::CellId(2), boundary(25.0)))],
            EventId(3), 1.17, 300.0, config(1.0)).unwrap();
        assert_eq!(sim.cells()[&CellId(2)].birth.environment.convection_w_m2_k, 10.0);
        assert_eq!(sim.cells()[&CellId(2)].environment.convection_w_m2_k, 25.0);
    });
}

#[test]
fn later_failure_rolls_back_earlier_births_thermal_aging_and_event_cursor() {
    let law = curve(0.0); let mut sim = simulation(&law, 0.0);
    let events = vec![born(0.0, 1, 500.0, 10.0), linked(1.0, 1, 1, 999, 1e-4)];
    with_cx(|cx| {
        assert!(sim.run_resolved_until(cx, &events, EventId(1), 2.0, 300.0, config(0.25)).is_err());
        assert_eq!(sim.time_s(), 0.0); assert!(sim.cells().is_empty()); assert!(sim.contacts().is_empty());
        sim.run_resolved_until(cx, &events[..1], EventId(1), 0.0, 300.0, config(0.25)).unwrap();
        assert_eq!(t(&sim, 1), 500.0);
        let before = sim.cells()[&CellId(1)].phase;
        let mut unsupported = cell(2, 500.0, 10.0);
        unsupported.boundary.reservoir = Some(resolved::Reservoir {
            area_m2: 1e-4, conductance_w_k: 0.01, temperature_k: 330.0 });
        assert!(sim.run_resolved_until(cx, &[at(1.0, resolved::Event::Birth(unsupported))],
            EventId(2), 2.0, 300.0, config(0.25)).is_err());
        assert_eq!(sim.cells()[&CellId(1)].phase, before); assert_eq!(sim.time_s(), 0.0);
    });
}

#[test]
fn checkpoint_restart_on_existing_step_boundaries_is_bit_identical() {
    let law = curve(0.0);
    let events = vec![born(0.0, 1, 400.0, 10.0), born(0.5, 2, 600.0, 10.0),
        linked(0.5, 1, 1, 2, 1e-4), at(1.5, resolved::Event::Separation(1))];
    let mut full = simulation(&law, 0.0); let mut resumed = full.clone();
    with_cx(|cx| {
        full.run_resolved_until(cx, &events, EventId(1), 2.0, 300.0, config(0.25)).unwrap();
        let first = resumed.run_resolved_until(cx, &events, EventId(1), 1.0, 300.0, config(0.25)).unwrap();
        assert_eq!(first.events, 3);
        let checkpoint = resumed.clone();
        resumed.run_resolved_until(cx, &events[3..], EventId(4), 2.0, 300.0, config(0.25)).unwrap();
        for id in [CellId(1), CellId(2)] {
            assert_eq!(full.cells()[&id].phase, resumed.cells()[&id].phase);
            assert_eq!(full.cells()[&id].birth.geometry, resumed.cells()[&id].birth.geometry);
        }
        assert_eq!(full.contacts()[&InterfaceId(1)].heat_a_to_b_j,
            resumed.contacts()[&InterfaceId(1)].heat_a_to_b_j);
        assert_eq!(checkpoint.time_s(), 1.0);
    });
}

#[test]
fn heterogeneous_material_enthalpy_datums_are_not_mixed() {
    let a = curve(0.0); let b = curve(1_000_000.0);
    let mut sim: Simulation<'_> = DepositionSimulation::new(vec![
        (MaterialId(7), DepositionMaterial { curve: &a, conductivity_w_m_k: 10.0, emissivity: 0.0 }),
        (MaterialId(99), DepositionMaterial { curve: &b, conductivity_w_m_k: 10.0, emissivity: 0.0 }),
    ], budget()).unwrap();
    let mut second = cell(2, 600.0, 0.0); second.material_id = resolved::MaterialId(99);
    second.specific_h_j_kg += 1_000_000.0;
    with_cx(|cx| {
        let r = sim.run_resolved_until(cx, &[born(0.0, 1, 300.0, 0.0),
            at(0.0, resolved::Event::Birth(second)), linked(0.0, 1, 1, 2, 1e-4)],
            EventId(1), 1.0, 300.0, config(1.0)).unwrap();
        close(t(&sim, 1), 300.0 + 300.0/102.0, 1e-8);
        close(t(&sim, 2), 600.0 - 300.0/102.0, 1e-8);
        close(r.incoming_enthalpy_j, 1500.0, 1e-8); close(sim.stored_enthalpy_j(), 1500.0, 1e-8);
    });
}

#[test]
fn nonlinear_radiation_closes_energy_and_buried_surfaces_do_not_radiate() {
    let law = curve(0.0); let mut sim = simulation(&law, 0.8);
    with_cx(|cx| {
        let r = sim.run_resolved_until(cx, &[born(0.0, 1, 500.0, 0.0)],
            EventId(1), 1.0, 300.0, config(1.0)).unwrap();
        let flux = 0.8 * fs_conduction::STEFAN_BOLTZMANN_W_M2_K4 * 1e-3
            * (300.0_f64.powi(4) - t(&sim, 1).powi(4));
        close(r.boundary_heat_j, flux, 1e-8); close(t(&sim, 1) - 500.0, flux, 1e-8);
        let mut buried = simulation(&law, 0.8);
        let mut powered = cell(1, 500.0, 10.0); powered.boundary.source_w = 2.0;
        let r = buried.run_resolved_until(cx, &[at(0.0, resolved::Event::Birth(powered)),
            born(0.0, 2, 500.0, 10.0), linked(0.0, 1, 1, 2, 1e-3)],
            EventId(1), 1.0, 300.0, config(1.0)).unwrap();
        assert_eq!(r.boundary_heat_j, 0.0); close(r.source_heat_j, 2.0, 1e-12);
        close(r.stored_energy_change_j - r.incoming_enthalpy_j, 2.0, 1e-8);
        assert_eq!(buried.exposed_areas_m2().unwrap()[&CellId(1)], 0.0);
    });
}

#[test]
fn shared_schedule_explicit_and_implicit_lanes_converge_to_the_same_cooling_law() {
    let law = curve(0.0); let events = vec![born(0.0, 1, 500.0, 10.0)];
    let exact = 300.0 + 200.0 * (-0.2_f64).exp();
    let mut gaps = Vec::new();
    with_cx(|cx| {
        for dt in [4.0, 2.0, 1.0] {
            let explicit = resolved::DepositionSimulation::new(events.clone(),
                &[(resolved::MaterialId(7), &law as &dyn resolved::ThermalLaw)],
                resolved::DepositionBudget { max_cells: 16, max_links: 16, max_events: 16, max_substeps: 1000 }).unwrap();
            let (explicit, _) = explicit.advance_to(20.0, resolved::StepConfig {
                max_dt_s: dt, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap();
            let mut implicit = simulation(&law, 0.0);
            implicit.run_resolved_until(cx, &events, EventId(1), 20.0, 300.0, config(dt)).unwrap();
            let lo = explicit.temperature_k(resolved::CellId(1)).unwrap(); let hi = t(&implicit, 1);
            assert!(lo < exact && exact < hi); gaps.push(hi-lo);
        }
    });
    assert!(gaps[0]/gaps[1] > 1.9 && gaps[1]/gaps[2] > 1.9);
}

#[test]
fn native_operator_work_is_linear_in_active_cells_and_links_per_step() {
    let law = curve(0.0);
    with_cx(|cx| {
        for n in [16_u64, 64] {
            let mut sim = simulation(&law, 0.0);
            let mut events: Vec<_> = (1..=n).map(|i| born(0.0, i, 500.0, 10.0)).collect();
            events.extend((1..n).map(|i| linked(0.0, i, i, i+1, 1e-4)));
            let r = sim.run_resolved_until(cx, &events, EventId(1), 10.0, 300.0, config(1.0)).unwrap();
            assert_eq!(r.steps, 10); assert_eq!(r.operator_solves, 10 * (2*n as usize - 1));
            assert_eq!(sim.cells().len(), n as usize); assert_eq!(sim.contacts().len(), n as usize-1);
            assert!(r.energy_residual_j.abs() <= 1e-7);
        }
    });
}

#[test]
fn isothermal_partition_preserves_mass_energy_and_external_cooling_rate() {
    // A conservative control-volume partition, not a spatial PDE convergence claim.
    let law = curve(0.0); let mut totals = Vec::new();
    with_cx(|cx| {
        for count in [1_u64, 2, 4] {
            let mut sim = simulation(&law, 0.0);
            let events: Vec<_> = (1..=count).map(|i| {
                let mut c = cell(i, 500.0, 10.0);
                c.mass_kg /= count as f64; c.geometry.volume_m3 /= count as f64;
                c.geometry.surface_area_m2 /= count as f64;
                at(0.0, resolved::Event::Birth(c))
            }).collect();
            let r = sim.run_resolved_until(cx, &events, EventId(1), 10.0, 300.0, config(1.0)).unwrap();
            close(sim.total_mass_kg(), 0.001, 1e-15);
            close(r.incoming_enthalpy_j, 300.0, 1e-9);
            totals.push(sim.stored_enthalpy_j());
        }
    });
    close(totals[0], totals[1], 1e-8); close(totals[0], totals[2], 1e-8);
}
