//! Native hydraulic -> original thermal inventory -> independently checked heat.
use super::*;
use crate::graph::{GraphSolution, LossGraph};
use crate::graph::tests::{boundary, config, edge, with_cx};
use fs_material::fluid_transport::FluidTransportState;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, ThermalLaw, TimedEvent,
};

struct Law { datum: f64, plateau: bool }
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let h = h - self.datum;
        if !h.is_finite() || !(-100.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + if self.plateau {
            if h < 20.0 { h } else if h < 80.0 { 20.0 } else { h - 60.0 }
        } else { h })
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
}

fn thermal<'a>(law: &'a Law, masses: [f64; 2], values: [f64; 2], heating: bool)
    -> DepositionSimulation<'a> {
    let events = (0..2).map(|i| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(i as u64), segment_id: i as u64, body_id: 1, material_id: MaterialId(1),
            mass_kg: masses[i], specific_h_j_kg: law.datum + values[i],
            geometry: CellGeometry { start_m: [i as f64, 0.0, 0.0], end_m: [i as f64 + 1.0, 0.0, 0.0],
                volume_m3: masses[i], surface_area_m2: 1.0 },
            boundary: Boundary { ambient_k: 400.0, htc_w_m2_k: if heating { 1.0 } else { 0.0 },
                reservoir: None, source_w: 0.0 },
        }),
    }).collect();
    DepositionSimulation::new(events, &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 8, max_events: 4, max_substeps: 10_000 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-8 }, || Ok(())).unwrap().0
}

fn loop_flow(cx: &Cx<'_>, reverse_declaration: bool) -> GraphSolution {
    let state = |rho| FluidTransportState::try_new(rho, 1.0, 0.1, 0.001).unwrap();
    let return_leg = if reverse_declaration { edge("return", 0, 1, 1.0) }
        else { edge("return", 1, 0, 1.0) };
    let mut controls = config();
    controls.absolute_flow_tolerance = fs_qty::VolumetricFlowRate::new(1e-15);
    controls.relative_flow_tolerance = 1e-15;
    LossGraph::new(2, vec![edge("riser", 0, 1, 1.0), return_leg]).unwrap()
        .with_boussinesq_gravity(&state(1.0), &[state(0.75), state(1.0)],
            &[[0.0; 3], [0.0, 0.0, 1.0]], [0.0, 0.0, -8.0], 0.3, cx).unwrap()
        .solve(&[boundary(0, 0.0)], controls, cx).unwrap()
}

fn open_flow(cx: &Cx<'_>, reverse_declaration: bool) -> GraphSolution {
    let edge = if reverse_declaration { edge("duct", 1, 0, 1.0) } else { edge("duct", 0, 1, 1.0) };
    LossGraph::new(2, vec![edge]).unwrap()
        .solve(&[boundary(0, 1.0), boundary(1, 0.0)], config(), cx).unwrap()
}

fn streams(datum: f64) -> Vec<MassBoundary> {
    vec![
        MassBoundary::Inlet { id: 20, cell: CellId(0), kg_s: 1.0, material_id: MaterialId(1),
            specific_h_j_kg: datum + 120.0 },
        MassBoundary::Outlet { id: 21, cell: CellId(1), kg_s: 1.0 },
    ]
}

fn binding(cx: &Cx<'_>, flow: &GraphSolution, state: &DepositionSimulation<'_>, boundaries: Vec<MassBoundary>)
    -> GraphEnthalpyTransport {
    GraphEnthalpyTransport::new(cx, flow, Density::new(1.0), &[CellId(0), CellId(1)], state, boundaries, 8).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.16e} != {expected:.16e}");
}

#[test]
fn g1_implicit_closed_loop_conserves_unequal_inventory_and_large_caloric_datum() {
    with_cx(|cx| {
        for datum in [0.0, 1_099_511_627_776.0] {
            for reverse in [false, true] {
                let law = Law { datum, plateau: false };
                let old = thermal(&law, [1.0, 3.0], [0.0, 44.0], false);
                let flow = loop_flow(cx, reverse);
                assert_eq!(flow.branches[0].flow.value(), 1.0);
                let transport = binding(cx, &flow, &old, vec![]);
                // Independent two-volume backward-Euler balance at dt=2:
                // 3 h0 - 2 h1 = 0; -2 h0 + 5 h1 = 132.
                // Courant=2 exceeds an explicit monotonic step's allowance.
                let (next, audit) = transport.advance_advection_implicit(cx, &old, 2.0, 1e-8, 1e-8).unwrap();
                near(next.cells()[&CellId(0)].specific_h_j_kg - datum, 24.0, 1e-8);
                near(next.cells()[&CellId(1)].specific_h_j_kg - datum, 36.0, 1e-8);
                near(audit.stored_energy_change_j, 0.0, 1e-8);
                assert_eq!(audit.external_heat_j, 0.0);
                assert_eq!(audit.incoming_mass_kg, 0.0);
                assert_eq!(old.time_s(), 0.0);
                assert_eq!(next.time_s(), 2.0);
            }
        }
    });
}

#[test]
fn g1_implicit_supply_uses_upwind_enthalpy_and_latent_storage_without_temperature_substitution() {
    with_cx(|cx| {
        for plateau in [false, true] {
            let law = Law { datum: 0.0, plateau };
            let old = thermal(&law, [1.0; 2], [0.0; 2], false);
            for reverse in [false, true] {
                let flow = open_flow(cx, reverse);
                let transport = binding(cx, &flow, &old, streams(0.0));
                let (next, audit) = transport.advance_advection_implicit(cx, &old, 1.0, 1e-8, 1e-8).unwrap();
                near(next.cells()[&CellId(0)].specific_h_j_kg, 60.0, 1e-8);
                near(next.cells()[&CellId(1)].specific_h_j_kg, 30.0, 1e-8);
                near(audit.external_heat_j, 90.0, 1e-8);
                near(audit.stored_energy_change_j, 90.0, 1e-8);
                if plateau {
                    assert_eq!(next.temperature_k(CellId(0)).unwrap(), 320.0);
                    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 320.0);
                }
            }
        }
    });
}

#[test]
fn g1_native_explicit_storage_converges_to_two_tank_residence_response() {
    with_cx(|cx| {
        let law = Law { datum: 0.0, plateau: false };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let flow = open_flow(cx, true);
        let transport = binding(cx, &flow, &old, streams(0.0));
        let exact = [120.0*(1.0-(-0.5_f64).exp()), 120.0*(1.0-1.5*(-0.5_f64).exp())];
        let mut errors = Vec::new();
        for dt in [0.01, 0.005] {
            let (next, energy, advection) = transport.advance_to(cx, &old, 0.5,
                StepConfig { max_dt_s: dt, energy_tolerance_j: 1e-8 }).unwrap();
            let error = (0..2).map(|i| (next.cells()[&CellId(i as u64)].specific_h_j_kg-exact[i]).abs()).fold(0.0, f64::max);
            errors.push(error);
            near(energy.external_heat_j, advection.external_enthalpy_j, 1e-8);
            near(energy.stored_energy_change_j, energy.external_heat_j, 1e-8);
            assert!(advection.maximum_courant <= 0.8);
            assert_eq!(energy.incoming_mass_kg, 0.0);
        }
        assert!(errors[1] < 0.55*errors[0], "first-order time refinement: {errors:?}");
    });
}

#[test]
fn g3_original_wall_exchange_composes_with_advection_and_direct_implicit_refuses_it() {
    with_cx(|cx| {
        let law = Law { datum: 0.0, plateau: false };
        let old = thermal(&law, [1.0; 2], [0.0; 2], true);
        let flow = loop_flow(cx, false);
        let transport = binding(cx, &flow, &old, vec![]);
        let (next, energy, advection) = transport.advance_to(cx, &old, 0.1,
            StepConfig { max_dt_s: 0.01, energy_tolerance_j: 1e-8 }).unwrap();
        assert!(next.temperature_k(CellId(0)).unwrap() > 300.0);
        assert_eq!(advection.external_enthalpy_j, 0.0);
        assert!(energy.external_heat_j > 0.0);
        near(energy.stored_energy_change_j, energy.external_heat_j, 1e-8);
        assert!(matches!(transport.advance_advection_implicit(cx, &old, 1.0, 1e-8, 1e-8),
            Err(TransientTransportError::Invalid(_))));
    });
}

#[test]
fn g0_pressure_boundaries_do_not_create_fluid_or_heat_and_binding_must_match_storage() {
    with_cx(|cx| {
        let law = Law { datum: 0.0, plateau: false };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let flow = open_flow(cx, false);
        assert!(GraphEnthalpyTransport::new(cx, &flow, Density::new(1.0), &[CellId(0), CellId(1)],
            &old, vec![], 8).is_err(), "fixed pressure does not authorize missing mass streams");
        assert!(GraphEnthalpyTransport::new(cx, &flow, Density::new(2.0), &[CellId(0), CellId(1)],
            &old, streams(0.0), 8).is_err(), "density must match declared inventory and volume");
        assert!(GraphEnthalpyTransport::new(cx, &flow, Density::new(1.0), &[CellId(0), CellId(0)],
            &old, streams(0.0), 8).is_err(), "two hydraulic nodes cannot alias one stored inventory");
        let mut closed = loop_flow(cx, false);
        closed.branches[0].flow = fs_qty::VolumetricFlowRate::new(1.001);
        assert!(GraphEnthalpyTransport::new(cx, &closed, Density::new(1.0), &[CellId(0), CellId(1)],
            &old, vec![], 8).is_err(), "a leaking pressure datum is not a closed loop");
    });
}

#[test]
fn g4_cancelled_steps_preserve_the_original_inventory_and_time() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 17, kernel_id: 711, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let law = Law { datum: 0.0, plateau: false };
        let old = thermal(&law, [1.0; 2], [0.0, 40.0], false);
        let flow = loop_flow(&cx, false);
        let transport = binding(&cx, &flow, &old, vec![]);
        gate.request();
        assert!(transport.advance_advection_implicit(&cx, &old, 1.0, 1e-8, 1e-8).is_err());
        assert!(transport.advance_to(&cx, &old, 1.0,
            StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-8 }).is_err());
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 40.0);
    });
}
