//! Liquid face traces, stored-energy tangents and native graph balances.
use super::*;
use crate::enthalpy::conduction::{ConductionObserver, ConductivityLaw};
use crate::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionState, ImplicitConductionWorkspace,
};

fn plan(state: &DepositionSimulation<'_>) -> LimitedEnthalpyReconstruction {
    LimitedEnthalpyReconstruction::try_new(
        state,
        &[AxialCellStencil {
            previous: CellId(0),
            cell: CellId(1),
            next: CellId(2),
        }],
        1,
        || Ok(()),
    )
    .unwrap()
}

fn liquid_flow(law: &PhaseLaw, rate: f64) -> Advected<()> {
    Advected::new((), flows(3, rate), 3)
        .unwrap()
        .with_liquid_carrier(carrier(law))
        .unwrap()
}

#[test]
fn liquid_reconstruction_does_not_advect_stationary_solid_latent_inventory() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = PhaseLaw { datum, identity: 7, melting_range_k: 0.0 };
        let initial = prepare(&law, &[20.0, 50.0, 80.0], 2.0);
        let reconstruction = plan(&initial);
        for rate in [-2.0, 2.0] {
            let flow = liquid_flow(&law, rate);
            let (endpoint, tangent) = flow
                .reconstructed_endpoint(&initial, &reconstruction, true, || Ok(()))
                .unwrap();
            assert!(endpoint.transfers.iter().all(|transfer| transfer.power_w == 0.0));
            assert!(tangent.is_empty());
            let step = BackwardEulerStep::new(&initial, 4.0, || Ok(())).unwrap();
            let trial = step
                .evaluate_with_transport(
                    &[0.0; 3],
                    |state| {
                        flow.endpoint_transport_with_reconstruction(
                            state, &reconstruction, || Ok(()),
                        )
                    },
                    || Ok(()),
                )
                .unwrap();
            let (next, (), audit) = trial.accept(1e-10, 1e-10, || Ok(())).unwrap();
            assert_eq!(next.cells(), initial.cells());
            assert_eq!(audit.external_heat_j, 0.0);
            assert_eq!(audit.stored_energy_change_j, 0.0);
        }
    }
}

#[test]
fn explicit_liquid_reconstruction_keeps_stationary_latent_inventory_and_bulk_limit() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = PhaseLaw { datum, identity: 7, melting_range_k: 0.0 };
        for values in [[20.0, 50.0, 80.0], [120.0, 160.0, 210.0]] {
            let initial = prepare(&law, &values, 2.0);
            for rate in [-2.0, 2.0] {
                let flow = liquid_flow(&law, rate).with_reconstruction(plan(&initial)).unwrap();
                let (next, observed, audit) = initial.advance_observed(
                    &flow, 0.125, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-3 },
                    || Ok(()),
                ).unwrap();
                assert_eq!(audit.external_heat_j, 0.0);
                assert!(observed.audit().maximum_courant <= 0.4);
                if values[0] < 100.0 {
                    assert_eq!(next.cells(), initial.cells());
                    assert_eq!(audit.stored_energy_change_j, 0.0);
                } else {
                    let bulk = Advected::new((), flows(3, rate), 3).unwrap()
                        .with_reconstruction(plan(&initial)).unwrap();
                    let (expected, _, expected_audit) = initial.advance_observed(
                        &bulk, 0.125, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-3 },
                        || Ok(()),
                    ).unwrap();
                    assert_eq!(next.cells(), expected.cells());
                    assert_eq!(audit, expected_audit);
                }
            }
        }
    }
}

#[test]
fn liquid_reconstruction_recovers_affine_bulk_traces_on_unequal_cells_and_reversal() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = PhaseLaw { datum, identity: 7, melting_range_k: 0.0 };
        // h_liquid(x)=110+20*x. The three centers are x=0.5, 2.5, 5.
        let mut state = prepare(&law, &[120.0, 160.0, 210.0], 2.0);
        for (i, [start, end]) in [[0.0, 1.0], [1.0, 4.0], [4.0, 6.0]]
            .into_iter()
            .enumerate()
        {
            let cell = state.cells.get_mut(&CellId(i as u64)).unwrap();
            cell.geometry.start_m = [start, 0.0, 0.0];
            cell.geometry.end_m = [end, 0.0, 0.0];
        }
        let reconstruction = plan(&state);
        for rate in [-2.0, 2.0] {
            let bulk = Advected::new((), flows(3, rate), 3).unwrap();
            let liquid = bulk.clone().with_liquid_carrier(carrier(&law)).unwrap();
            let (bulk_endpoint, bulk_tangent) = bulk
                .reconstructed_endpoint(&state, &reconstruction, true, || Ok(()))
                .unwrap();
            let (liquid_endpoint, liquid_tangent) = liquid
                .reconstructed_endpoint(&state, &reconstruction, true, || Ok(()))
                .unwrap();
            assert_eq!(liquid_endpoint.transfers, bulk_endpoint.transfers);
            assert_eq!(liquid_tangent, bulk_tangent);
            let face = liquid_endpoint.transfers.iter()
                .find(|face| face.a == CellId(1)).unwrap();
            // Native internal powers subtract the same component reference h0=120.
            close(
                face.power_w / rate.abs() + 120.0,
                if rate > 0.0 { 190.0 } else { 130.0 },
                1e-12,
            );
        }
    }
}

fn heat_rates(
    state: &DepositionSimulation<'_>,
    flow: &Advected<()>,
    reconstruction: &LimitedEnthalpyReconstruction,
) -> [f64; 3] {
    let endpoint = flow
        .endpoint_transport_with_reconstruction(state, reconstruction, || Ok(()))
        .unwrap();
    assert!(endpoint.heat.bulk.sources_w.is_empty());
    let mut rates = [0.0; 3];
    for transfer in endpoint.transfers {
        rates[transfer.a.0 as usize] -= transfer.power_w;
        rates[transfer.b.0 as usize] += transfer.power_w;
    }
    rates
}

#[test]
fn liquid_reconstruction_tangent_uses_each_samples_actual_phase() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    for values in [
        [20.0, 160.0, 260.0], // Left active neighbor is in coexistence.
        [260.0, 160.0, 20.0], // Right active neighbor is in coexistence.
        [130.0, 170.0, 250.0],
        [20.0, 50.0, 80.0],
    ] {
        let state = prepare(&law, &values, 2.0);
        let reconstruction = plan(&state);
        for rate in [-0.5, 0.5] {
            let flow = liquid_flow(&law, rate);
            let (_, corrections) = flow
                .reconstructed_endpoint(&state, &reconstruction, true, || Ok(()))
                .unwrap();
            let mut jacobian = [[0.0; 3]; 3];
            for (row, col, value) in flow.implicit_enthalpy_tangent(|| Ok(())).unwrap() {
                jacobian[row.0 as usize][col.0 as usize] += value
                    * flow.implicit_carried_enthalpy_slope(&state.cells()[&col]).unwrap();
            }
            for (row, col, value) in corrections {
                jacobian[row.0 as usize][col.0 as usize] += value;
            }
            // Independent central differences of the complete physical heat
            // rates, away from liquidus and limiter switches.
            let epsilon = 1e-4;
            for col in 0..3 {
                let mut plus = state.clone();
                let mut minus = state.clone();
                plus.cells.get_mut(&CellId(col as u64)).unwrap().specific_h_j_kg += epsilon;
                minus.cells.get_mut(&CellId(col as u64)).unwrap().specific_h_j_kg -= epsilon;
                let plus = heat_rates(&plus, &flow, &reconstruction);
                let minus = heat_rates(&minus, &flow, &reconstruction);
                for row in 0..3 {
                    close(jacobian[row][col], (plus[row] - minus[row]) / (2.0 * epsilon), 1e-8);
                }
                close(jacobian.iter().map(|row| row[col]).sum(), 0.0, 1e-14);
            }
        }
    }
}

#[test]
fn liquid_reconstruction_falls_back_beside_dry_samples_without_admitting_solid_flow() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let state = prepare(&law, &[-20.0, 160.0, 240.0], 2.0);
    let reconstruction = plan(&state);
    for rate in [-0.5, 0.5] {
        // Cell 0 is a geometric sample with zero flow; 1 and 2 recirculate.
        let flow = Advected::new((), vec![
            MassFlow { id: 0, a: CellId(0), b: CellId(1), kg_s: 0.0 },
            MassFlow { id: 1, a: CellId(1), b: CellId(2), kg_s: rate },
            MassFlow { id: 2, a: CellId(2), b: CellId(1), kg_s: rate },
        ], 3).unwrap().with_liquid_carrier(carrier(&law)).unwrap();
        let upwind = flow.endpoint_transport(&state, || Ok(())).unwrap();
        let (limited, tangent) = flow
            .reconstructed_endpoint(&state, &reconstruction, true, || Ok(()))
            .unwrap();
        assert_eq!(limited.transfers, upwind.transfers);
        assert!(tangent.is_empty());
        assert!(matches!(
            liquid_flow(&law, rate).endpoint_transport_with_reconstruction(
                &state, &reconstruction, || Ok(()),
            ),
            Err(DepositionError::InvalidInput("nonzero liquid mass flow crosses a fully solid cell"))
        ));
    }
}

struct ImplicitPhase<'a>(&'a PhaseLaw);
impl ConductivityLaw for ImplicitPhase<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.0.temperature_k(h)?;
        Ok(1.0)
    }
}
impl ImplicitConductionLaw for ImplicitPhase<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw { self.0 }
    fn conductivity_law(&self) -> &dyn ConductivityLaw { self }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [self.0.datum - 100.0, self.0.datum + 300.0]
    }
    fn state_at_specific_enthalpy(&self, h: f64) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.0.temperature_k(h)?,
            temperature_slope_kg_k_j: if (0.0..=100.0).contains(&(h - self.0.datum)) {
                0.0
            } else {
                0.5
            },
            conductivity_w_m_k: 1.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

#[test]
fn liquid_reconstruction_native_graph_crosses_liquidus_with_both_dynamic_tangents() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let initial = prepare(&law, &[280.0, 20.0, 10.0], 2.0);
    let reconstruction = plan(&initial);
    let flow = liquid_flow(&law, 1.0);
    let implicit = ImplicitPhase(&law);
    let conduction = ConductionObserver::new(
        &initial, &[(MAT, &implicit)], &[], &[], || Ok(()),
    ).unwrap();
    let solver = ImplicitConductionGraph::new(
        &initial, &conduction, &[(MAT, &implicit)], Default::default(), || Ok(()),
    ).unwrap().with_implicit_advection_reconstruction(
        &initial, &flow, &reconstruction, || Ok(()),
    ).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    let (next, report) = solver.advance_to(
        &initial, 2.0,
        ImplicitConductionControls {
            maximum_dt_s: 2.0,
            maximum_newton_iterations: 8,
            maximum_halvings: 0,
            cell_tolerance_j: 1e-10,
            total_tolerance_j: 1e-9,
            ..ImplicitConductionControls::default()
        },
        &mut workspace, || Ok(()),
    ).unwrap();
    // Independent BE balances at lambda=q*dt/m=2:
    // h0=280+2*(100-h0), h1=20+2*(h0-(h1+100)/2),
    // h2=10+2*((h1+100)/2-100).
    // The endpoint's right minmod branch is consistent: 20 < 40.
    for (i, expected) in [160.0, 120.0, 30.0].into_iter().enumerate() {
        close(next.cells()[&CellId(i as u64)].specific_h_j_kg, expected, 1e-9);
    }
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    assert!(report.maximum_cell_residual_j <= 1e-10);
    close(report.energy.external_heat_j, 0.0, 0.0);
    close(report.energy.stored_energy_change_j, 0.0, 1e-9);
    close(next.energy_j(), initial.energy_j(), 1e-9);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(flow.audit(), Default::default());
}
