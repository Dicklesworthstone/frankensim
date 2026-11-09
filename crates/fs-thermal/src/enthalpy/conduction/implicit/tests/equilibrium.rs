//! G1/G3/G4: analytical stationary balances, original-source response and refusal.
use super::super::coupled_source::{ImplicitCoupledHeatSourceLaw, ImplicitCoupledHeatSourceState};
use super::super::equilibrium::{ThermalOperatingControls, ThermalOperatingError};
use super::super::harmonic::ThermalImpedanceBudget;
use super::*;

fn operating_controls() -> ThermalOperatingControls {
    ThermalOperatingControls {
        maximum_newton_iterations: 32,
        maximum_backtracks: 24,
        maximum_residual_evaluations: 512,
        cell_tolerance_w: 1e-10,
        total_tolerance_w: 1e-10,
    }
}

#[test]
fn stationary_chain_and_graph_close_each_heat_row_and_feed_native_impedance() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (mut thermal, conduction) = setup(
        &law,
        &[500.0, 200.0],
        Some(300.0),
        vec![TimedEvent {
            time_s: 2.0,
            event: Event::Boundary(CellId(0), boundary(Some(600.0))),
        }],
    );
    thermal.cells.get_mut(&CellId(1)).unwrap().boundary.source_w = 25.0;
    let before = thermal.cells().clone();
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let point = chain
        .solve_operating_point(&thermal, operating_controls(), &mut workspace, || Ok(()))
        .unwrap();
    let sparse = graph
        .solve_operating_point(&thermal, operating_controls(), &mut workspace, || Ok(()))
        .unwrap();
    // G_link=G_bath=1 W/K: both the 25 W internal flux and the bath debit
    // must balance their individual cells, not merely cancel in the sum.
    for (i, expected) in [325.0, 350.0].into_iter().enumerate() {
        near(point.temperatures_k()[i], expected, 1e-11);
        near(sparse.temperatures_k()[i], expected, 1e-11);
        near(point.specific_enthalpies_j_kg()[i], expected, 1e-11);
        assert!(point.net_powers_w()[i].abs() <= 1e-10);
    }
    assert_eq!(point.report().nonlinear_iterations, 1);
    assert!(sparse.report().linear_iterations > 0);
    let model = point
        .prepare_thermal_impedance(
            ThermalImpedanceBudget {
                maximum_cells: 2,
                maximum_dense_bytes: 4096,
                maximum_dense_work: 1000,
                equilibrium_tolerance_w: 1e-10,
                relative_residual_tolerance: 1e-12,
            },
            || Ok(()),
        )
        .unwrap();
    let response = model
        .factor(0.0, || Ok(()))
        .unwrap()
        .impedance_column(CellId(1), || Ok(()))
        .unwrap();
    near(response[0].re, 1.0, 1e-12);
    near(response[1].re, 2.0, 1e-12);
    assert_eq!(point.operating_time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.cursor, 3); // Future bath event was not applied.
}

#[test]
fn stationary_nonlinear_conductivity_uses_original_material_tangent() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let (mut thermal, conduction) = setup(&law, &[150.0], Some(300.0), Vec::new());
    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 100.0;
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    // At the minimum of the quadratic heat loss the exact Jacobian vanishes;
    // this is a real singular Newton seed, not permission to perturb the model.
    assert!(
        chain
            .solve_operating_point(&thermal, operating_controls(), &mut workspace, || Ok(()))
            .is_err()
    );
    thermal.cells.get_mut(&CellId(0)).unwrap().specific_h_j_kg = 180.0;
    let point = chain
        .solve_operating_point(&thermal, operating_controls(), &mut workspace, || Ok(()))
        .unwrap();
    let exact = (300.0 + 130_000.0_f64.sqrt()) / 2.0;
    near(point.temperatures_k()[0], exact, 1e-10);
    assert!(point.report().nonlinear_iterations > 1);
    assert!(point.report().maximum_cell_residual_w <= 1e-10);
    near(
        point.temperatures_k()[0] * (point.temperatures_k()[0] - 300.0) / 100.0,
        100.0,
        1e-10,
    );
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 180.0);
}

struct CoupledSource<'a> {
    calorics: [&'a dyn ThermalLaw; 2],
}
impl ImplicitCoupledHeatSourceLaw for CoupledSource<'_> {
    fn caloric_laws(&self) -> &[&dyn ThermalLaw] {
        &self.calorics
    }
    fn power_dependencies(&self) -> &[(usize, usize)] {
        &[(0, 1), (1, 0)]
    }
    fn state_at_specific_enthalpies_and_time(
        &self,
        h: &[f64],
        _: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitCoupledHeatSourceState, DepositionError> {
        checkpoint()?;
        Ok(ImplicitCoupledHeatSourceState {
            powers_w: vec![10.0 + 0.1 * (h[1] - 300.0), 5.0 + 0.2 * (h[0] - 300.0)],
            power_enthalpy_derivatives_w_kg_j: vec![0.1, 0.2],
        })
    }
}

#[test]
fn stationary_graph_keeps_nonlocal_heat_feedback_in_the_same_sparse_solve() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[280.0, 480.0, 200.0], Some(300.0), Vec::new());
    let source = CoupledSource {
        calorics: [&law, &law],
    };
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_coupled_heat_source(&thermal, &[CellId(0), CellId(2)], &source, 2, || Ok(()))
    .unwrap();
    let point = graph
        .solve_operating_point(
            &thermal,
            operating_controls(),
            &mut ImplicitConductionWorkspace::new(3).unwrap(),
            || Ok(()),
        )
        .unwrap();
    // Eliminate the middle node and both native source laws analytically.
    let x = 16.0 / 0.66;
    for (actual, expected) in
        point
            .temperatures_k()
            .iter()
            .zip([300.0 + x, 305.0 + 1.2 * x, 310.0 + 1.4 * x])
    {
        near(*actual, expected, 1e-10);
    }
    assert_eq!(point.report().nonlinear_iterations, 1);
    assert!(point.report().maximum_cell_residual_w <= 1e-10);
}

#[test]
fn stationary_solve_refuses_unconstrained_isolated_inventory_despite_zero_total_power() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[400.0, 200.0], None, Vec::new());
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        graph.solve_operating_point(
            &thermal,
            operating_controls(),
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        ),
        Err(ThermalOperatingError::Unanchored { .. })
    ));
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, 200.0);
}

#[test]
fn stationary_latent_plateau_is_preserved_and_never_regularized() {
    let law = Law {
        nonlinear_k: false,
        plateau: true,
    };
    let (mut thermal, conduction) = setup(&law, &[410.0], Some(300.0), Vec::new());
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let point = chain
        .solve_operating_point(&thermal, operating_controls(), &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(point.specific_enthalpies_j_kg(), &[410.0]);
    assert_eq!(point.temperatures_k(), &[300.0]);
    assert_eq!(point.report().nonlinear_iterations, 0);
    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 1.0;
    assert!(
        chain
            .solve_operating_point(&thermal, operating_controls(), &mut workspace, || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 410.0);
}

#[test]
fn stationary_domain_work_and_cancellation_refuse_without_changing_history() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (mut thermal, conduction) = setup(&law, &[400.0], Some(300.0), Vec::new());
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut budget = operating_controls();
    budget.maximum_residual_evaluations = 1;
    assert!(matches!(
        chain.solve_operating_point(&thermal, budget, &mut workspace, || Ok(())),
        Err(ThermalOperatingError::Original(
            ImplicitConductionError::Thermal(DepositionError::Budget(_))
        ))
    ));
    let mut calls = 0;
    assert!(matches!(
        chain.solve_operating_point(&thermal, operating_controls(), &mut workspace, || {
            calls += 1;
            if calls >= 25 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ThermalOperatingError::Original(
            ImplicitConductionError::Thermal(DepositionError::Cancelled)
        ))
    ));
    // The only mathematical root is beyond the original material domain.
    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 900.0;
    budget = operating_controls();
    budget.maximum_newton_iterations = 4;
    assert!(matches!(
        chain.solve_operating_point(&thermal, budget, &mut workspace, || Ok(())),
        Err(ThermalOperatingError::NotConverged { .. })
    ));
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
    assert_eq!(thermal.time_s(), 0.0);
}
