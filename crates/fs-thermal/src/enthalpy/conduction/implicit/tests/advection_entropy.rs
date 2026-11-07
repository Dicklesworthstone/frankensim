//! Entropy on the original implicit upwind/conduction trajectory.
use super::*;

fn accounted(audit: EntropyAudit) -> f64 {
    let flow = audit.advection.unwrap();
    audit.reservoir_entropy_input_j_k
        + audit.conduction_generation_j_k
        + audit.unclassified_source_at_cell_j_k
        + flow.boundary_entropy_input_j_k
        + flow.upwind_mixing_j_k
        + flow.balance_residual_j_k
        + audit.temporal_defect_j_k
}

#[test]
fn implicit_entropy_closed_loop_separates_spatial_mixing_and_temporal_defect() {
    for datum in [0.0, 1e9] {
        let law = ShiftedLaw {
            base: Law { nonlinear_k: false, plateau: false },
            datum,
        };
        let (thermal, conduction) = setup(&law, &[400.0 + datum, 200.0 + datum], None, vec![]);
        let flow = circulating((), 30.5);
        let chain = ImplicitConductionChain::new(
            &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
        ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap()
            .with_implicit_advection_entropy(&thermal, || Ok(())).unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal, &conduction, &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(), || Ok(()),
        ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap()
            .with_implicit_advection_entropy(&thermal, || Ok(())).unwrap();
        let mut work = ImplicitConductionWorkspace::new(2).unwrap();
        for graph_path in [false, true] {
            let (next, report) = if graph_path {
                graph.advance_to(&thermal, 1.0, controls(1.0), &mut work, || Ok(()))
            } else {
                chain.advance_to(&thermal, 1.0, controls(1.0), &mut work, || Ok(()))
            }.unwrap();
            // Unit masses/Cp and G=1 give contrast 200/(1+2*31.5)=3.125.
            let (hot, cold) = (301.5625_f64, 298.4375_f64);
            near(next.temperature_k(CellId(0)).unwrap(), hot, 1e-9);
            near(next.temperature_k(CellId(1)).unwrap(), cold, 1e-9);
            let conduction_entropy = (hot - cold) * (1.0 / cold - 1.0 / hot);
            let storage_entropy = (hot / 400.0).ln() + (cold / 200.0).ln();
            let entropy = report.entropy;
            let carried = entropy.advection.unwrap();
            near(carried.boundary_entropy_input_j_k, 0.0, 0.0);
            near(carried.upwind_mixing_j_k, 30.5 * conduction_entropy, 1e-12);
            near(carried.balance_residual_j_k, 0.0, 1e-12);
            near(entropy.conduction_generation_j_k, conduction_entropy, 1e-12);
            near(entropy.thermal_storage_change_j_k, storage_entropy, 1e-12);
            near(accounted(entropy), storage_entropy, 1e-12);
            assert!(entropy.temporal_defect_j_k > 0.0);
            assert!(!entropy.has_unclassified_sources);
            assert_eq!(entropy.unclassified_source_at_cell_j_k, 0.0);
            near(report.energy.stored_energy_change_j, 0.0, 1e-9);
            near(report.energy.external_heat_j, 0.0, 0.0);
        }
        assert_eq!(flow.audit().substeps, 0);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn implicit_entropy_open_flow_uses_actual_be_and_sdirk_stage_temperatures() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (mut thermal, conduction) = setup(&law, &[400.0], Some(500.0), vec![]);
    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 20.0;
    let flow = open_flow(1, 2.0, 300.0);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap()
        .with_implicit_advection_entropy(&thermal, || Ok(())).unwrap();
    let mut work = ImplicitConductionWorkspace::new(1).unwrap();
    let dt = 0.4;
    let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
    // The independent scalar equation is T'=2*(300-T)+(500-T)+20.
    let first = (400.0 + gamma * dt * 1120.0) / (1.0 + gamma * dt * 3.0);
    for second_order in [false, true] {
        let expected = if second_order {
            (400.0 + (1.0 - gamma) * dt * (1120.0 - 3.0 * first)
                + gamma * dt * 1120.0) / (1.0 + gamma * dt * 3.0)
        } else {
            (400.0 + dt * 1120.0) / (1.0 + dt * 3.0)
        };
        let (next, report) = if second_order {
            chain.advance_sdirk2_to(&thermal, dt, controls(dt), &mut work, || Ok(()))
        } else {
            chain.advance_to(&thermal, dt, controls(dt), &mut work, || Ok(()))
        }.unwrap();
        near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-10);
        let samples = if second_order {
            [((1.0 - gamma) * dt, first), (gamma * dt, expected)]
        } else {
            [(dt, expected), (0.0, expected)]
        };
        let weighted = |rate: &dyn Fn(f64) -> f64| -> f64 {
            samples.iter().map(|&(weight, t)| weight * rate(t)).sum()
        };
        let boundary = weighted(&|t| 2.0 * (300.0 / t).ln());
        let mixing = weighted(&|t| 2.0 * ((300.0 - t) / t - (300.0 / t).ln()));
        let bath = weighted(&|t| (500.0 - t) / 500.0);
        let generation = weighted(&|t| (500.0 - t) * (1.0 / t - 1.0 / 500.0));
        let source = weighted(&|t| 20.0 / t);
        let storage = (expected / 400.0).ln();
        let entropy = report.entropy;
        let carried = entropy.advection.unwrap();
        near(carried.boundary_entropy_input_j_k, boundary, 1e-12);
        near(carried.upwind_mixing_j_k, mixing, 1e-12);
        near(carried.balance_residual_j_k, 0.0, 1e-12);
        near(entropy.reservoir_entropy_input_j_k, bath, 1e-12);
        near(entropy.conduction_generation_j_k, generation, 1e-12);
        near(entropy.unclassified_source_at_cell_j_k, source, 1e-12);
        near(entropy.temporal_defect_j_k, storage - boundary - mixing - bath - generation - source, 1e-12);
        near(entropy.thermal_storage_change_j_k, storage, 1e-12);
        near(accounted(entropy), storage, 1e-12);
        assert!(entropy.has_unclassified_sources);
        assert!(carried.upwind_mixing_j_k > 0.0);
        near(report.energy.stored_energy_change_j, expected - 400.0, 1e-10);
        near(report.energy.external_heat_j, expected - 400.0, 1e-9);
        assert_eq!(report.energy.substeps, 1);
    }
}

#[test]
fn implicit_entropy_carries_latent_energy_without_isothermal_mixing() {
    let law = Law { nonlinear_k: true, plateau: true };
    let (mut thermal, conduction) = setup(&law, &[350.0, 400.0], None, vec![]);
    thermal.cells.get_mut(&CellId(1)).unwrap().mass_kg = 2.0;
    let flow = open_flow(2, 4.0, 450.0);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap()
        .with_implicit_advection_entropy(&thermal, || Ok(())).unwrap();
    let mut work = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, report) = chain.advance_to(&thermal, 2.0, controls(2.0), &mut work, || Ok(())).unwrap();
    let first = (350.0 + 8.0 * 450.0) / 9.0;
    let second = (2.0 * 400.0 + 8.0 * first) / 10.0;
    let external = 8.0 * (450.0 - second);
    near(next.cells()[&CellId(0)].specific_h_j_kg, first, 1e-10);
    near(next.cells()[&CellId(1)].specific_h_j_kg, second, 1e-10);
    for id in [CellId(0), CellId(1)] {
        assert_eq!(next.temperature_k(id).unwrap(), 300.0);
    }
    let entropy = report.entropy;
    let carried = entropy.advection.unwrap();
    near(carried.boundary_entropy_input_j_k, external / 300.0, 1e-12);
    near(carried.upwind_mixing_j_k, 0.0, 1e-12);
    near(carried.balance_residual_j_k, 0.0, 1e-12);
    near(entropy.temporal_defect_j_k, 0.0, 1e-12);
    near(entropy.thermal_storage_change_j_k, external / 300.0, 1e-12);
    assert_eq!(entropy.conduction_generation_j_k, 0.0);
    assert!(!entropy.has_unclassified_sources);
    near(report.energy.external_heat_j, external, 1e-9);
}

#[test]
fn implicit_entropy_adaptive_account_contains_only_the_fine_trajectory() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (thermal, conduction) = setup(&law, &[400.0], None, vec![]);
    let flow = open_flow(1, 2.0, 300.0);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap()
        .with_implicit_advection_entropy(&thermal, || Ok(())).unwrap();
    let mut work = ImplicitConductionWorkspace::new(1).unwrap();
    let policy = adaptive::AdaptiveImplicitConductionControls {
        step: controls(1.0),
        initial_dt_s: 1.0,
        maximum_trials: 4,
        enthalpy_absolute_tolerance_j_kg: 1000.0,
        temperature_absolute_tolerance_k: 1000.0,
        phase_fraction_absolute_tolerance: None,
        relative_tolerance: 0.0,
    };
    let (next, report) = chain.advance_adaptive_to(
        &thermal, 1.0, policy, &mut work, || Ok(()),
    ).unwrap();
    assert_eq!(report.accepted_comparisons, 1);
    assert_eq!(report.solve.energy.substeps, 2);
    near(next.temperature_k(CellId(0)).unwrap(), 325.0, 1e-10);
    // Two accepted dt=1/2 endpoints are 350 K and 325 K. The coarse 333 1/3 K
    // endpoint is discarded; including its entropy would fail these formulas.
    let boundary = (300.0_f64 / 350.0).ln() + (300.0_f64 / 325.0).ln();
    let mixing = (300.0 - 350.0) / 350.0 + (300.0 - 325.0) / 325.0 - boundary;
    let carried = report.solve.entropy.advection.unwrap();
    near(carried.boundary_entropy_input_j_k, boundary, 1e-12);
    near(carried.upwind_mixing_j_k, mixing, 1e-12);
    near(accounted(report.solve.entropy), (325.0_f64 / 400.0).ln(), 1e-12);
    assert!(!report.solve.entropy.has_unclassified_sources);
}

#[derive(Clone)]
struct CancellingSource;
impl DepositionObserver for CancellingSource {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(BTreeMap::from([(CellId(0), -20.0)]))
    }
}

#[test]
fn implicit_entropy_keeps_real_sources_classified_when_their_watts_cancel() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (mut thermal, conduction) = setup(&law, &[400.0], None, vec![]);
    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.source_w = 20.0;
    let flow = open_flow(1, 2.0, 300.0);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap()
        .with_implicit_advection_entropy(&thermal, || Ok(())).unwrap();
    let mut work = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, _, report) = chain.advance_observed_to(
        &thermal, &CancellingSource, 0.5, controls(0.5), &mut work, || Ok(()),
    ).unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), 350.0, 1e-10);
    assert!(report.entropy.has_unclassified_sources);
    near(report.entropy.unclassified_source_at_cell_j_k, 0.0, 0.0);
    near(report.energy.external_heat_j, -50.0, 1e-9);
    near(accounted(report.entropy), (350.0_f64 / 400.0).ln(), 1e-12);
}

#[test]
fn implicit_entropy_late_cancellation_and_reconstruction_refuse_atomically() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (thermal, conduction) = setup(&law, &[400.0], None, vec![]);
    let flow = open_flow(1, 2.0, 300.0);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap()
        .with_implicit_advection_entropy(&thermal, || Ok(())).unwrap();
    let initial_cells = thermal.cells().clone();
    let mut work = ImplicitConductionWorkspace::new(1).unwrap();
    let mut polls = 0;
    let (expected, report) = chain.advance_sdirk2_to(
        &thermal, 0.5, controls(0.5), &mut work, || { polls += 1; Ok(()) },
    ).unwrap();
    let mut visited = 0;
    assert!(matches!(chain.advance_sdirk2_to(
        &thermal, 0.5, controls(0.5), &mut work, || {
            visited += 1;
            if visited == polls - 1 { Err(DepositionError::Cancelled) } else { Ok(()) }
        },
    ), Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))));
    assert_eq!(thermal.cells(), &initial_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(flow.audit().substeps, 0);
    let (replayed, retry) = chain.advance_sdirk2_to(
        &thermal, 0.5, controls(0.5), &mut work, || Ok(()),
    ).unwrap();
    assert_eq!(replayed.cells(), expected.cells());
    assert_eq!(retry.entropy, report.entropy);

    let (thermal, conduction) = setup(&law, &[200.0, 250.0, 300.0], None, vec![]);
    let flow = open_flow(3, 0.1, 180.0);
    let plan = axial_plan(&thermal);
    let graph = ImplicitConductionGraph::new(
        &thermal, &conduction, &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(), || Ok(()),
    ).unwrap().with_implicit_advection_reconstruction(
        &thermal, &flow, &plan, || Ok(()),
    ).unwrap();
    assert!(matches!(graph.with_implicit_advection_entropy(&thermal, || Ok(())),
        Err(ImplicitConductionError::Invalid(
            "reconstructed enthalpy transport requires its own entropy closure"
        ))));
}
