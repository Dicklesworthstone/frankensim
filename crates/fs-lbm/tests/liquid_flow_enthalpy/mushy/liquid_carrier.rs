//! Moving liquid carries its phase enthalpy through the original solid skeleton.
use super::*;
use fs_lbm::conjugate::ViscousHeatControls;
#[cfg(feature = "implicit-enthalpy-coupling")]
use fs_lbm::conjugate::enthalpy::advection::feedback::implicit::ImplicitFlowControls;

#[path = "liquid_carrier/open.rs"]
mod open;

fn state<'a>(domain: &VoxelDomain, curve: &'a EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'a> {
    let events = [20.0, 40.0, 60.0, 80.0].into_iter().enumerate().map(|(i, h)| {
        let [x, y, z] = domain.coords(i);
        let center = domain.center(x, y, z);
        let mut start = center;
        let mut end = center;
        start[0] -= 0.5;
        end[0] += 0.5;
        TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(i as u64), segment_id: i as u64, body_id: 0, material_id: MAT,
            mass_kg: 1.0, specific_h_j_kg: h,
            geometry: CellGeometry { start_m: start, end_m: end, volume_m3: 1.0, surface_area_m2: 6.0 },
            boundary: Boundary { ambient_k: 500.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }) }
    }).collect();
    DepositionSimulation::new(events, &[(MAT, curve as &dyn ThermalLaw)], DepositionBudget {
        max_cells: 4, max_links: 4, max_events: 4, max_substeps: 100,
    }).unwrap().advance_to(0.0, thermal(), || Ok(())).unwrap().0
}

fn cavity() -> SimpleConfig {
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall { velocity: [1.0, 0.0, 0.0] };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-3;
    config
}
fn heating() -> ViscousHeatControls {
    ViscousHeatControls { maximum_reynolds: 0.1, power_tolerance_w: 1e-7 }
}
fn thermal() -> StepConfig { StepConfig { max_dt_s: 10.0, energy_tolerance_j: 1e-8 } }
fn flow_input<'a>(ids: &'a [Option<CellId>], curve: &'a EquilibriumEnthalpyPhaseCurve) -> FlowEnthalpyInput<'a> {
    FlowEnthalpyInput { cell_ids: ids, material_id: MAT, caloric: curve,
        inlet_specific_h_j_kg: [None; 6], max_cells: 4, max_connections: 64 }
}

// Independent constrained four-face Stokes solution for mu=100 and unit lid
// speed. Every nonzero face flux is +/-q; no phase fraction is transported.
fn circulation(h: [f64; 4]) -> (f64, [f64; 4]) {
    let beta = h.map(|h| { let f = h / 100.0; 100.0 * (1.0 - f).powi(2) / f.powi(3) });
    let q = 200.0 / (2400.0 + beta.iter().sum::<f64>());
    let watts = std::array::from_fn(|i| {
        let shear = if i < 2 { 600.0 * q * q }
            else { 500.0 * q * q + 100.0 * (1.0 - q).powi(2) };
        shear + beta[i] * q * q
    });
    (q, watts)
}

#[test]
fn native_liquid_circulation_only_melts_stationary_solid_by_physical_wall_heat() {
    let card = continuity::source(100.0);
    let curve = continuity::chart(&card);
    let phase = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = cavity();
    let gate = CancelGate::new();
    let ids = std::array::from_fn::<_, 4, _>(|i| Some(CellId(i as u64)));
    let initial = state(&domain, &curve);
    let base = MaterialFlowFeedback::try_new_phase(
        &initial, &domain, flow_input(&ids, &curve), &config, (), &gate, 20, &phase,
        |_, cell: &DepositedCell, _| continuity::resolve(&card, &phase, cell.specific_h_j_kg),
    ).unwrap();
    assert!(base.clone().with_liquid_carrier().is_err());
    let mut observer = base.clone().with_viscous_heating(heating()).unwrap().with_liquid_carrier().unwrap();
    assert!(observer.clone().with_unsteady_momentum().is_err());
    assert!(base.with_unsteady_momentum().unwrap().with_liquid_carrier().is_err());
    assert_eq!(observer.liquid_carrier().unwrap().caloric_identity(), curve.caloric_identity().unwrap());
    let mut current = initial.clone();
    let mut h = [20.0, 40.0, 60.0, 80.0];
    for target in [0.05, 0.1] {
        let (q, watts) = circulation(h);
        let (next, accepted, audit) = current.advance_observed(&observer, target, thermal(), || Ok(())).unwrap();
        assert_eq!(accepted.momentum_solves(), 1);
        near(accepted.applied().unwrap().flow.field.flux_x(1, 1, 0), q);
        assert!(accepted.audit().transported_mass_kg > 0.0);
        assert_eq!(accepted.audit().external_enthalpy_j, 0.0);
        for i in 0..4 {
            h[i] += 0.05 * watts[i];
            near(next.cells()[&CellId(i as u64)].specific_h_j_kg, h[i]);
            assert_eq!(next.temperature_k(CellId(i as u64)).unwrap(), 500.0);
        }
        let work = accepted.viscous_heating_audit().unwrap();
        near(work.generated_heat_j, 0.05 * watts.iter().sum::<f64>());
        near(work.supplied_wall_work_j, work.generated_heat_j);
        near(audit.stored_energy_change_j, work.generated_heat_j);
        near(audit.external_heat_j, work.generated_heat_j);
        assert_eq!(next.mass_kg(), initial.mass_kg());
        assert!(accepted.liquid_carrier().is_some());
        current = next;
        observer = accepted;
    }
    assert!(observer.with_liquid_carrier().is_err());
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
#[cfg(feature = "implicit-enthalpy-coupling")]
fn implicit_phase_transport_retains_liquid_carriage_at_every_flow_endpoint() {
    let card = continuity::source(100.0);
    let curve = continuity::chart(&card);
    let phase = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = cavity();
    let gate = CancelGate::new();
    let ids = std::array::from_fn::<_, 4, _>(|i| Some(CellId(i as u64)));
    let initial = state(&domain, &curve);
    let observer = MaterialFlowFeedback::try_new_phase(
        &initial, &domain, flow_input(&ids, &curve), &config, (), &gate, 2048, &phase,
        |_, cell: &DepositedCell, _| continuity::resolve(&card, &phase, cell.specific_h_j_kg),
    ).unwrap().with_viscous_heating(heating()).unwrap().with_liquid_carrier().unwrap();
    let controls = ImplicitFlowControls {
        max_unknowns: 4, max_newton_iterations: 32, difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-7, total_energy_tolerance_j: 4e-7,
        newton: fs_solver::NewtonKrylovConfig { absolute_tolerance: 2e-8,
            relative_tolerance: 0.0, linear_restart: 4, ..Default::default() },
    };
    let dt = 0.05;
    let (next, report) = observer.step_implicit(&initial, dt, controls, &[0.0; 4], || Ok(())).unwrap();
    let h = std::array::from_fn(|i| next.cells()[&CellId(i as u64)].specific_h_j_kg);
    let (q, watts) = circulation(h);
    near(report.endpoint.flow.field.flux_x(1, 1, 0), q);
    for i in 0..4 {
        let increment = h[i] - initial.cells()[&CellId(i as u64)].specific_h_j_kg;
        assert!((increment - dt * watts[i]).abs() < 1e-7);
        assert_eq!(next.temperature_k(CellId(i as u64)).unwrap(), 500.0);
    }
    near(report.mechanical.generated_heat_j, dt * watts.iter().sum::<f64>());
    near(report.thermal.external_heat_j, report.mechanical.generated_heat_j);
    assert_eq!(report.net_stream_enthalpy_j, 0.0);
    assert!(report.combined_residual_j.abs() < 4e-7);
    assert_eq!(initial.time_s(), 0.0);
    assert!(observer.applied().is_none());
}
