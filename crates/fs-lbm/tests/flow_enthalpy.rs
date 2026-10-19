#![cfg(feature = "enthalpy-coupling")]
//! Actual material-driven momentum, caloric storage and shared source history.
//! Synthetic constants exercise the coupling, not measured material validation.

use std::collections::BTreeMap;

use fs_blake3::ContentHash;
use fs_conduction::transient::deposition::*;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    FacePatch, FanCurve, FvBoundary, InternalFan, SimpleConfig, VoxelDomain,
};
use fs_lbm::conjugate::enthalpy::advection::FlowEnthalpyInput;
use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

const MAT: MaterialId = MaterialId(3);

fn curve(rho: f64, latent: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0xAD; 32]),
        [(0.0, 300.0, 0.0), (20.0, 310.0, 0.0),
         (20.0 + latent, 310.0, 1.0), (1020.0 + latent, 810.0, 1.0)]
            .map(|(h, t, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h, temperature_k: t,
                liquid_mass_fraction: f, bulk_density_kg_m3: rho,
            }).to_vec()).unwrap()
}

fn config() -> StepConfig { StepConfig { max_dt_s: 0.125, energy_tolerance_j: 1e-9 } }

fn initial<'a>(domain: &VoxelDomain, rho: f64, law: &'a dyn ThermalLaw) -> DepositionSimulation<'a> {
    let events = (0..domain.cell_count()).map(|c| {
        let [x, y, z] = domain.coords(c);
        let centre = domain.center(x, y, z);
        let widths = domain.widths(c);
        let mut start = centre;
        let mut end = centre;
        start[0] -= 0.5 * widths[0];
        end[0] += 0.5 * widths[0];
        TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(c as u64), segment_id: c as u64, body_id: 1,
            material_id: MAT, mass_kg: rho * domain.volume(c), specific_h_j_kg: 0.0,
            geometry: CellGeometry { start_m: start, end_m: end,
                volume_m3: domain.volume(c), surface_area_m2: 6.0 },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0,
                reservoir: None, source_w: 0.0 },
        }) }
    }).collect();
    DepositionSimulation::new(events, &[(MAT, law)], DepositionBudget {
        max_cells: domain.cell_count(), max_links: 64,
        max_events: domain.cell_count(), max_substeps: 10000,
    }).unwrap().advance_to(0.0, config(), || Ok(())).unwrap().0
}

fn input<'a>(ids: &'a [Option<CellId>], law: &'a dyn ThermalLaw) -> FlowEnthalpyInput<'a> {
    FlowEnthalpyInput { cell_ids: ids, material_id: MAT, caloric: law,
        inlet_specific_h_j_kg: [Some(0.0); 6], max_cells: ids.len(), max_connections: 64 }
}

fn throughflow() -> SimpleConfig {
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet { velocity: [0.25, 0.0, 0.0] };
    faces[1] = FvBoundary::Outlet;
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config
}

#[derive(Clone, Default, Debug, PartialEq)]
struct Heater {
    durations: Vec<f64>,
    supplied_j: f64,
    reject_step: Option<usize>,
}
impl DepositionObserver for Heater {
    type Error = DepositionError;
    fn sources(&mut self, _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<BTreeMap<CellId, f64>, Self::Error>
    {
        Ok(BTreeMap::from([(CellId(1), 100.0)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> { Ok(()) }
    fn step(&mut self, _: f64, dt: f64, _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>, _: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error>
    {
        self.durations.push(dt);
        self.supplied_j += dt * 100.0;
        if self.reject_step == Some(self.durations.len()) {
            return Err(DepositionError::InvalidInput("independent history refused endpoint"));
        }
        Ok(())
    }
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-8 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn actual_solved_face_flux_transports_latent_enthalpy_and_keeps_source_work_separate() {
    // This test exercises a declared passive, equal-density caloric plateau.
    // It does not assert a partly solid material's physical flow mobility.
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(2.0, 2.0, 0.01, 0.1).unwrap();
    let fluid = MaterialChtFluid::try_new(&source).unwrap();
    let gate = CancelGate::new();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let law = curve(2.0, 100.0);
    let old = initial(&domain, 2.0, &law);
    let ids: Vec<_> = (0..3).map(|c| Some(CellId(c))).collect();
    let observer = solved.advect_enthalpy(&old, input(&ids, &law), Heater::default(), &gate).unwrap();
    let (next, accepted, audit) = old.advance_observed(&observer, 1.0, config(), || Ok(())).unwrap();
    assert!(core::ptr::eq(accepted.material(), &source));
    assert!(core::ptr::eq(accepted.field(), &solved.flow().field));
    assert_eq!(accepted.transport(), solved.transport());
    for flow in accepted.advection().flows() { near(flow.kg_s, 0.5); }
    // Independent donor-cell recursion, driven by the SAME accepted durations.
    // A Kelvin update would stall at the plateau and miss its stored energy.
    let mut h = [0.0; 3];
    let mut exported_j = 0.0;
    for &dt in &accepted.advection().inner().durations {
        exported_j += dt * 0.5 * h[2];
        let old_h = h;
        for c in 0..3 {
            let upstream = if c == 0 { 0.0 } else { old_h[c - 1] };
            let heat = if c == 1 { 100.0 } else { 0.0 };
            h[c] += dt * (heat + 0.5 * (upstream - old_h[c])) / 2.0;
        }
    }
    for (c, expected) in h.into_iter().enumerate() {
        near(next.cells()[&CellId(c as u64)].specific_h_j_kg, expected);
    }
    let flow_audit = accepted.advection().audit();
    near(flow_audit.external_enthalpy_j, -exported_j);
    near(flow_audit.inlet_mass_kg, 0.5);
    near(flow_audit.outlet_mass_kg, 0.5);
    near(audit.external_heat_j, accepted.advection().inner().supplied_j - exported_j);
    near(next.energy_j() - old.energy_j(), audit.external_heat_j);
    near(next.temperature_k(CellId(1)).unwrap(), 310.0);
    let phase = law.state_at_specific_enthalpy(next.cells()[&CellId(1)].specific_h_j_kg).unwrap();
    assert!(phase.liquid_mass_fraction() > 0.0 && phase.liquid_mass_fraction() < 1.0);
    near(next.mass_kg(), old.mass_kg());
    assert_eq!(observer.advection().audit().substeps, 0);
    assert_eq!(observer.advection().inner(), &Heater::default());
}

#[test]
fn material_viscosity_changes_solved_flow_and_advected_heat_without_a_second_flux_model() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::wall(); 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut config_flow = SimpleConfig::new(faces);
    config_flow.tolerance = 1e-13;
    config_flow.pressure_tolerance = 1e-12;
    config_flow.internal_fans.push(InternalFan {
        patch: FacePatch { axis: 0, index: 2, lo: [0, 0, 0], hi: [4, 1, 1] },
        blows_positive: true,
        curve: FanCurve::new(&[(0.0, 8.0), (100.0, 8.0)]).unwrap(),
    });
    let law = curve(1.0, 100.0);
    let old = initial(&domain, 1.0, &law);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let mut outputs = Vec::new();
    for mu in [1.0, 4.0] {
        let source = FluidTransportState::try_new(1.0, 2.0, 0.01, mu).unwrap();
        let fluid = MaterialChtFluid::try_new(&source).unwrap();
        let solved = fluid.solve_flow(&domain, &config_flow, &gate).unwrap();
        let observer = solved.advect_enthalpy(&old, input(&ids, &law), Heater::default(), &gate).unwrap();
        let (next, accepted, audit) = old.advance_observed(&observer, 1.0, config(), || Ok(())).unwrap();
        let q = solved.flow().field.flux_x(2, 0, 0);
        assert!(q > 0.0);
        for edge in accepted.advection().flows() { near(edge.kg_s, q); }
        near(next.energy_j() - old.energy_j(), audit.external_heat_j);
        near(audit.external_heat_j, 100.0 + accepted.advection().audit().external_enthalpy_j);
        outputs.push((q, next.cells()[&CellId(1)].specific_h_j_kg));
    }
    assert!(outputs[1].0 < outputs[0].0); // More viscous fluid flows less.
    assert!(outputs[1].1 > outputs[0].1); // Same heater loses less advected heat.
}

#[test]
fn foreign_caloric_law_density_map_and_rejected_history_preserve_accepted_state() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(1.0, 2.0, 0.01, 0.1).unwrap();
    let fluid = MaterialChtFluid::try_new(&source).unwrap();
    let gate = CancelGate::new();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let law = curve(1.0, 100.0);
    let foreign = curve(1.0, 200.0); // Same initial T and cell identity are insufficient.
    let old = initial(&domain, 1.0, &law);
    let ids: Vec<_> = (0..3).map(|c| Some(CellId(c))).collect();
    assert!(solved.advect_enthalpy(&old, input(&ids, &foreign), (), &gate).is_err());
    let wrong_mass = initial(&domain, 2.0, &law);
    assert!(solved.advect_enthalpy(&wrong_mass, input(&ids, &law), (), &gate).is_err());
    let duplicated = [Some(CellId(0)), Some(CellId(0)), Some(CellId(2))];
    assert!(solved.advect_enthalpy(&old, input(&duplicated, &law), (), &gate).is_err());
    let history = Heater { reject_step: Some(2), ..Heater::default() };
    let observer = solved.advect_enthalpy(&old, input(&ids, &law), history.clone(), &gate).unwrap();
    assert!(old.advance_observed(&observer, 1.0, config(), || Ok(())).is_err());
    assert_eq!(observer.advection().inner(), &history);
    assert_eq!(observer.advection().audit().substeps, 0);
    near(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
    let result = old.advance_observed(&observer, 1.0, config(), || Err(DepositionError::Cancelled));
    assert!(matches!(result, Err(DepositionError::Cancelled)));
    assert_eq!(observer.advection().inner(), &history);
}
