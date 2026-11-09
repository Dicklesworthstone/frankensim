//! G1/G3/G4 physical input response, with independent original-equation oracles.
use super::*;
use super::super::equilibrium::ThermalOperatingControls;
use super::super::harmonic::{ThermalImpedanceBudget, driven::ThermalBathChannel};
use crate::enthalpy::advection::{Advected, MassBoundary, MassFlow};
use crate::enthalpy::advection::schedule::{AdvectionSchedule, AdvectionScheduleKnot};
use fs_math::c64::C64;

fn budget() -> ThermalImpedanceBudget {
    ThermalImpedanceBudget {
        maximum_cells: 8, maximum_dense_bytes: 100_000, maximum_dense_work: 100_000,
        equilibrium_tolerance_w: 1e-9, relative_residual_tolerance: 1e-12,
    }
}
fn flow(cells: usize, rate: f64, h: f64) -> Advected<()> {
    Advected::new((), (0..cells-1).map(|i| MassFlow {
        id: i as u64, a: CellId(i as u64), b: CellId(i as u64+1), kg_s: rate,
    }).collect(), cells+1).unwrap().with_boundaries(vec![
        MassBoundary::Inlet {
            id: cells as u64, cell: CellId(0), kg_s: rate,
            material_id: MaterialId(0), specific_h_j_kg: h,
        },
        MassBoundary::Outlet {
            id: cells as u64+1, cell: CellId(cells as u64-1), kg_s: rate,
        },
    ], cells+1).unwrap()
}
fn close(a: C64, b: C64) {
    assert!((a-b).abs() <= 2e-12*(1.0+b.abs()), "{a:?} != {b:?}");
}
fn powers(thermal: &mut DepositionSimulation<'_>, values: &[f64]) {
    for (i, &value) in values.iter().enumerate() {
        thermal.cells.get_mut(&CellId(i as u64)).unwrap().boundary.source_w = value;
    }
}

#[test]
fn g1_physical_inlet_and_common_pump_share_the_original_directed_resolvent() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (mut thermal, conduction) = setup(&law, &[310.0, 320.0], None, vec![]);
    powers(&mut thermal, &[10.0, 30.0]);
    let before = thermal.cells().clone();
    let flow = flow(2, 2.0, 300.0);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap();
    let model = chain.prepare_driven_thermal_impedance(
        &thermal, budget(), 1, || Ok(()),
    ).unwrap();
    assert_eq!(model.common_flow_input_w(), &[-20.0, -20.0]);
    assert_eq!(model.inlets()[0].boundary_id, 2);
    assert_eq!(model.inlets()[0].mass_rate_kg_s, 2.0);
    let input = BTreeMap::from([(2, C64::new(3.0, 2.0))]);
    for omega in [0.0, 0.2, 4.0] {
        let result = model.factor(omega, || Ok(())).unwrap().response(
            &[C64::ONE, C64::ZERO], &BTreeMap::new(), &input,
            C64::new(0.1, 0.05), || Ok(()),
        ).unwrap();
        // K=[[3,-1],[-3,3]]. Inlet gives (6+4i,0), pump
        // gives (-2-i,-2-i), and only one additional watt is supplied.
        let a = C64::new(3.0, -omega);
        let determinant = a*a-C64::new(3.0, 0.0);
        let p = [C64::new(5.0, 3.0), C64::new(-2.0, -1.0)];
        assert_eq!(result.forcing_w, p);
        close(result.field.energy_j[0], (a*p[0]+p[1])/determinant);
        close(result.field.energy_j[1], (p[0].scale(3.0)+a*p[1])/determinant);
        assert!(result.field.relative_residual <= 1e-12);
        assert_eq!(result.field.temperature_k, result.field.energy_j);
    }
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(flow.audit().substeps, 0);
}

#[test]
fn g1_bath_and_source_drives_use_original_mass_and_surface_conductance() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (mut thermal, conduction) = setup(&law, &[300.0], Some(300.0), vec![]);
    thermal.cells.get_mut(&CellId(0)).unwrap().mass_kg = 2.0;
    // One square metre remains exposed: its original htc is 2 W/(m2 K).
    thermal.cells.get_mut(&CellId(0)).unwrap().boundary.htc_w_m2_k = 2.0;
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap();
    let model = chain.prepare_driven_thermal_impedance(
        &thermal, budget(), 0, || Ok(()),
    ).unwrap();
    let baths = BTreeMap::from([
        ((CellId(0), ThermalBathChannel::Ambient), C64::new(1.0, 0.5)),
        ((CellId(0), ThermalBathChannel::Reservoir), C64::new(2.0, -1.0)),
    ]);
    for omega in [0.0, 2.0] {
        let result = model.factor(omega, || Ok(())).unwrap().response(
            &[C64::ONE], &baths, &BTreeMap::new(), C64::ZERO, || Ok(()),
        ).unwrap();
        // Total G=3, m*Cp=2. Forcing=1 + 2*(1+0.5i) + (2-i)=5.
        assert_eq!(result.forcing_w, [C64::new(5.0, 0.0)]);
        let expected = C64::new(5.0, 0.0)/C64::new(1.5, -omega);
        close(result.field.energy_j[0], expected);
        close(result.field.temperature_k[0], expected.scale(0.5));
    }
}

#[test]
fn g1_scheduled_flow_uses_the_same_frozen_clock_inlet_and_rate() {
    let law = Law { nonlinear_k: false, plateau: false };
    let (mut thermal, conduction) = setup(&law, &[450.0], None, vec![]);
    powers(&mut thermal, &[50.0]);
    thermal.cells.get_mut(&CellId(0)).unwrap().mass_kg = 2.0;
    let flow = flow(1, 2.0, 300.0);
    let schedule = AdvectionSchedule::new(vec![
        AdvectionScheduleKnot {
            time_s: 0.0, flow_multiplier: 0.5,
            inlet_enthalpies_j_kg: BTreeMap::from([(1, 400.0)]),
        },
        AdvectionScheduleKnot {
            time_s: 1.0, flow_multiplier: 2.0,
            inlet_enthalpies_j_kg: BTreeMap::from([(1, 600.0)]),
        },
    ], 2, 1, || Ok(())).unwrap();
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_scheduled_implicit_advection(
        &thermal, &flow, &schedule, || Ok(()),
    ).unwrap();
    let point = chain.solve_operating_point(
        &thermal, ThermalOperatingControls {
            maximum_newton_iterations: 8, maximum_backtracks: 8,
            maximum_residual_evaluations: 32, cell_tolerance_w: 1e-10,
            total_tolerance_w: 1e-10,
        }, &mut ImplicitConductionWorkspace::new(1).unwrap(), || Ok(()),
    ).unwrap();
    let model = point.prepare_driven_thermal_impedance(budget(), 1, || Ok(())).unwrap();
    assert_eq!(model.original().operating_time_s(), 0.0);
    assert_eq!(model.inlets()[0].mass_rate_kg_s, 1.0);
    assert_eq!(model.inlets()[0].specific_enthalpy_j_kg, 400.0);
    assert_eq!(model.common_flow_input_w(), &[-50.0]);
    let result = model.factor(3.0, || Ok(())).unwrap().response(
        &[C64::ZERO], &BTreeMap::new(),
        &BTreeMap::from([(1, C64::new(6.0, -2.0))]),
        C64::new(0.2, 0.1), || Ok(()),
    ).unwrap();
    let energy = C64::new(-4.0, -7.0)/C64::new(0.5, -3.0);
    close(result.field.energy_j[0], energy);
    close(result.field.temperature_k[0], energy.scale(0.5));
    assert_eq!(flow.audit().substeps, 0);
    assert_eq!(schedule.next_transition_time_s(0.0), Some(1.0));
}

#[test]
fn g1_inlet_and_pump_can_change_latent_energy_with_zero_temperature_response() {
    let law = Law { nonlinear_k: false, plateau: true };
    let (mut thermal, conduction) = setup(&law, &[360.0, 380.0], None, vec![]);
    powers(&mut thermal, &[20.0, 40.0]);
    let flow = flow(2, 2.0, 350.0);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap();
    let model = chain.prepare_driven_thermal_impedance(
        &thermal, budget(), 1, || Ok(()),
    ).unwrap();
    assert_eq!(model.common_flow_input_w(), &[-20.0, -40.0]);
    for omega in [0.0, 4.0] {
        let result = model.factor(omega, || Ok(())).unwrap().response(
            &[C64::ZERO; 2], &BTreeMap::new(),
            &BTreeMap::from([(2, C64::new(4.0, 0.0))]),
            C64::new(0.1, 0.0), || Ok(()),
        ).unwrap();
        let a = C64::new(2.0, -omega);
        let first = C64::new(6.0, 0.0)/a;
        close(result.field.energy_j[0], first);
        close(result.field.energy_j[1], (first.scale(2.0)-C64::new(4.0, 0.0))/a);
        assert_eq!(result.field.temperature_k, [C64::ZERO; 2]);
    }
}

struct ShiftedLaw { original: Law, datum: f64 }
impl ThermalLaw for ShiftedLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.original.temperature_k(h-self.datum)
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
}
impl ConductivityLaw for ShiftedLaw {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.original.conductivity_w_m_k(h-self.datum)
    }
}
impl ImplicitConductionLaw for ShiftedLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw { self }
    fn conductivity_law(&self) -> &dyn ConductivityLaw { self }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] { [100.0+self.datum, 1000.0+self.datum] }
    fn state_at_specific_enthalpy(&self, h: f64)
        -> Result<ImplicitConductionState, DepositionError> {
        self.original.state_at_specific_enthalpy(h-self.datum)
    }
}
#[test]
fn g3_common_rate_and_inlet_columns_do_not_depend_on_formation_datum() {
    let mut results = Vec::new();
    for datum in [0.0, 1_099_511_627_776.0] {
        let law = ShiftedLaw {
            original: Law { nonlinear_k: false, plateau: false }, datum,
        };
        let (mut thermal, conduction) = setup(&law, &[310.0+datum, 320.0+datum], None, vec![]);
        powers(&mut thermal, &[10.0, 30.0]);
        let flow = flow(2, 2.0, 300.0+datum);
        let chain = ImplicitConductionChain::new(
            &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
        ).unwrap().with_implicit_advection(&thermal, &flow, || Ok(())).unwrap();
        let model = chain.prepare_driven_thermal_impedance(
            &thermal, budget(), 1, || Ok(()),
        ).unwrap();
        results.push(model.factor(0.2, || Ok(())).unwrap().response(
            &[C64::ZERO; 2], &BTreeMap::new(),
            &BTreeMap::from([(2, C64::new(3.0, -1.0))]),
            C64::new(0.125, 0.0), || Ok(()),
        ).unwrap());
    }
    assert_eq!(results[0], results[1]);
}

struct IntegralBoundary<'a>(&'a Law);
impl boundary::ImplicitBoundaryLaw for IntegralBoundary<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw { self.0 }
    fn state_at_specific_enthalpy(
        &self, h: f64, bath: f64, geometry: boundary::ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<boundary::ImplicitBoundaryState, DepositionError> {
        checkpoint()?;
        let t = self.0.temperature_k(h)?;
        // Exact integral of k(T)=T/100 across this finite half-cell.
        let scale = geometry.area_m2/(100.0*geometry.half_distance_m);
        Ok(boundary::ImplicitBoundaryState {
            conductance_w_k: scale*(t+bath)*0.5,
            power_enthalpy_derivative_w_kg_j: scale*t,
        })
    }
}
#[test]
fn g4_original_domains_cancellation_budget_and_unknown_input_derivatives_refuse() {
    let law = Law { nonlinear_k: true, plateau: false };
    let (thermal, conduction) = setup(&law, &[300.0], Some(300.0), vec![]);
    let boundary_law = IntegralBoundary(&law);
    let chain = ImplicitConductionChain::new(
        &thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_boundary_laws(
        &thermal, &[(CellId(0), &boundary_law)], &[], || Ok(()),
    ).unwrap();
    let model = chain.prepare_driven_thermal_impedance(
        &thermal, budget(), 0, || Ok(()),
    ).unwrap();
    let baths = BTreeMap::from([(
        (CellId(0), ThermalBathChannel::Reservoir), C64::ONE,
    )]);
    assert!(model.forcing(&[C64::ZERO], &baths, &BTreeMap::new(), C64::ZERO, || Ok(())).is_err());
    let response = model.factor(0.0, || Ok(())).unwrap().response(
        &[C64::ONE], &BTreeMap::new(), &BTreeMap::new(), C64::ZERO, || Ok(()),
    ).unwrap();
    close(response.field.temperature_k[0], C64::new(1.0/3.0, 0.0));
    assert!(model.forcing(&[C64::ZERO], &BTreeMap::new(), &BTreeMap::new(), C64::ONE, || Ok(())).is_err());

    let flow = flow(1, 2.0, 300.0);
    let chain = chain.with_implicit_advection(&thermal, &flow, || Ok(())).unwrap();
    assert!(chain.prepare_driven_thermal_impedance(&thermal, budget(), 0, || Ok(())).is_err());
    let model = chain.prepare_driven_thermal_impedance(&thermal, budget(), 1, || Ok(())).unwrap();
    for inlet in [
        BTreeMap::from([(9, C64::ONE)]),
        BTreeMap::from([(1, C64::new(500.0, 0.0))]),
    ] {
        assert!(model.forcing(&[C64::ZERO], &BTreeMap::new(), &inlet, C64::ZERO, || Ok(())).is_err());
    }
    assert!(model.forcing(&[C64::ZERO], &BTreeMap::new(), &BTreeMap::new(), C64::ONE, || Ok(())).is_err());
    let mut calls = 0;
    let inputs = BTreeMap::from([(1, C64::ONE)]);
    let expected = model.forcing(&[C64::ZERO], &BTreeMap::new(), &inputs, C64::ZERO, || {
        calls += 1; Ok(())
    }).unwrap();
    for stop in [1, calls/2, calls] {
        let mut count = 0;
        assert!(model.forcing(&[C64::ZERO], &BTreeMap::new(), &inputs, C64::ZERO, || {
            count += 1;
            if count == stop { Err(DepositionError::Cancelled) } else { Ok(()) }
        }).is_err());
    }
    assert_eq!(model.forcing(&[C64::ZERO], &BTreeMap::new(), &inputs, C64::ZERO, || Ok(())).unwrap(), expected);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), 300.0);
    assert_eq!(flow.audit().substeps, 0);
}
