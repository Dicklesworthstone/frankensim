//! Authored numerical component/load fixtures; not measured ceramic calibration.
use super::{near, rod};
#[path = "../circuit_enthalpy/common.rs"]
mod thermal_source;

#[path = "circuit/thermal.rs"]
mod ceramic_thermal;

use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve};
use fs_couple::circuit_enthalpy::material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::circuit_enthalpy::piezoelectric::{self as heated, PiezoelectricEnthalpyInput};
use fs_couple::material_circuit::{MaterialCircuitError, MaterialCircuitSolve, MaterialResistor};
use fs_couple::material_circuit::piezoelectric::PiezoelectricCircuit;
use fs_material::conductor::{ElectricalPropertyBasis, resolve_uniform_conductor};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_thermal::enthalpy::{CellId, MaterialId};

fn source(branch: usize, resistance: f64) -> MaterialResistor {
    let (card, _) = thermal_source::material(resistance, 1.0, 0.0);
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor(&card, &thermal_source::point(300.0),
            MaterialPropertySelection::SingleClaimOnly, 1.0, 1.0).unwrap(),
    }
}
fn solve() -> MaterialCircuitSolve {
    MaterialCircuitSolve { energy_tolerance_j: 1e-10, power_tolerance_w: 1e-10 }
}
fn loaded(d: f64, density: f64, resistance: f64, x: f64, v: f64, q: f64)
    -> PiezoelectricCircuit {
    PiezoelectricCircuit::new(CircuitGraph { node_count: 2, branches: vec![
        (1,0,Branch::Capacitor { farads: f64::NAN }),
        (1,0,Branch::Resistor { ohms: f64::NAN }),
    ], transformers: vec![] }, vec![source(1,resistance)], vec![], vec![q],
        0, rod(d,density), x, v, &[], || false).unwrap()
}

// Hand elimination of the three midpoint constitutive equations:
// Qdot=-V/R, m vdot=F-k_E*x+theta*V, Q=C*V+theta*x.
// Coefficients below come independently from the authored source values.
fn load_oracle(x: f64, v: f64, q: f64, d: f64, density: f64,
    resistance: f64, force: f64, h: f64) -> (f64,f64,f64,f64) {
    let mass = density/3.0;
    let stiffness = 2.5;
    let capacitance = 0.25*(2.0-10.0*d*d);
    let theta = 2.5*d;
    let denominator = 1.0+h*h*stiffness/(4.0*mass);
    let a = (x+0.5*h*v+h*h*force/(4.0*mass))/denominator;
    let b = h*h*theta/(4.0*mass*denominator);
    let voltage = (q-theta*a)/(capacitance+h/(2.0*resistance)+theta*b);
    let midpoint_x = a+b*voltage;
    (2.0*midpoint_x-x, 4.0*(midpoint_x-x)/h-v, q-h*voltage/resistance, voltage)
}
fn energy(x: f64, v: f64, q: f64, d: f64, density: f64) -> f64 {
    let capacitance = 0.25*(2.0-10.0*d*d);
    0.5*(density/3.0)*v*v + 1.25*x*x + (q-2.5*d*x).powi(2)/(2.0*capacitance)
}

#[test]
fn g1_source_substitution_drives_actual_resistive_current_motion_and_heat() {
    let mut endpoints = Vec::new();
    for (d,density,resistance) in [
        (0.2,3.0,1.0), (-0.2,3.0,1.0), (0.4,3.0,1.0),
        (0.2,6.0,1.0), (0.2,3.0,4.0), (0.0,3.0,1.0),
    ] {
        let (mut x,mut v,mut q) = (0.02,0.01,if d < 0.0 { -0.02 } else { 0.02 });
        let mut body = loaded(d,density,resistance,x,v,q);
        let initial = body.clone();
        let mut total_heat = 0.0;
        let mut total_work = 0.0;
        let initial_energy = energy(x,v,q,d,density);
        for _ in 0..8 {
            let expected = load_oracle(x,v,q,d,density,resistance,0.03,0.05);
            let old = body.clone();
            let next = body.step_candidate(vec![source(1,resistance)], &[], 0.03,
                0.05, solve(), || false).unwrap();
            let report = next.step();
            near(next.model().extension_m(),expected.0,2e-10);
            near(next.model().velocity_m_s(),expected.1,2e-9);
            near(next.model().charge_c(),expected.2,2e-10);
            near(report.resistors[0].voltage_v,expected.3,2e-10);
            near(report.resistors[0].current_a,expected.3/resistance,2e-10);
            near(report.heat_j,0.05*expected.3*expected.3/resistance,2e-10);
            near(report.stored_after_j,energy(expected.0,expected.1,expected.2,d,density),2e-10);
            near(report.mechanical_supply_j,0.03*(expected.0-x),2e-10);
            near(report.electrical_supply_j,0.0,0.0);
            near(report.other_electrical_energy_j[0],0.0,0.0);
            near(report.other_electrical_energy_j[1],0.0,0.0);
            total_heat += report.heat_j; total_work += report.mechanical_supply_j;
            near(report.stored_after_j+total_heat-initial_energy-total_work,0.0,2e-10);
            assert_eq!(body,old);
            assert_eq!(next.model().rod().material(),initial.rod().material());
            (x,v,q) = (expected.0,expected.1,expected.2);
            body = next.into_parts().0;
        }
        endpoints.push((body.extension_m(),body.charge_c(),total_heat));
    }
    // Reversing poling and charge reverses voltage/current, retaining mechanics
    // and positive Joule heat on the SAME native physical trajectory.
    near(endpoints[0].0,endpoints[1].0,2e-10);
    near(endpoints[0].1,-endpoints[1].1,2e-10);
    near(endpoints[0].2,endpoints[1].2,2e-10);
    assert!((endpoints[0].0-endpoints[3].0).abs()>1e-6);
    assert!((endpoints[0].1-endpoints[4].1).abs()>1e-5);
}

#[test]
fn g1_other_charge_and_flux_stores_retain_native_kirchhoff_exchange() {
    let (x,q,flux) = (0.01,0.05,0.01);
    let voltage = (q-0.5*x)/0.4;
    let graph = CircuitGraph { node_count: 2, branches: vec![
        (1,0,Branch::Inductor { henries: 0.5 }),
        (1,0,Branch::Capacitor { farads: 0.3 }),
        (1,0,Branch::Resistor { ohms: f64::NAN }),
        (1,0,Branch::Capacitor { farads: f64::NAN }),
    ], transformers: vec![] };
    let mut body = PiezoelectricCircuit::new(graph,vec![source(2,2.0)],
        vec![flux],vec![0.3*voltage,q],3,rod(0.2,3.0),x,0.0,&[],||false).unwrap();
    for _ in 0..5 {
        let next = body.step_candidate(vec![source(2,2.0)],&[],0.02,0.03,solve(),||false).unwrap();
        let proposed = next.model();
        let report = next.step();
        let midpoint_voltage = report.resistors[0].voltage_v;
        let before_phi = body.circuit().fluxes_wb()[0];
        let after_phi = proposed.circuit().fluxes_wb()[0];
        let before_other_q = body.circuit().charges_c()[0];
        let after_other_q = proposed.circuit().charges_c()[0];
        near((after_phi-before_phi)/0.03,midpoint_voltage,2e-10);
        near((proposed.charge_c()-body.charge_c()+after_other_q-before_other_q)/0.03
            +(after_phi+before_phi)/(2.0*0.5)+midpoint_voltage/2.0,0.0,2e-9);
        near(after_other_q/0.3,proposed.voltage_v().unwrap(),2e-10);
        near(report.other_electrical_energy_j[1],
            after_phi*after_phi/(2.0*0.5)+after_other_q*after_other_q/(2.0*0.3),2e-10);
        near(report.stored_after_j+report.heat_j-report.stored_before_j,
            0.02*(proposed.extension_m()-body.extension_m()),2e-10);
        body = next.into_parts().0;
    }
}

#[test]
fn g1_native_open_and_short_limits_use_actual_coupled_initial_voltage() {
    for short in [false,true] {
        let x = 0.02;
        let q = if short { 0.5*x } else { 0.0 };
        let mut branches = vec![(1,0,Branch::Capacitor { farads: f64::NAN })];
        if short { branches.push((1,0,Branch::VoltageSource { port: 0 })); }
        let ports: &[f64] = if short { &[0.0] } else { &[] };
        let body = PiezoelectricCircuit::new(CircuitGraph {
            node_count: 2, branches, transformers: vec![],
        },vec![],vec![],vec![q],0,rod(0.2,3.0),x,0.0,ports,||false).unwrap();
        // The nonzero shorted charge is compatible with the coupled law even
        // though Q/C alone would report a nonzero voltage.
        let next = body.step_candidate(vec![],ports,0.0,0.1,solve(),||false).unwrap();
        let stiffness = if short { 2.5 } else { 3.125 };
        let a = stiffness*0.1*0.1/4.0;
        let expected = x*(1.0-a)/(1.0+a);
        near(next.model().extension_m(),expected,2e-10);
        near(next.model().charge_c(),if short { 0.5*expected } else { 0.0 },2e-10);
        near(next.step().stored_after_j,next.step().stored_before_j,2e-10);
        near(next.step().heat_j,0.0,0.0);
        near(next.step().electrical_supply_j,0.0,0.0);
    }
}

fn thermal_parts<'a>(cards: &'a [MaterialCard;2], curves: &'a [ResolvedHeatCapacityCurve;2],
    context: &'a QueryPoint, key: &'a PropertyKey)
    -> ([MaterialCell<'a>;2],[CellResistor;2],[ConductionLink;1]) {
    let cells = [0,1].map(|index| MaterialCell {
        cell: thermal_source::DESTINATIONS[index].cell, material_id: MaterialId(index as u64),
        curve: curves[index].curve(),
        point: ThermalMaterialPoint { material: &cards[index], context,
            temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((key,MaterialPropertySelection::SingleClaimOnly)),
    });
    let resistors = thermal_source::DESTINATIONS.map(|destination| CellResistor {
        destination, selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0, area_m2: 1.0,
    });
    (cells,resistors,[ConductionLink { link: 7, distance_a_m: 0.5, distance_b_m: 0.5 }])
}
fn thermal_controls() -> CircuitEnthalpySolve {
    let mut controls = thermal_source::config();
    controls.cell_energy_tolerance_j = 1e-10;
    controls.thermal_energy_tolerance_j = 1e-10;
    controls.total_energy_tolerance_j = 2e-10;
    controls.difference_step_j = 1e-5;
    controls.newton.absolute_tolerance = 1e-12;
    controls.circuit = solve();
    controls
}

#[test]
fn g1_finite_load_heating_changes_the_same_piezo_motion_and_current() {
    let run = |alpha: f64| {
        let (a,ca) = thermal_source::material(1.0,0.001,alpha);
        let (b,cb) = thermal_source::material(2.0,0.002,alpha);
        let cards = [a,b]; let curves = [ca,cb];
        let context = thermal_source::point(300.0);
        let key = PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let (cells,resistors,links) = thermal_parts(&cards,&curves,&context,&key);
        let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let mut thermal = thermal_source::initial(&curves,1.0,None);
        let prepared = material.prepare_with_bases(&thermal,&[
            (1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity),
        ],||false).unwrap();
        let mut body = PiezoelectricCircuit::new(CircuitGraph { node_count: 3, branches: vec![
            (1,0,Branch::Capacitor { farads: f64::NAN }),
            (1,2,Branch::Resistor { ohms: f64::NAN }),
            (2,0,Branch::Resistor { ohms: f64::NAN }),
        ], transformers: vec![] },prepared.resolve(&thermal,||false).unwrap().resistors,
            vec![],vec![0.1],0,rod(0.2,3.0),0.02,0.0,&[],||false).unwrap();
        let initial = thermal.energy_j()+body.piezoelectric_energy_j().unwrap();
        let mut total_work = 0.0;
        for index in 1..=6 {
            let target = f64::from(index)*0.1;
            let h = target-thermal.time_s();
            let next = heated::step_candidate(&thermal,target,PiezoelectricEnthalpyInput {
                model: &body, electrical_ports: &[], force_n: 0.01,
                destinations: &thermal_source::DESTINATIONS,
            },thermal_controls(),&[0.0;2],|state|prepared.resolve(state,||false),||false).unwrap();
            let t1 = next.thermal().temperature_k(CellId(30)).unwrap();
            let t2 = next.thermal().temperature_k(CellId(10)).unwrap();
            let r1 = 1.0+alpha*(t1-300.0);
            let r2 = 2.0*(1.0+alpha*(t2-300.0));
            let expected = load_oracle(body.extension_m(),body.velocity_m_s(),body.charge_c(),
                0.2,3.0,r1+r2,0.01,h);
            near(next.model().extension_m(),expected.0,2e-9);
            near(next.model().charge_c(),expected.2,2e-9);
            let report = next.report();
            let current = expected.3/(r1+r2);
            near(report.electrical.resistors[0].current_a,current,2e-9);
            near(report.electrical.resistors[1].current_a,current,2e-9);
            near(report.electrical.resistors[0].heat_j,h*r1*current*current,2e-10);
            near(report.electrical.resistors[1].heat_j,h*r2*current*current,2e-10);
            near(report.deposited_joule_j,h*(r1+r2)*current*current,2e-10);
            near(report.thermal.stored_energy_change_j,report.deposited_joule_j,2e-10);
            near(report.total_energy_residual_j,0.0,2e-10);
            assert_eq!(report.electrical.resistors[0].material.conductor.material().card_identity(),
                cards[0].content_hash());
            total_work += report.electrical.mechanical_supply_j;
            near(next.thermal().energy_j()+next.model().piezoelectric_energy_j().unwrap(),
                initial+total_work,2e-9);
            let (next_thermal,next_body,_) = next.into_parts();
            thermal = next_thermal; body = next_body;
        }
        assert!(thermal.temperature_k(CellId(30)).unwrap()>300.1);
        assert!(thermal.temperature_k(CellId(10)).unwrap()>300.1);
        (body.extension_m(),body.charge_c())
    };
    let fixed = run(0.0);
    let heated = run(0.005);
    assert!((heated.0-fixed.0).abs()>1e-7);
    assert!((heated.1-fixed.1).abs()>1e-6);
}

#[test]
fn g4_refused_native_and_thermal_candidates_preserve_all_accepted_histories() {
    let body = loaded(0.2,3.0,1.0,0.02,0.0,0.02);
    let before = body.clone();
    assert!(body.step_candidate(vec![],&[],0.0,0.1,solve(),||false).is_err());
    assert!(body.step_candidate(vec![source(1,1.0)],&[],1e8,0.1,solve(),||false).is_err());
    let mut count = 0;
    body.step_candidate(vec![source(1,1.0)],&[],0.01,0.1,solve(),||{
        count+=1;false
    }).unwrap();
    let mut calls = 0;
    assert!(matches!(body.step_candidate(vec![source(1,1.0)],&[],0.01,0.1,solve(),||{
        calls+=1;calls==count
    }),Err(MaterialCircuitError::Cancelled)));
    assert_eq!(body,before);
    assert!(PiezoelectricCircuit::new(CircuitGraph { node_count: 2, branches: vec![
        (1,0,Branch::Capacitor { farads: f64::NAN }),
        (1,0,Branch::VoltageSource { port: 0 }),
    ],transformers:vec![]},vec![],vec![],vec![0.0],0,rod(0.2,3.0),0.02,0.0,
        &[0.0],||false).is_err());

    let (a,ca) = thermal_source::material(1.0,0.001,0.005);
    let (b,cb) = thermal_source::material(2.0,0.002,0.005);
    let cards = [a,b]; let curves = [ca,cb];
    let context = thermal_source::point(300.0);
    let key = PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let (cells,resistors,links) = thermal_parts(&cards,&curves,&context,&key);
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = thermal_source::initial(&curves,1.0,None);
    let prepared = material.prepare_with_bases(&thermal,&[
        (1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity),
    ],||false).unwrap();
    let body = PiezoelectricCircuit::new(CircuitGraph { node_count: 3, branches: vec![
        (1,0,Branch::Capacitor { farads: f64::NAN }),
        (1,2,Branch::Resistor { ohms: f64::NAN }),
        (2,0,Branch::Resistor { ohms: f64::NAN }),
    ],transformers:vec![]},prepared.resolve(&thermal,||false).unwrap().resistors,
        vec![],vec![0.1],0,rod(0.2,3.0),0.02,0.0,&[],||false).unwrap();
    let before_body = body.clone(); let before_cells = thermal.cells().clone();
    let input = PiezoelectricEnthalpyInput {
        model: &body,electrical_ports: &[],force_n: 0.01,destinations: &thermal_source::DESTINATIONS,
    };
    let mut count = 0;
    let accepted = heated::step_candidate(&thermal,0.1,input,thermal_controls(),&[0.0;2],
        |state|prepared.resolve(state,||false),||{count+=1;false}).unwrap();
    let mut calls = 0;
    assert!(matches!(heated::step_candidate(&thermal,0.1,input,thermal_controls(),&[0.0;2],
        |state|prepared.resolve(state,||false),||{calls+=1;calls==count}),
        Err(CircuitEnthalpyError::Cancelled)));
    let mut budget = thermal_controls(); budget.maximum_circuit_trials = 1;
    assert!(heated::step_candidate(&thermal,0.1,input,budget,&[0.0;2],
        |state|prepared.resolve(state,||false),||false).is_err());
    assert_eq!(body,before_body);
    assert_eq!(thermal.cells(),&before_cells);
    assert_eq!(thermal.time_s(),0.0);
    let retry = heated::step_candidate(&thermal,0.1,input,thermal_controls(),&[0.0;2],
        |state|prepared.resolve(state,||false),||false).unwrap();
    assert_eq!(accepted.model(),retry.model());
    assert_eq!(accepted.thermal().cells(),retry.thermal().cells());
}
