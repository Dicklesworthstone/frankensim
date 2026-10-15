//! Shared source fixtures, actual circuit/gas/phase owners, independent balances.
use super::*;
use crate::material_circuit::{MaterialCircuit, MaterialCircuitSolve, MaterialResistor};
use crate::material_circuit::thermal::{CircuitHeatingError, CircuitThermalSolve, CircuitThermalSource};
use crate::cavity_phs::piston::thermal::ThermalPistonError;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_phs::circuit::{Branch, CircuitGraph};

fn electrical() -> CircuitThermalSolve {
    CircuitThermalSolve { circuit:MaterialCircuitSolve {energy_tolerance_j:1e-8,power_tolerance_w:1e-9},
        total_energy_tolerance_j:1e-6,maximum_circuit_trials:130 }
}
fn capacitor(heater:&MaterialOhmicBody<'_>,thermal:UniformEnthalpyStepInput<'_>,charge:f64) -> MaterialCircuit {
    MaterialCircuit::new(CircuitGraph {node_count:2,branches:vec![
        (1,0,Branch::Capacitor {farads:0.002}),(1,0,Branch::Resistor {ohms:1.0}),
    ],transformers:vec![]},vec![MaterialResistor {branch:1,
        conductor:heater.resolve_conductor(thermal.initial,thermal.volume_m3).unwrap()}],vec![],vec![charge],&[]).unwrap()
}
fn link(_: &ResolvedThermochemicalGasState,b:EquilibriumPhaseState)
    -> Result<ConductanceSample<EquilibriumPhaseState>,JouleHeatingError> {
    Ok(ConductanceSample {conductance_w_k:1.0,source:b})
}
fn old_piston() -> GasPiston { GasPiston::new(gas(0.75,0.001),1.0,0.01,0.0).unwrap() }

#[test]
fn source_bound_capacitor_discharge_heats_body_gas_and_motion_once() {
    let (card,caloric,key)=fixture(Some(1.0),"T"); let context=point(300.0);
    let ideal=heater(&card,&context,OhmicDrive::Current(0.0)); let ohmic=ideal.ohmic_body();
    let thermal=body(caloric.curve(),100000.0); let circuit=capacitor(&ohmic,thermal,0.2);
    let saved=circuit.clone(); let piston=old_piston();
    let layers=ideal.conduction_layers(&key,MaterialPropertySelection::SingleClaimOnly,0.02,0.001,0.02);
    let result=piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),
        CircuitThermalSource {circuit:&circuit,resistor_branch:1,body:ohmic.clone(),inputs:&[]},electrical(),
        |g,b| layers.resolve(g,b).map_err(JouleHeatingError::Thermal),||false).unwrap();
    let e=result.electrical().step();
    assert!(e.heat_j>0.0 && e.stored_after_j<10.0);
    assert!((e.stored_after_j+e.heat_j-10.0).abs()<1e-8);
    assert_eq!(e.native.supplied,0.0);
    assert!(result.coupled().thermal().state.temperature_k()>400.0);
    assert!(result.coupled().piston().gas().state().state().temperature>400.0);
    assert!(result.coupled().piston().velocity_m_s()>0.0);
    let midpoint=result.coupled().step().thermal_midpoint;
    assert_eq!(e.resistors[0].material.conductor,ohmic.resolve_conductor(midpoint,thermal.volume_m3).unwrap());
    assert_eq!(circuit,saved);
    assert!(result.total_energy_residual_j().abs()<electrical().total_energy_tolerance_j);
    let (_,_,next,report)=result.into_parts();
    assert!(report.circuit_trials>1); // Many trial solves, ONE accepted debit.
    assert!(next.charges_c()[0]<saved.charges_c()[0]);
    assert!((report.power.evaluated_energy_j-report.electrical.heat_j).abs()<1e-9);
}

fn trajectory(scale:f64) -> (GasPiston,MaterialCircuit,f64) {
    let (card,caloric,_)=fixture(Some(scale),"T"); let context=point(300.0);
    let ohmic=heater(&card,&context,OhmicDrive::Voltage(0.0)).ohmic_body();
    let mut thermal=body(caloric.curve(),100000.0); let h0=thermal.initial.specific_enthalpy_j_kg();
    let mut circuit=capacitor(&ohmic,thermal,0.2);
    let mut piston=old_piston(); let initial=piston.clone(); let mut heat=0.0;
    for _ in 0..80 {
        let result=piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),
            CircuitThermalSource {circuit:&circuit,resistor_branch:1,body:ohmic.clone(),inputs:&[]},
            electrical(),link,||false).unwrap();
        let (next,accepted,storage,report)=result.into_parts();
        heat+=report.electrical.heat_j;
        assert!(report.total_energy_residual_j.abs()<electrical().total_energy_tolerance_j);
        assert_eq!(accepted.state.material_card_identity(),card.content_hash());
        assert!((report.electrical.stored_after_j+heat-10.0).abs()<80.0*electrical().circuit.energy_tolerance_j);
        let delta_gas=next.gas().mass_kg()*(next.gas().state().specific_internal_energy_j_kg()
            -initial.gas().state().specific_internal_energy_j_kg());
        let kinetic=0.5*next.mass_kg()*next.velocity_m_s().powi(2);
        let wall=thermal.mass_kg*(accepted.state.specific_enthalpy_j_kg()-h0);
        let exterior=100000.0*(next.gas().volume_m3()-initial.gas().volume_m3());
        assert!((report.electrical.stored_after_j-10.0+delta_gas+kinetic+wall+exterior).abs()<1e-5);
        piston=next; circuit=storage; thermal.initial=accepted.state;
    }
    (piston,circuit,heat)
}
#[test]
fn resistance_only_data_change_alters_finite_supply_depletion_and_the_real_trajectory() {
    let (a,ca,qa)=trajectory(1.0); let (b,cb,qb)=trajectory(2.0);
    assert!(qa>qb && qb>0.0 && qa<10.0);
    assert!(ca.charges_c()[0]<cb.charges_c()[0]);
    assert!((a.position_m()-b.position_m()).abs()>1e-8);
    assert!((a.gas().state().state().sound_speed-b.gas().state().state().sound_speed).abs()>1e-4);
}

#[test]
fn latent_phase_changes_resistance_and_charge_drain_without_temperature_animation() {
    let (card,_,_)=fixture(Some(1.0),"melt"); let context=point(300.0);
    let knot=|h,t,x| EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
        liquid_mass_fraction:x,bulk_density_kg_m3:4000.0};
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),vec![knot(0.0,200.0,0.0),
        knot(100000.0,400.0,0.0),knot(160000.0,400.0,1.0),knot(400000.0,800.0,1.0)]).unwrap();
    let thermal=body(&curve,130000.0);
    let mut ohmic=heater(&card,&context,OhmicDrive::Current(0.0)).ohmic_body();
    ohmic.material.liquid_fraction_axis=Some("melt");
    let circuit=capacitor(&ohmic,thermal,0.2); let piston=old_piston();
    let result=piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),
        CircuitThermalSource {circuit:&circuit,resistor_branch:1,body:ohmic,inputs:&[]},electrical(),
        |_,b| Ok::<_,JouleHeatingError>(ConductanceSample {conductance_w_k:0.0,source:b}),||false).unwrap();
    let heat=result.electrical().step().heat_j;
    let r=15.0+10.0*heat/(2.0*thermal.mass_kg*60000.0);
    let tau2=2.0*0.002*r; let dt=thermal.duration_s;
    // Independently derived implicit RC/latent equation, not another heat step.
    assert!((heat*(tau2+dt).powi(2)-10.0*8.0*0.002*r*dt).abs()<1e-10);
    assert_eq!(result.coupled().thermal().state.temperature_k(),400.0);
    assert!(result.coupled().thermal().state.liquid_mass_fraction()>0.5);
    assert!(result.electrical().step().resistors[0].material.conductor.resistance_ohm()>15.0);
    assert_eq!(result.coupled().piston(),&piston);
}

#[test]
fn zero_charge_preserves_all_history_but_does_not_hide_source_support() {
    let (card,caloric,_)=fixture(Some(1.0),"T"); let context=point(300.0);
    let thermal=body(caloric.curve(),100000.0);
    let ohmic=heater(&card,&context,OhmicDrive::Voltage(0.0)).ohmic_body();
    let circuit=capacitor(&ohmic,thermal,0.0); let piston=old_piston();
    let result=piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),
        CircuitThermalSource {circuit:&circuit,resistor_branch:1,body:ohmic.clone(),inputs:&[]},electrical(),link,||false).unwrap();
    assert_eq!(result.electrical().step().heat_j,0.0);
    assert_eq!(result.coupled().piston(),&piston);
    assert_eq!(result.coupled().thermal().state,thermal.initial);
    let outside=context.clone().with("lot",9.0).unwrap();
    let mut bad=ohmic; bad.material.context=&outside;
    assert!(piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),
        CircuitThermalSource {circuit:&circuit,resistor_branch:1,body:bad,inputs:&[]},electrical(),link,||false).is_err());
}

#[test]
fn circuit_budget_and_final_cancellation_never_partially_drain_the_store() {
    let (card,caloric,_)=fixture(Some(1.0),"T"); let context=point(300.0);
    let thermal=body(caloric.curve(),100000.0);
    let ohmic=heater(&card,&context,OhmicDrive::Voltage(0.0)).ohmic_body();
    let circuit=capacitor(&ohmic,thermal,0.2); let saved=circuit.clone();
    let piston=old_piston(); let saved_piston=piston.clone();
    let source=|| CircuitThermalSource {circuit:&circuit,resistor_branch:1,body:ohmic.clone(),inputs:&[]};
    assert!(matches!(piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),source(),
        CircuitThermalSolve {maximum_circuit_trials:1,..electrical()},link,||false),
        Err(ThermalPistonError::Interface(CircuitHeatingError::Budget))));
    let mut total=0;
    piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),source(),electrical(),link,
        || {total+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),source(),electrical(),link,
        || {actual+=1;actual==total}),Err(ThermalPistonError::Cancelled)));
    assert_eq!(circuit,saved); assert_eq!(piston,saved_piston);
}

#[test]
fn unbound_resistive_sinks_cannot_be_counted_as_whole_body_heat() {
    let (card,caloric,_)=fixture(Some(1.0),"T"); let context=point(300.0);
    let thermal=body(caloric.curve(),100000.0);
    let ohmic=heater(&card,&context,OhmicDrive::Voltage(0.0)).ohmic_body();
    let resistor=ohmic.resolve_conductor(thermal.initial,thermal.volume_m3).unwrap();
    let circuit=MaterialCircuit::new(CircuitGraph {node_count:3,branches:vec![
        (1,0,Branch::Capacitor {farads:0.002}),(1,2,Branch::Resistor {ohms:1.0}),
        (2,0,Branch::Resistor {ohms:1.0}),
    ],transformers:vec![]},vec![MaterialResistor {branch:1,conductor:resistor.clone()},
        MaterialResistor {branch:2,conductor:resistor}],vec![],vec![0.2],&[]).unwrap();
    assert!(matches!(old_piston().circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),
        CircuitThermalSource {circuit:&circuit,resistor_branch:1,body:ohmic,inputs:&[]},electrical(),link,||false),
        Err(ThermalPistonError::Invalid(_))));
}
