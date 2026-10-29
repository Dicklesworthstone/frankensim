//! Shared source fixtures, actual circuit/gas/phase owners, independent balances.
use super::*;
use crate::material_circuit::{MaterialCircuit, MaterialCircuitSolve, MaterialResistor};
use crate::material_circuit::thermal::{CircuitHeatingError, CircuitThermalSolve, CircuitThermalSource};
use crate::cavity_phs::piston::thermal::ThermalPistonError;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_phs::circuit::{Branch, CircuitGraph};
use crate::material_circuit::conductivity::{ConductionAxis, ConductionCoordinate,
    ConductivitySolve};

fn electrical() -> CircuitThermalSolve {
    CircuitThermalSolve { circuit:MaterialCircuitSolve {energy_tolerance_j:1e-8,power_tolerance_w:1e-9},
        total_energy_tolerance_j:1e-6,maximum_circuit_trials:130 }
}
fn capacitor(heater:&MaterialOhmicBody<'_>,thermal:UniformEnthalpyStepInput<'_>,charge:f64) -> MaterialCircuit {
    capacitor_with_basis(heater,thermal,charge,ElectricalPropertyBasis::Resistivity)
}
fn capacitor_with_basis(heater:&MaterialOhmicBody<'_>,thermal:UniformEnthalpyStepInput<'_>,charge:f64,
    basis:ElectricalPropertyBasis) -> MaterialCircuit {
    MaterialCircuit::new(CircuitGraph {node_count:2,branches:vec![
        (1,0,Branch::Capacitor {farads:0.002}),(1,0,Branch::Resistor {ohms:1.0}),
    ],transformers:vec![]},vec![MaterialResistor {branch:1,
        conductor:heater.resolve_conductor_with_basis(thermal.initial,thermal.volume_m3,basis).unwrap()}],vec![],vec![charge],&[]).unwrap()
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
    // One step from rest the initial pressure imbalance sets the sign of the
    // velocity; the discharge must add outward motion to that same step.
    let cold=capacitor(&ohmic,thermal,0.0);
    let unheated=piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),
        CircuitThermalSource {circuit:&cold,resistor_branch:1,body:ohmic.clone(),inputs:&[]},electrical(),
        |g,b| layers.resolve(g,b).map_err(JouleHeatingError::Thermal),||false).unwrap();
    assert!(result.coupled().piston().velocity_m_s()>unheated.coupled().piston().velocity_m_s(),
        "{} vs {}",result.coupled().piston().velocity_m_s(),unheated.coupled().piston().velocity_m_s());
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

#[test]
fn conductivity_only_card_drives_actual_circuit_body_gas_and_piston_transaction() {
    let basis=ElectricalPropertyBasis::Conductivity;
    // This authored source is sigma(T)=0.1*T S/m, independently evaluated
    // before deriving rho=1/sigma. The card supplies no resistivity claim.
    let (card,caloric,key)=fixture_with_basis(Some(1600.0),"T",basis);let context=point(300.0);
    assert!(card.claims_for("electrical_resistivity").is_empty());
    let ideal=heater(&card,&context,OhmicDrive::Voltage(50.0));let ohmic=ideal.ohmic_body();
    let thermal=body(caloric.curve(),100000.0);let piston=old_piston();
    assert!(ideal.resolve(thermal.initial,thermal.volume_m3).is_err());
    let power=ideal.resolve_with_basis(thermal.initial,thermal.volume_m3,basis).unwrap();
    assert!((power.power_w-250.0).abs()<1e-10);
    let circuit=capacitor_with_basis(&ohmic,thermal,0.2,basis);let saved=circuit.clone();
    let layers=ideal.conduction_layers(&key,MaterialPropertySelection::SingleClaimOnly,0.02,0.001,0.02);
    let supply=||CircuitThermalSource {circuit:&circuit,resistor_branch:1,body:ohmic.clone(),inputs:&[]};
    assert!(piston.circuit_thermal_candidate(drive(0.0),thermal,controls(),bounds(),bounds(),supply(),electrical(),
        |g,b|layers.resolve(g,b).map_err(JouleHeatingError::Thermal),||false).is_err());
    let next=piston.circuit_thermal_candidate_with_basis(drive(0.0),thermal,controls(),bounds(),bounds(),
        supply(),electrical(),basis,|g,b|layers.resolve(g,b).map_err(JouleHeatingError::Thermal),||false).unwrap();
    let midpoint=next.coupled().step().thermal_midpoint;
    let electrical=next.electrical().step();let conductor=&electrical.resistors[0].material.conductor;
    assert_eq!(conductor.source_basis(),basis);
    assert_eq!(conductor, &ohmic.resolve_conductor_with_basis(midpoint,thermal.volume_m3,basis).unwrap());
    assert_eq!(conductor.source_property().answer().receipt.selected,card.claims_for(basis.property_name())[0].0);
    assert!((conductor.resistance_ohm()-4000.0/midpoint.temperature_k()).abs()<1e-10);
    assert!(electrical.heat_j>0.0 && electrical.stored_after_j<10.0);
    assert!(next.coupled().thermal().state.temperature_k()>400.0);
    assert!(next.coupled().piston().gas().state().state().temperature>400.0);
    assert!((electrical.stored_after_j+electrical.heat_j-10.0).abs()<1e-8);
    assert!(next.total_energy_residual_j().abs()<1e-6);
    assert_eq!(circuit,saved);
}

const FIELD_DIMS: Dims = Dims([1, 1, -3, 0, -1, 0]);
const FIELD_AXES: [ConductionAxis<'static>; 1] = [ConductionAxis {
    name: "field", coordinate: ConductionCoordinate::ElectricFieldMagnitude,
}];

fn field_point(t: f64) -> QueryPoint {
    point(t).with_quantity("field", QuantitySpec::dimensional(FIELD_DIMS), 0.0).unwrap()
}

fn field_fixture(slope: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (cp.clone(), PropertyValue::Scalar { value: 500.0, dims: SPECIFIC_HEAT_CAPACITY_DIMS }),
        (density.clone(), PropertyValue::Scalar { value: 4000.0, dims: Density::DIMS }),
        (PropertyKey::new(basis.property_name(), basis.dimensions()), PropertyValue::Surface2D {
            axes: ["T".into(), "field".into()], axis_dims: [Temperature::DIMS, FIELD_DIMS],
            samples: [[200.0, 0.0], [800.0, 0.0], [200.0, 20000.0], [800.0, 20000.0]]
                .into_iter().map(|at| (at, 0.00005 * at[0] + slope * at[1])).collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]], dims: basis.dimensions(),
        }),
    ] {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else { InterpolationPolicy::LinearInside };
        claims.insert_claim(PropertyClaim { key, value, interpolation,
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 200.0, 800.0)
                .with_quantity("field", QuantitySpec::dimensional(FIELD_DIMS), 0.0, 20000.0),
            uncertainty: UncertaintyModel::Unstated, observations: vec![],
            provenance: Provenance { source: "manufactured rho(T,E), not a measured material".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "authored nonlinear resistor".into(),
        phase: "solid".into(), process: "analytic thermal/field fixture".into(), revision: 0 },
        claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &density, temperature_axis: "T",
        lower: &field_point(200.0), upper: &field_point(800.0),
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-5).unwrap();
    (card, caloric)
}

fn field_capacitor(heater: &MaterialOhmicBody<'_>, thermal: UniformEnthalpyStepInput<'_>) -> MaterialCircuit {
    let conductor = resolve_uniform_conductor(heater.material.material,
        &heater.material.at(thermal.initial).unwrap(), MaterialPropertySelection::SingleClaimOnly,
        heater.length_m, heater.area_m2).unwrap();
    MaterialCircuit::new(CircuitGraph { node_count: 2, branches: vec![
        (1, 0, Branch::Capacitor { farads: 0.002 }), (1, 0, Branch::Resistor { ohms: 1.0 }),
    ], transformers: vec![] }, vec![MaterialResistor { branch: 1, conductor }],
        vec![], vec![0.2], &[]).unwrap()
}

fn nonlinear() -> ConductivitySolve {
    ConductivitySolve { maximum_iterations: 32, absolute_resistance_tolerance_ohm: 1e-10,
        relative_resistance_tolerance: 1e-12 }
}

fn nonlinear_electrical() -> CircuitThermalSolve {
    CircuitThermalSolve { maximum_circuit_trials: 4096, ..electrical() }
}

#[test]
fn operating_field_and_temperature_close_inside_the_native_gas_piston_transaction() {
    let mut results = Vec::new();
    for slope in [0.0, 1e-6] {
        let (card, caloric) = field_fixture(slope);
        let context = field_point(300.0);
        let ohmic = heater(&card, &context, OhmicDrive::Current(0.0)).ohmic_body();
        let thermal = body(caloric.curve(), 100000.0);
        let circuit = field_capacitor(&ohmic, thermal);
        let saved = circuit.clone();
        let piston = old_piston();
        let supply = || CircuitThermalSource { circuit: &circuit, resistor_branch: 1,
            body: ohmic.clone(), inputs: &[] };
        assert!(piston.circuit_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
            supply(), electrical(), link, || false).is_err());
        let candidate = piston.circuit_conductivity_thermal_candidate(
            drive(0.0), thermal, controls(), bounds(), bounds(), supply(), nonlinear_electrical(),
            &FIELD_AXES, nonlinear(), link, || false,
        ).unwrap();
        let midpoint = candidate.coupled().step().thermal_midpoint.temperature_k();
        let electrical = candidate.electrical().step();
        let sample = &electrical.evaluated[0];
        let native = &electrical.electrical;
        let field = sample.point.electric_field_v_m;
        let resistance = 400.0 * (0.00005 * midpoint + slope * field.abs());
        assert!((sample.conductor.resistance_ohm() - resistance).abs() < 1e-10);
        // Independent RC midpoint identities, with the SAME solved R(T,E).
        let current = 0.2 / (0.002 * resistance + 0.5 * thermal.duration_s);
        assert!((sample.point.current_a - current).abs() < 1e-8);
        assert!((sample.point.voltage_v / ohmic.length_m - field).abs() < 1e-8);
        assert!((native.heat_j - thermal.duration_s * current.powi(2) * resistance).abs() < 1e-9);
        assert!((native.stored_after_j + native.heat_j - 10.0).abs() < 1e-8);
        assert!(sample.power_defect_w.abs() <= nonlinear_electrical().circuit.power_tolerance_w);
        let point = sample.conductor.material().query_point();
        assert_eq!(point.iter().find(|(key, _)| key == "T").map(|(_, value)| *value), Some(midpoint));
        assert_eq!(point.iter().find(|(key, _)| key == "field").map(|(_, value)| *value), Some(field.abs()));
        assert_eq!(circuit, saved);
        let (next_piston, next_body, next_circuit, report) = candidate.into_parts();
        assert!(report.total_energy_residual_j.abs() < 1e-6);
        assert!(report.circuit_trials >= report.electrical.iterations);
        assert!((report.power.evaluated_energy_j - report.electrical.electrical.heat_j).abs() < 1e-12);
        results.push((report.electrical.electrical.heat_j, next_body.state.temperature_k(),
            next_piston.gas().state().state().temperature, next_piston.velocity_m_s(),
            next_circuit.charges_c()[0]));
    }
    // A material-only operating-field response changes the ACTUAL current,
    // finite supply drain, body/gas temperatures and pressure-driven motion.
    assert!(results[0].0 > results[1].0);
    assert!(results[0].1 > results[1].1);
    assert!(results[0].2 > results[1].2);
    assert!(results[0].3 > results[1].3);
    assert!(results[0].4 < results[1].4);
}

#[test]
fn nonlinear_thermal_work_failure_and_late_cancellation_preserve_all_accepted_owners() {
    let (card, caloric) = field_fixture(1e-6);
    let context = field_point(300.0);
    let ohmic = heater(&card, &context, OhmicDrive::Current(0.0)).ohmic_body();
    let thermal = body(caloric.curve(), 100000.0);
    let circuit = field_capacitor(&ohmic, thermal);
    let saved = circuit.clone();
    let piston = old_piston();
    let saved_piston = piston.clone();
    let supply = || CircuitThermalSource { circuit: &circuit, resistor_branch: 1,
        body: ohmic.clone(), inputs: &[] };
    let run = |work, solve, axes: &[ConductionAxis<'_>], cancel: &mut dyn FnMut() -> bool| {
        piston.circuit_conductivity_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
            supply(), work, axes, solve, link, cancel)
    };
    assert!(matches!(run(CircuitThermalSolve { maximum_circuit_trials: 1, ..nonlinear_electrical() },
        nonlinear(), &FIELD_AXES, &mut || false),
        Err(ThermalPistonError::Interface(CircuitHeatingError::Budget))));
    assert!(matches!(run(nonlinear_electrical(), ConductivitySolve { maximum_iterations: 1, ..nonlinear() },
        &FIELD_AXES, &mut || false),
        Err(ThermalPistonError::Interface(CircuitHeatingError::Circuit(
            crate::material_circuit::MaterialCircuitError::ConductivityConvergence { .. })))));
    assert!(run(nonlinear_electrical(), nonlinear(), &[], &mut || false).is_err());
    let outside = context.clone().with("field", 30000.0).unwrap();
    let bad = CircuitThermalSource { circuit: &circuit, resistor_branch: 1,
        body: MaterialOhmicBody { material: ThermalMaterialPoint { context: &outside, ..ohmic.material },
            ..ohmic.clone() }, inputs: &[] };
    assert!(piston.circuit_conductivity_thermal_candidate(drive(0.0), thermal, controls(), bounds(), bounds(),
        bad, nonlinear_electrical(), &FIELD_AXES, nonlinear(), link, || false).is_err());
    let mut count = 0;
    run(nonlinear_electrical(), nonlinear(), &FIELD_AXES, &mut || { count += 1; false }).unwrap();
    let mut calls = 0;
    assert!(matches!(run(nonlinear_electrical(), nonlinear(), &FIELD_AXES,
        &mut || { calls += 1; calls == count }), Err(ThermalPistonError::Cancelled)));
    assert_eq!(circuit, saved);
    assert_eq!(piston, saved_piston);
}
