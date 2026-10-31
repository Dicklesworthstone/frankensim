//! G1/G3/G4: original field dynamics, full source bundles and finite receivers.
use super::*;
#[path = "../electrical_contacts/source.rs"]
mod scalar;
use fs_couple::circuit_enthalpy::motion::{FieldMotionEnthalpyInput, conductivity};
use fs_couple::circuit_enthalpy::contacts::{ContactHeatCells, ContactThermalCoefficients};
use fs_couple::material_circuit::contacts::SeriesElectricalContact;
use fs_couple::material_circuit::conductivity::contact::{ContactAxis, ContactCoordinate,
    ContactDirection, ContactResponse};
use fs_couple::material_circuit::motion::MaterialFieldMotion;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_ga::Vec3;
use fs_material::field_storage::{ResolvedCapacitor, ELECTRIC_PERMITTIVITY_DIMS,
    ELECTRIC_PERMITTIVITY_PROPERTY};
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};

fn moving(resistors: Vec<MaterialResistor>, driven: bool) -> MaterialFieldMotion {
    let dielectric = scalar::scalar_card(ELECTRIC_PERMITTIVITY_PROPERTY,ELECTRIC_PERMITTIVITY_DIMS,1.0);
    let capacitor = ResolvedCapacitor::resolve(&dielectric,&QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,1.0,1.0).unwrap();
    let material = scalar::scalar_card("density",Density::DIMS,100.0);
    let body = resolve_body(&[RegionBinding {key:"guided electrode",material:&material,
        point:&QueryPoint::new(),selection:MaterialPropertySelection::SingleClaimOnly,
        geometry:VolumeMoments::cuboid(Vec3::new(1.0,1.0,1.0)).unwrap()}],1,||false).unwrap();
    let branches = if driven {
        vec![(1,0,Branch::Capacitor {farads:f64::NAN}),
            (2,0,Branch::Resistor {ohms:f64::NAN}),
            (2,0,Branch::Resistor {ohms:f64::NAN}),
            (2,0,Branch::VoltageSource {port:0})]
    } else {
        vec![(1,0,Branch::Capacitor {farads:f64::NAN}),
            (1,2,Branch::Resistor {ohms:f64::NAN}),
            (2,0,Branch::Resistor {ohms:f64::NAN})]
    };
    let graph = CircuitGraph {node_count:3,branches,transformers:vec![]};
    let charge = if driven {0.8} else {10.0};
    let ports: &[f64] = if driven {&[3.0]} else {&[]};
    let circuit = MaterialCircuit::new_with_storage(graph,resistors,
        vec![MaterialStorage {branch:0,element:MaterialStorageElement::Capacitor(capacitor)}],
        vec![],vec![charge],ports).unwrap();
    MaterialFieldMotion::new(circuit,0,body,if driven {0.02} else {0.0},[0.5,1.5]).unwrap()
}
fn stored(motion: &MaterialFieldMotion) -> f64 {
    let q = motion.circuit().charges_c()[0];
    0.5*q*q*motion.gap_m()+0.5*100.0*motion.velocity_m_s().powi(2)
}

#[test]
fn independent_constant_charge_motion_and_nonlinear_heating_share_external_work() {
    let beta = 0.05;
    let (a,ca) = field_material(beta,2000.0);
    let (b,cb) = common::material(3.0,2.0,0.005);
    let cards = [a,b]; let curves = [ca,cb]; let thermal = initial(&curves,0.0,None);
    let model = moving(circuit(&cards).resistors().to_vec(),true);
    let dt = 0.2; let force = 0.1;
    let input = FieldMotionEnthalpyInput {motion:&model,ports:&[3.0],
        external_force_n:force,destinations:&DESTINATIONS};
    let result = conductivity::step_candidate(&thermal,dt,input,&responses(&cards),config(),
        conductivity_controls(),&[0.0;2],|s|coefficients(s,model.circuit(),&cards,&curves),||false).unwrap();
    let heat_a = dt*9.0/(1.0+beta*3.0);
    let rhs = dt*9.0/(2.0*3.0);
    let rise_b = 2.0*rhs/(1.0+(1.0+4.0*0.005*rhs).sqrt());
    let acceleration = (force-0.5*0.8*0.8)/100.0;
    let displacement = 0.02*dt+0.5*acceleration*dt*dt;
    near(result.thermal().temperature_k(CellId(30)).unwrap(),300.0+heat_a);
    near(result.thermal().temperature_k(CellId(10)).unwrap(),300.0+rise_b);
    near(result.motion().gap_m(),1.0+displacement);
    near(result.motion().velocity_m_s(),0.02+acceleration*dt);
    near(result.motion().circuit().charges_c()[0],0.8);
    near(result.report().electrical.evaluated[0].point.electric_field_v_m,3.0);
    near(result.report().deposited_joule_j,heat_a+2.0*rise_b);
    let field = &result.report().electrical.field;
    near(field.electrical_supply_j,heat_a+2.0*rise_b);
    near(field.mechanical_supply_j,force*displacement);
    near(result.thermal().energy_j()-thermal.energy_j()+stored(result.motion())-stored(&model),
        field.electrical_supply_j+force*displacement);
}

#[test]
fn prepared_material_swaps_change_heat_transport_charge_and_actual_gap_without_callbacks() {
    let mut endpoints = Vec::new();
    for beta in [0.02,0.08] {
        let (a,ca) = thermo_field_material(beta); let (b,cb) = common::material(3.0,2.0,0.005);
        let cards = [a,b]; let curves = [ca,cb]; let mut thermal = initial(&curves,1.0,None);
        let context = [query(300.0),point(300.0)];
        let key = PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let cells = network_cells(&cards,&curves,&context,&key); let resistors = network_resistors();
        let network = MaterialNetwork {cells:&cells,resistors:&resistors,links:&NETWORK_LINKS};
        let responses = network_responses();
        let prepared = network.prepare_conductivity(&thermal,&responses,||false).unwrap();
        let seeds = prepared.resolve_conductivity_seeds(&thermal,&responses,||false).unwrap();
        let mut motion = moving(seeds.resistors,false); let total = thermal.energy_j()+stored(&motion);
        for interval in 1..=3 {
            let target = f64::from(interval)*0.1; let dt = target-thermal.time_s();
            let result = prepared.field_motion_conductivity_candidate(&thermal,target,&motion,&[],0.0,
                &responses,config(),conductivity_controls(),&[0.0;2],||false).unwrap();
            let ta = result.thermal().temperature_k(CellId(30)).unwrap();
            let tb = result.thermal().temperature_k(CellId(10)).unwrap();
            let electrical = &result.report().electrical;
            let operating = electrical.evaluated[0].point.electric_field_v_m.abs();
            let resistance = 1.0+0.004*(ta-300.0)+beta*operating;
            near(electrical.evaluated[0].conductor.resistance_ohm(),resistance);
            let k = 1.0+0.002*(ta-300.0)+0.3*operating;
            let conductance = 1.0/(0.5/k+0.5);
            near(result.report().source.links[0].conductance_w_k,conductance);
            let heat_a = result.thermal().cells()[&CellId(30)].specific_h_j_kg
                -thermal.cells()[&CellId(30)].specific_h_j_kg;
            near(heat_a,electrical.field.resistors[0].heat_j+dt*conductance*(tb-ta));
            near(result.report().deposited_joule_j,electrical.field.heat_j);
            near(result.thermal().energy_j()+stored(result.motion()),total);
            assert!(result.report().circuit_trials>electrical.iterations);
            let (t,m,_) = result.into_parts(); thermal=t; motion=m;
        }
        endpoints.push((motion.gap_m(),motion.circuit().charges_c()[0],
            thermal.temperature_k(CellId(30)).unwrap()));
    }
    assert!(endpoints[1].1>endpoints[0].1,"more field-dependent resistance retains charge");
    assert!(endpoints[1].0<endpoints[0].0,"retained charge strengthens actual attraction");
    assert!((endpoints[1].2-endpoints[0].2).abs()>0.01);
}

#[test]
fn all_native_trials_count_and_refusal_preserves_thermal_and_moving_histories() {
    let (a,ca) = field_material(0.05,2000.0); let (b,cb) = common::material(3.0,2.0,0.005);
    let cards = [a,b]; let curves = [ca,cb]; let thermal = initial(&curves,1.0,None);
    let model = moving(circuit(&cards).resistors().to_vec(),false);
    let input = FieldMotionEnthalpyInput {motion:&model,ports:&[],external_force_n:0.0,destinations:&DESTINATIONS};
    let accepted_model = model.clone(); let accepted_cells = thermal.cells().clone();
    let mut polls = 0;
    let result = conductivity::step_candidate(&thermal,0.1,input,&responses(&cards),config(),
        conductivity_controls(),&[0.0;2],|s|coefficients(s,model.circuit(),&cards,&curves),
        ||{polls+=1;false}).unwrap();
    let mut actual = 0;
    assert!(matches!(conductivity::step_candidate(&thermal,0.1,input,&responses(&cards),config(),
        conductivity_controls(),&[0.0;2],|s|coefficients(s,model.circuit(),&cards,&curves),
        ||{actual+=1;actual==polls}),Err(CircuitEnthalpyError::Cancelled)));
    let limit = fs_couple::circuit_enthalpy::CircuitEnthalpySolve {
        maximum_circuit_trials:result.report().circuit_trials-1,..config() };
    assert!(matches!(conductivity::step_candidate(&thermal,0.1,input,&responses(&cards),limit,
        conductivity_controls(),&[0.0;2],|s|coefficients(s,model.circuit(),&cards,&curves),||false),
        Err(CircuitEnthalpyError::CircuitBudget)));
    assert_eq!(model,accepted_model); assert_eq!(thermal.cells(),&accepted_cells); assert_eq!(thermal.time_s(),0.0);
}

const CONTACT_V: Dims = Dims([2,1,-3,0,-1,0]);
const CONTACT_ROUTES: [ContactHeatCells;1] = [ContactHeatCells {id:7,
    cells:[CellId(30),CellId(10)],heat_fraction_to_first:0.25}];

fn contact_motion(cards: &[MaterialCard;2], interface: &InterfaceSystemCard) -> MaterialFieldMotion {
    let seed = moving(circuit(cards).resistors().to_vec(),false);
    let circuit = seed.circuit().rebind_contacts_candidate(vec![
        scalar::response_contact(7,1,interface,1.0,"V",CONTACT_V)],&[]).unwrap();
    MaterialFieldMotion::new(circuit,0,seed.body().clone(),0.0,[0.5,1.5]).unwrap()
}
fn contact_coefficients(state: &DepositionSimulation<'_>, motion: &MaterialFieldMotion,
    cards: &[MaterialCard;2], curves: &[ResolvedHeatCapacityCurve;2], interface: &InterfaceSystemCard)
    -> Result<ContactThermalCoefficients<Vec<EquilibriumPhaseState>>,String>
{
    let bulk = coefficients(state,motion.circuit(),cards,curves)?;
    let old = &motion.circuit().contacts()[0];
    let point = old.source().resolved().query()
        .with("Ta",bulk.source[0].temperature_k()).map_err(|e|e.to_string())?
        .with("Tb",bulk.source[1].temperature_k()).map_err(|e|e.to_string())?;
    let contact = SeriesElectricalContact::new(old.id(),old.branch(),
        old.source().at_query_point(interface,&point).map_err(|e|e.to_string())?,old.area_m2())
        .map_err(|e|e.to_string())?;
    Ok(ContactThermalCoefficients {bulk,contacts:vec![contact]})
}

// Independent continuous five-state reference: capacitor charge, physical gap,
// velocity, and two finite temperatures. No material or coupled solver is called.
// R_a=1+0.05 V_a; R_c=1.2[1+0.001(T_a-300)+delta V_c]; R_b=3.
fn contact_reference(delta: f64) -> [f64;5] {
    let derivative = |x: [f64;5]| {
        let laws = |i: f64| (1.0/(1.0-0.05*i),
            1.2*(1.0+0.001*(x[3]-300.0))/(1.0-1.2*delta*i));
        let (mut lo,mut hi) = (0.0,4.0);
        for _ in 0..70 {
            let i=0.5*(lo+hi); let (ra,rc)=laws(i);
            if i*(ra+3.0+rc)>x[0]*x[1] {hi=i;} else {lo=i;}
        }
        let i=0.5*(lo+hi); let (ra,rc)=laws(i);
        let heat_a=i*i*ra; let heat_b=i*i*3.0; let contact=i*i*rc;
        let conductance=1.0+0.05*i*rc;
        let exchange=conductance*(x[4]-x[3]);
        [-i,x[2],-0.5*x[0]*x[0]/100.0,
            heat_a+0.25*contact+exchange,(heat_b+0.75*contact-exchange)/2.0]
    };
    let add = |x: [f64;5],k: [f64;5],h: f64| std::array::from_fn(|i|x[i]+h*k[i]);
    let mut x=[10.0,1.0,0.0,300.0,300.0]; let h=0.0001;
    for _ in 0..2000 {
        let a=derivative(x); let b=derivative(add(x,a,h*0.5));
        let c=derivative(add(x,b,h*0.5)); let d=derivative(add(x,c,h));
        x=std::array::from_fn(|i|x[i]+h*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0);
    }
    x
}
fn contact_trajectory(delta: f64, steps: usize) -> [f64;5] {
    let (a,ca)=field_material(0.05,2000.0); let (b,cb)=common::material(3.0,2.0,0.0);
    let cards=[a,b]; let curves=[ca,cb]; let mut thermal=initial(&curves,1.0,None);
    let interface=scalar::response_card(&cards[0],&cards[1],"V",CONTACT_V,1.2,0.001,delta,10.0);
    let mut motion=contact_motion(&cards,&interface);
    let axes=[ContactAxis {name:"V",coordinate:ContactCoordinate::VoltageMagnitude}];
    let contacts=[ContactResponse {id:7,card:&interface,direction:ContactDirection::AlongBranch,axes:&axes}];
    let total=thermal.energy_j()+stored(&motion); let dt=0.2/steps as f64;
    for n in 1..=steps {
        let input=FieldMotionEnthalpyInput {motion:&motion,ports:&[],external_force_n:0.0,destinations:&DESTINATIONS};
        let result=conductivity::step_with_contact_response_candidate(&thermal,n as f64*dt,input,
            &responses(&cards),Some(&contacts),&CONTACT_ROUTES,config(),conductivity_controls(),&[0.0;2],
            |state|contact_coefficients(state,&motion,&cards,&curves,&interface),
            |_,field,mut transport| {
                let voltage=field.evaluated_contacts[0].point.voltage_v.abs();
                transport.conductances_w_k.insert(7,1.0+0.05*voltage);
                Ok(transport)
            },||false).unwrap();
        let ta=result.thermal().temperature_k(CellId(30)).unwrap();
        let tb=result.thermal().temperature_k(CellId(10)).unwrap();
        let report=result.report(); let electrical=&report.electrical; let field=&electrical.field;
        let source=&electrical.evaluated_contacts[0]; let contact=&field.contacts[0];
        near(source.contact.resistance_ohm(),1.2*(1.0+0.001*(ta-300.0)+delta*contact.voltage_v.abs()));
        near(source.contact.source().resolved().state_coordinate("Ta").unwrap(),ta);
        near(source.contact.source().resolved().state_coordinate("Tb").unwrap(),tb);
        near(source.contact.source().resolved().state_coordinate("V").unwrap(),contact.voltage_v.abs());
        near(electrical.evaluated[0].conductor.resistance_ohm(),1.0+0.05*field.resistors[0].voltage_v.abs());
        let exchange=dt*(1.0+0.05*contact.voltage_v.abs())*(tb-ta);
        near(result.thermal().cells()[&CellId(30)].specific_h_j_kg-thermal.cells()[&CellId(30)].specific_h_j_kg,
            field.resistors[0].heat_j+0.25*contact.heat_j+exchange);
        near(result.thermal().cells()[&CellId(10)].specific_h_j_kg-thermal.cells()[&CellId(10)].specific_h_j_kg,
            field.resistors[1].heat_j+0.75*contact.heat_j-exchange);
        near(result.thermal().energy_j()+stored(result.motion()),total);
        near(report.deposited_joule_j,field.heat_j);
        assert_eq!(report.source.destinations,CONTACT_ROUTES);
        assert_eq!(source.contact.source().resolved().card_identity(),interface.content_hash());
        assert_eq!(electrical.source_evaluations,3*electrical.iterations);
        assert!(report.circuit_trials>electrical.iterations);
        let (t,m,_)=result.into_parts(); thermal=t; motion=m;
    }
    [motion.circuit().charges_c()[0],motion.gap_m(),motion.velocity_m_s(),
        thermal.temperature_k(CellId(30)).unwrap(),thermal.temperature_k(CellId(10)).unwrap()]
}

#[test]
fn g1_joint_contact_temperature_voltage_and_motion_refine_to_independent_five_state_dynamics() {
    let expected=contact_reference(0.03);
    let mut errors=Vec::new(); let mut last=[0.0;5];
    for steps in [4,8,16] {
        last=contact_trajectory(0.03,steps);
        errors.push((0..5).map(|i|(last[i]-expected[i]).abs()).fold(0.0,f64::max));
    }
    assert!(errors[1]<0.65*errors[0] && errors[2]<0.65*errors[1],"{errors:?}");
    assert!(errors[2]<0.005,"bounded first-order thermal/material lag");
    let fixed=contact_trajectory(0.0,16);
    assert!((last[0]-fixed[0]).abs()>1e-3);
    assert!((last[1]-fixed[1]).abs()>1e-7);
    assert!((last[3]-fixed[3]).abs()>1e-3);
}

#[test]
fn g4_contact_motion_heat_refusals_and_retry_preserve_all_three_histories() {
    let (a,ca)=field_material(0.05,2000.0); let (b,cb)=common::material(3.0,2.0,0.0);
    let cards=[a,b]; let curves=[ca,cb]; let thermal=initial(&curves,1.0,None);
    let cells=thermal.cells().clone();
    let axes=[ContactAxis {name:"V",coordinate:ContactCoordinate::VoltageMagnitude}];
    for limit in [0.01,10.0] {
        let interface=scalar::response_card(&cards[0],&cards[1],"V",CONTACT_V,1.2,0.001,0.03,limit);
        let motion=contact_motion(&cards,&interface); let saved=motion.clone();
        let contact_responses=[ContactResponse {id:7,card:&interface,direction:ContactDirection::AlongBranch,axes:&axes}];
        let input=FieldMotionEnthalpyInput {motion:&motion,ports:&[],external_force_n:0.0,destinations:&DESTINATIONS};
        let run=|routes: &[ContactHeatCells],controls,check: &mut dyn FnMut()->bool|
            conductivity::step_with_contact_response_candidate(&thermal,0.1,input,&responses(&cards),
                Some(&contact_responses),routes,controls,conductivity_controls(),&[0.0;2],
                |state|contact_coefficients(state,&motion,&cards,&curves,&interface),
                |_,_,transport|Ok(transport),check);
        if limit<1.0 {assert!(run(&CONTACT_ROUTES,config(),&mut ||false).is_err());}
        else {
            assert!(run(&[],config(),&mut ||false).is_err());
            let mut polls=0; let result=run(&CONTACT_ROUTES,config(),&mut ||{polls+=1;false}).unwrap();
            let mut actual=0;
            assert!(matches!(run(&CONTACT_ROUTES,config(),&mut ||{actual+=1;actual==polls}),
                Err(CircuitEnthalpyError::Cancelled)));
            let limit=fs_couple::circuit_enthalpy::CircuitEnthalpySolve {
                maximum_circuit_trials:result.report().circuit_trials-1,..config()};
            assert!(matches!(run(&CONTACT_ROUTES,limit,&mut ||false),Err(CircuitEnthalpyError::CircuitBudget)));
            let retry=run(&CONTACT_ROUTES,config(),&mut ||false).unwrap();
            assert_eq!(retry.motion(),result.motion());
            assert_eq!(retry.thermal().cells(),result.thermal().cells());
        }
        assert_eq!(motion,saved); assert_eq!(thermal.cells(),&cells); assert_eq!(thermal.time_s(),0.0);
    }
}
